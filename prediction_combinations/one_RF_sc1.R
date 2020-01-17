library(tidyverse)
library(ranger)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Participating team names
submissions <-  readRDS("./submission_data/intermediate_data/sc1_ranked_teams.rds") %>% as.character()
# Median prediction per team per condition (cell line, treatment, time, marker) and the median values of 32 non predicted markers\)
sub_data_median <- readRDS("./submission_data/intermediate_data/sc1_median_conditions_np.rds")
# Error per team per condition (cell line, treatment, time, marker) and the median values of 32 non predicted markers. 
# Mean of error is final scorre
sub_data_err <- readRDS( "./submission_data/intermediate_data/sc1_condErr_medianVals.rds")  %>%
  mutate_if(is.character, as.factor)
# All single cell predictions of all teams and true values (standard) as well as 32 marker values of non predicted markers
sub_data_all <- readRDS("./submission_data/intermediate_data/sc1_all_NP_predictions.rds")  %>%
  mutate_if(is.character, as.factor) %>%
  group_by(treatment, time, marker) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  mutate(n = ifelse(n<500, n, 500)) %>%
  mutate(time = case_when(time==14 ~ 13,
                          time==18 ~ 17,
                          TRUE ~ time))

np_markers <- c("b.CATENIN", "cleavedCas", "CyclinB", "GAPDH", "IdU", "Ki.67", "p.4EBP1", 
                "p.AKT.Thr308.", "p.AMPK", "p.BTK", "p.CREB", "p.FAK", "p.GSK3b", "p.H3", "p.JNK",
                "p.MAP2K3", "p.MAPKAPK2", "p.MEK", "p.MKK3.MKK6", "p.MKK4", "p.NFkB", "p.p38", "p.p53",
                "p.p90RSK", "p.PDPK1", "p.RB", "p.S6K", "p.SMAD23", "p.SRC", "p.STAT1", "p.STAT3", 
                "p.STAT5")

cell_lines <- unique(sub_data_err$cell_line)

# Keep track of scores per CV loop
scores <- tibble("CV_loop" = NA,
                 "val_CL" = NA,
                 "asv_RF_32" = NA)

# Save predicted values 
asv_RF_32_predictions <- tibble()

# Save error per condition of predictions
asv_RF_32_pred_errors <- tibble()
for (i in 1:length(cell_lines)) {
  
  print(paste("Iteration ", i, " out of ", length(cell_lines), ".", sep = ""))

  train_data_all <- sub_data_all %>% 
    filter(!cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
  val_data_all <- sub_data_all %>%
    filter(cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker) 
  
  # Random forest
  # True value as function of all teams and always treatment and time
  RF_formula <- paste0("standard ~", paste(c(submissions, "treatment", "time", np_markers), collapse = "+")) 
  RF_train <- train_data_all %>%
    group_by(treatment, time, marker) %>%
    sample_n(n) %>%
    ungroup()
  RF <- RF_train %>%
    group_by(marker) %>%
    nest() %>%
    mutate(RF = map(data, ~ranger(as.formula(RF_formula), data = ., importance = "impurity", 
                                  max.depth = 6, always.split.variables = c("treatment", "time"))))


  # Predict and score RF
  RF_pred <- val_data_all %>%
    group_by(marker) %>%
    nest() %>%
    left_join(select(RF, -data)) %>%
    mutate(prediction = map2(RF, data, predict))%>%
    mutate(asv_RF_32_pred =  map2(data, prediction, function(x, y) {add_column(x, "asv_RF_32_pred" = y$predictions)})) %>%
    select(marker, asv_RF_32_pred) %>%
    unnest(cols = c(asv_RF_32_pred)) %>%
    ungroup() %>%
    select(-c(submissions, np_markers, n))
  RF_error <- RF_pred %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(asv_RF_32_error = sqrt(sum((standard - asv_RF_32_pred)^2) / n())) %>%
    select(cell_line, treatment, time, marker, asv_RF_32_error) %>%
    ungroup()
  RF_val_score <- RF_error %>%
    pull(asv_RF_32_error) %>%
    mean()
  
  asv_RF_32_predictions <- bind_rows(asv_RF_32_predictions, RF_pred)
  asv_RF_32_pred_errors <- bind_rows(asv_RF_32_pred_errors, RF_error)
  
  # Keep track of scores per CV loop
  scores <- scores %>% add_row("CV_loop" = i,
                   "val_CL" = cell_lines[i],
                   "asv_RF_32" = RF_val_score)
  
}
scores <- filter(scores, !is.na(CV_loop)) 
saveRDS(scores, "./prediction_combinations/SC1/LOO_CV_asvRF_incl32_scores.rds")

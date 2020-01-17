## SC1; DREAM paper
# Do crossvalidation to find how to combine the predictions, leaving two cell lines out as validation set
# A meta decision tree, arbiter,  random forest and linear modelare built
# MDT predicts bet team per condition, with feature Tr, Ti, Ma
# Linear combination SC error predicted by arbiter with features Tr, Ti, Ma, 32 marker values 
# arbiter trained on condition/median level and predicts single cells
# RF directly predicts marker values with features Tr, Ti, Ma and the teams predictions
# A linear model per Tr, Ti, Ma (condition) directly predicts marker values from the teams predictions

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

# To keep count how often features are selected by MDT and the arbiter
feature_importance <- tibble("feature" = c("treatment", "time", "marker", np_markers), 
                             "importance" = 0)
# Keep track of scores per CV loop
scores <- tibble("CV_loop" = NA,
                 "val_CL" = NA,
                 "ttm_MDT" = NA,
                 "cond_lm" = NA,
                 "cond_RF" = NA,
                 "RF" = NA)

# Save predicted values 
all_predictions <- tibble()

# Save error per condition of predictions
all_pred_errors <- tibble()
for (i in 1:length(cell_lines)) {
  
  print(paste("Iteration ", i, " out of ", length(cell_lines), ".", sep = ""))
  
  # select training and validation data
  train_data_err <- sub_data_err %>% # Used for arbiter, MDT and best single
    filter(!cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker) 
  train_data_median <- sub_data_median %>% 
    filter(!cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)  
  train_data_all <- sub_data_all %>% 
    filter(!cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
  train_data_nest <- train_data_all %>% 
    group_by(treatment, time, marker) %>% 
    nest()
  
  val_data_err <- sub_data_err %>%
    filter(cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
  val_data_median <- sub_data_median %>%
    filter(cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
  val_data_all <- sub_data_all %>%
    filter(cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker) 
  val_data_nest <- val_data_all %>% 
    group_by(treatment, time, marker) %>% 
    nest()
  
  ## Train models
  # formula for condition + 32 markers
  features <- paste(feature_importance$feature, collapse = "+")
  
  # Meta decision tree
  MDT <- ranger(best_sub ~ treatment + time + marker, data=train_data_err, num.trees = 100, 
                importance = "impurity")
  
  # True value as function of all teams for lm and RF trained seperately for each condition
  sc_formula <- paste0("standard ~", paste(submissions, collapse = "+"))
  # Linear model
  linModel <- train_data_nest  %>%
    mutate(lm = map(data, ~lm(sc_formula, data = .)))
  # Random forest per condition
  cond_RF_train <- train_data_nest %>%
    mutate(data = map(data, ~sample_n(., n)))
    
  cond_RF <- cond_RF_train %>%
    mutate(RF = map(data, ~ranger(as.formula(sc_formula), data = ., importance = "impurity")))
  
  # Random forest
  # True value as function of all teams, treatment, time and marker
  RF_formula <- paste0("standard ~", paste(c(submissions, "treatment", "time", "marker"), collapse = "+")) 
  RF_train <- train_data_all %>%
    group_by(treatment, time, marker) %>%
    sample_n(n)
  RF <- ranger(as.formula(RF_formula), data = RF_train, importance = "impurity")
  

  
  ## ------------- Assess model performanceson condition level -------------------------------
  # Select for each condition the MDT estimated best predictor
  # Use error of the best predictors over the entire condition to score
  MDT_pred <-  val_data_err %>% 
    add_column(pred = as.character(predict(MDT, val_data_err)$predictions)) %>%
    select(-c(np_markers, best_sub)) %>%
    gather(team, MDT_error, submissions) %>%
    filter(pred==team) %>%
    select(cell_line, treatment, time, marker, MDT_error)
  
  MDT_val_score <- MDT_pred %>%
    pull(MDT_error) %>%
    mean()
  
  ## ------------------ Assess models on single cell level ---------------------------------
  # Linear model
  #Predicted values by linear model
  lm_pred <- val_data_nest %>%
    left_join(select(linModel, -data)) %>%
    filter(treatment == "EGF", time==23) %>%
    mutate(prediction = map2(lm, data, predict)) %>%
    mutate(lm_pred =  map2(data, prediction, function(x, y) {add_column(x, "lm_pred" = y)})) %>%
    select(treatment, time, marker, lm_pred) %>%
    unnest(cols=c(lm_pred)) %>%
    ungroup() %>%
    select(-c(submissions, np_markers, n))
  
  # Error of predicted values by lm  
  lm_pred_error <- lm_pred %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(lm_error = sqrt(sum((standard - lm_pred)^2) / n())) %>%
    select(cell_line, treatment, time, marker, lm_error) %>%
    ungroup()
  
  lm_val_score <- lm_pred_error %>%
    pull(lm_error) %>%
    mean()
  
  # RF trained by condition
  cond_RF_pred <- val_data_nest %>%
    left_join(select(cond_RF, -data)) %>%
    mutate(prediction = map2(RF, data, predict)) %>%
    mutate(cond_RF_pred =  map2(data, prediction, function(x, y) {add_column(x, "cond_RF_pred" = y$predictions)})) %>%
    select(treatment, time, marker, cond_RF_pred) %>%
    unnest(cols = c(cond_RF_pred)) %>%
    ungroup() %>%
    select(-c(submissions, np_markers, n))
  cond_RF_error <- cond_RF_pred %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(cond_RF_error = sqrt(sum((standard - cond_RF_pred)^2) / n())) %>%
    select(cell_line, treatment, time, marker, cond_RF_error) %>%
    ungroup()
  cond_RF_val_score <- cond_RF_error %>%
    pull(cond_RF_error) %>%
    mean()
  
  #Random Forest
  RF_pred <- val_data_all %>%
    add_column(RF_pred = predict(RF, data=val_data_all)$predictions)  %>%
    select(-c(submissions, np_markers, n))
  RF_error <- RF_pred %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(RF_error = sqrt(sum((standard - RF_pred)^2) / n())) %>%
    select(cell_line, treatment, time, marker, RF_error)
  RF_val_score <- RF_error %>%
    pull(RF_error) %>%
    mean()
  
  ## Update score table 
  scores <- scores %>% add_row("CV_loop" = i,
                               "val_CL" = cell_lines[i],
                               "ttm_MDT" = MDT_val_score,
                               "cond_lm" = lm_val_score,
                               "cond_RF" = cond_RF_val_score,
                               "RF" = RF_val_score)
  # Save all predictions
  predictions <- plyr::join_all(list(lm_pred, cond_RF_pred, RF_pred)) %>%
    as_tibble()
  all_predictions <- bind_rows(all_predictions, predictions)
  
  # save all predictions errors
  errors <- plyr::join_all(list(MDT_pred, lm_pred_error, cond_RF_error, RF_error)) %>%
    as_tibble()
  all_pred_errors <- bind_rows(all_pred_errors, errors)
  
  ## Update feature importance of only treatment, timepoint and marker
  MDT_importance <- importance(MDT) %>% 
    enframe(name="feature", value="MDT_imp") %>%
    mutate(MDT_imp = MDT_imp/sum(MDT_imp)) 
  
  RF_importance <- importance(RF) %>% 
    enframe(name="feature", value= "RF_imp") %>%
    mutate(RF_imp = RF_imp/sum(RF_imp)) %>%
    filter(feature %in% feature_importance$feature)
    
  feature_importance <- feature_importance %>%
    left_join(MDT_importance) %>% 
    left_join(RF_importance) %>%
    mutate(importance = importance + MDT_imp + RF_imp) %>%
    select(feature, importance) 
  
}
scores <- filter(scores, !is.na(CV_loop))  
if (TRUE) {
  saveRDS(scores, "./prediction_combinations/SC1/LOO_CV_DREAM_scores.rds")
  saveRDS(feature_importance, "./prediction_combinations/SC1/LOO_CV_DREAM_feature_importance.rds")
  saveRDS(all_predictions, "./prediction_combinations/SC1/LOO_CV_DREAM_all_predictions.rds")
  
}
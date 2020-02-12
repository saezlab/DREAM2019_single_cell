## SC1; ranfom forest
# Do crossvalidation to find how to combine the predictions, leaving one cell line out as validation set
# A meta decision tree, arbiter and linear model are built
# Predictions made by selecting values of best submission as predicted by MDT per condition, 
# values taken from submission with lowest error predicted by the arbiter per condition
# linear combination of all submission based on error predicted by arbiter for single cells 
# and a linear combination of all predicted single cells by linear model that was trained on predicted and true values. 

library(tidyr)
library(purrr)
library(tibble)
library(dplyr)
library(randomForest)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Participating team names
submissions <-  readRDS("./submission_data/intermediate_data/sc1_ranked_teams.rds") %>% as.character()
# Error per team per condition (cell line, treatment, time, marker) and the median values of 32 non predicted markers
sub_data_err <- readRDS( "./submission_data/intermediate_data/sc1_condErr_medianVals.rds") 
# All single cell predictions of all teams and true values (standard) as well as 32 marker values of non predicted markers
sub_data_all <- readRDS("./submission_data/intermediate_data/sc1_all_NP_predictions.rds")
# Format all predictions of teams as CL, Tr, Ti, Ma, Nested predictions of cells
sub_data_nested <- readRDS("./submission_data/intermediate_data/sc1_all_predictions_nested.rds")
  

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
                 "MDT" = NA, 
                 "arbiter" = NA,
                 "lc_arbiter" = NA,
                 "lm" = NA,
                 "single" = NA,
                 "SC_MDT" = NA)
# Save errors of all predictions made by each model
all_pred_error <- tibble()

# All predictions
all_pred <- tibble()

for (i in 1:length(cell_lines)) {
  
  print(paste("Iteration ", i, " out of ", length(cell_lines), ".", sep = ""))
  
  # select training and validation data
  train_data_err <- sub_data_err %>% # Used for arbiter, MDT and best single
    filter(!cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker) %>%
    mutate_if(is.character, as.factor)
  train_data_all <- sub_data_all %>% # Used for lm
    filter(!cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
  
  val_data_err <- sub_data_err %>%
    filter(cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker) %>%
    mutate_if(is.character, as.factor)
  val_data_all <- sub_data_all %>%
    filter(cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker) %>%
    mutate_if(is.character, as.factor)
  
  ## Train models
  # formula for condition + 32 markers
  features <- paste(feature_importance$feature, collapse = "+")
  # Meta Decision Tree
  if (FALSE) {
    MDT_train <- train_data_err %>% 
      group_by(treatment, time, marker, best_sub) %>% 
      summarise(n=n()) %>%
      arrange(treatment, time, marker, n) %>% 
      top_n(1) %>%
      group_by(treatment, time, marker)  %>%
      mutate(n=n()) %>%
      nest(best_sub)  %>%
      mutate(sel_sub = map(data, ~.[sample(1:n, 1),])) %>% 
      select(-c(n, data)) %>% 
      unnest() } # If there is no distinguishing cell lines
  MDT <- randomForest(as.formula(paste0("best_sub ~", features)), data=train_data_err, ntree = 100, 
                      localImp = TRUE)
  
  #Arbiter
  arbiter <- lapply(submissions, function(x)
  {randomForest(as.formula(paste0(x, "~", features)), data=train_data_err, ntree = 100, 
                localImp = TRUE)})
  names(arbiter) <- submissions
  
  # Linear model
  lm_formula <- paste0("standard ~", paste(submissions, collapse = "+"))
  linModel <- lm(lm_formula, data=train_data_all)
  
  # Select best single model
  single_model_scores <- train_data_err  %>%
    select(submissions) %>%
    colMeans()
  single_best <- names(single_model_scores)[which.min(single_model_scores)]
  
  
  ## ------------- Assess model performanceson condition level -------------------------------
  # Select for each condition the MDT estimated best predictor
  # Use error of the best predictors over the entire condition to score
  MDT_pred_best_team <-  val_data_err %>% 
    add_column(pred = as.character(predict(MDT, val_data_err, type = "class"))) %>%
    select(-c(np_markers, best_sub)) %>%
    gather(team, MDT_error, submissions) %>%
    filter(pred==team) 
  
  MDT_pred_values <- MDT_pred_best_team %>%
    left_join(sub_data_nested) %>%
    select(cell_line, treatment, time, marker, data) %>%
    unnest(data) %>%
    arrange(glob_cellID, treatment, time, cellID, fileID, marker) 
  
  MDT_pred_error <- MDT_pred_best_team %>%
    select(cell_line, treatment, time, marker, MDT_error)
  
  MDT_val_score <- MDT_pred_error %>%
    pull(MDT_error) %>%
    mean()
  
  # Predict error of each submission on the validation data with the arbiter
  # Error per condition, not per single cell
  arbiter_pred <- lapply(arbiter, FUN = function(a) {predict(a, val_data_err)}) %>% as_tibble()
  
  # Select submission with lowest predicted error and score the combination
  arbiter_pred_best_team <- val_data_err %>%
    add_column(pred = names(arbiter_pred)[max.col(-arbiter_pred)]) %>%
    select(-c(np_markers, best_sub)) %>%
    gather(team, arbiter_error, submissions) %>%
    filter(pred==team)
  
  arbiter_pred_values <- arbiter_pred_best_team  %>%
    left_join(sub_data_nested) %>%
    select(cell_line, treatment, time, marker, data) %>%
    unnest(data) %>%
    arrange(glob_cellID, treatment, time, cellID, fileID, marker) 
  
  arbiter_pred_error <- arbiter_pred_best_team %>%
    select(cell_line, treatment, time, marker, arbiter_error)
  
  arbiter_val_score <- arbiter_pred_error  %>%
    pull(arbiter_error) %>%
    mean()
  
  ## ------------------ Assess models on single cell level ---------------------------------
  # Linear combination of prediction based on predicted error by the arbiter W is the weights
  SC_arbiter_pred <- lapply(arbiter, FUN = function(a) {predict(a, val_data_all)}) %>% as_tibble()
  W <- lapply(SC_arbiter_pred, function(v) {(1/rowSums(SC_arbiter_pred)/v)}) %>% as_tibble()
  W <- lapply(W, function(x){x/rowSums(W)}) %>% as_tibble()
  
  lc_arbiter_pred <- val_data_all %>% 
    add_column(prediction = rowSums(select(val_data_all, submissions) * W)) %>%
    arrange(glob_cellID, treatment, time, cellID, fileID, marker) 
    
  lc_arbiter_error <- lc_arbiter_pred %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(lc_arbiter_error = sqrt(sum((standard - prediction)^2) / n())) %>%
    ungroup() %>%
    select(cell_line, treatment, time, marker, lc_arbiter_error)
  
  lc_arbiter_val_score <- lc_arbiter_error %>%
    pull(lc_arbiter_error) %>%
    mean()
  
  # Linear model
  lm_pred <- val_data_all %>% 
    add_column(prediction = predict(linModel, val_data_all)) %>%
    arrange(glob_cellID, treatment, time, cellID, fileID, marker) 
  
  lm_error <- lm_pred %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(lm_error = sqrt(sum((standard - prediction)^2) / n())) %>%
    select(cell_line, treatment, time, marker, lm_error)
  
  lm_val_score <- lm_error %>%
    pull(lm_error) %>%
    mean()
  
  # Score predictions of single best model
  SB_pred <- sub_data_nested %>%
    ungroup() %>%
    filter(team == single_best) %>%
    right_join(val_data_err) %>%
    select(cell_line, treatment, time, marker, data) %>%
    add_column(model = "single best") %>%
    unnest(data) %>%
    arrange(model, glob_cellID, treatment, time, cellID, fileID, marker) 
  
  SB_error <- val_data_err %>% 
    select(cell_line, treatment, time, marker, single_best) %>%
    rename(SB_error = single_best)
  
  SB_val_score <- SB_error %>%
    pull(SB_error) %>%
    mean()
  
  # MDT on single cells
  SC_MDT_pred <-  val_data_all %>% 
    add_column(pred = predict(MDT, val_data_all, type = "class")) %>%
    select(-np_markers) %>%
    gather(team, prediction, submissions) %>%
    filter(pred==team) %>%
    arrange(glob_cellID, treatment, time, cellID, fileID, marker) 
  
  SC_MDT_error <-SC_MDT_pred %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(SC_MDT_error = sqrt(sum((standard - prediction)^2) / n())) %>%
    select(cell_line, treatment, time, marker, SC_MDT_error)
  
  SC_MDT_val_score <- SC_MDT_error %>%
    pull(SC_MDT_error) %>%
    mean()
  
  ## Update score table 
  scores <- scores %>% add_row("CV_loop" = i,
                               "val_CL" = cell_lines[i],
                               "MDT" = MDT_val_score, 
                               "arbiter" = arbiter_val_score,
                               "lc_arbiter" = lc_arbiter_val_score,
                               "lm" = lm_val_score,
                               "single" = SB_val_score,
                               "SC_MDT" = SC_MDT_val_score)
  
  ## Update feature importance
  MDT_importance <- importance(MDT)[,"MeanDecreaseGini"] %>% 
    enframe(name="feature", value="MDT_imp") %>%
    mutate(MDT_imp = MDT_imp/sum(MDT_imp))

  arbiter_importance <- lapply(arbiter, function(x){importance(x)[,"IncNodePurity"]}) %>%
    bind_cols() %>%
    mutate_at(submissions, ~./sum(.)) %>%
    rowSums() 
  
  feature_importance <- feature_importance %>%
    left_join(MDT_importance) %>% 
    add_column(arbiter_imp = arbiter_importance) %>%
    mutate(importance = importance + MDT_imp + arbiter_imp) %>%
    select(feature, importance) 
  
  # Save all predictions and all errors of predictions
  predictions <- val_data_all %>%
    select(glob_cellID, treatment, time, cellID, fileID, marker) %>%
    arrange(glob_cellID, treatment, time, cellID, fileID, marker) %>%
    add_column("MDT" = MDT_pred_values$prediction,
               "arbiter" = arbiter_pred_values$prediction,
               "lc_arbiter" = lc_arbiter_pred$prediction,
               "lm" = lm_pred$prediction,
               "SB" = SB_pred$prediction,
               "SC_MDT" = SC_MDT_pred$prediction)
  all_pred <- bind_rows(all_pred, predictions)
  
  errors <- plyr::join_all(list(MDT_pred_error, arbiter_pred_error, lc_arbiter_error, lm_error, 
                                SB_error, SC_MDT_error)) %>%
    as_tibble()
  
  all_pred_error <- bind_rows(all_pred_error, errors)
}
scores <- filter(scores, !is.na(CV_loop))  
if (FALSE) {
  saveRDS(scores, "./prediction_combinations/SC1/LOO_CV_RF_scores.rds")
  saveRDS(feature_importance, "./prediction_combinations/SC1/LOO_CV_RF_feature_importance.rds")
  saveRDS(all_pred_error, "./prediction_combinations/SC1/LOO_CV_RF_pred_error.rds")
  saveRDS(all_pred, "./prediction_combinations/SC1/LOO_CV_RF_all_predictions.rds")
  
}

CV_scores <- readRDS("./prediction_combinations/SC1/LOO_CV_RF_scores.rds")
CV_features <- readRDS("./prediction_combinations/SC1/LOO_CV_RF_feature_importance.rds")

colMeans(select(CV_scores, -c(CV_loop, val_CL)))

CV_scores %>%
  gather(model, score, -c(CV_loop, val_CL)) %>%
  ggplot(aes(model, score, fill=model)) +
  geom_boxplot() +
  labs(title = "SC1; leave one out CV") +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

CV_scores %>%
  gather(model, score, -c(CV_loop, val_CL)) %>%
  ggplot(aes(val_CL, score, fill = model)) +
  geom_col(position = "dodge") +
  labs(title = "SC1; leave one out CV") +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

CV_features  %>%
  arrange(desc(importance)) %>%
  ggplot(aes(reorder(feature, -importance), importance/(6*23), fill=feature)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90), 
        axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.position = "none") +
  labs(y="feature importance", x="feature", title="SC1; leave one out CV")

# Check why linear model performs so poorly
if (FALSE) {
teams <- sample(submissions, 10)

score <- sub_data_all %>%
  select(glob_cellID, cell_line, treatment, time, marker, cellID, fileID, standard, teams) %>%
  gather(team, value, teams) %>%
  group_by(glob_cellID, cell_line, treatment, time, marker, cellID, fileID, standard) %>%
  summarise(prediction = mean(value)) 

condition_error <- score %>%
  group_by(cell_line, treatment, time, marker) %>%
  summarise(error = sqrt(sum((standard - prediction)^2) / n()))
  
condition_error %>% group_by(cell_line) %>% summarise(score = mean(error))%>% pull(score) %>% mean()

lm <- lm(lm_formula, data=sub_data_all)
predictions <- predict(lm, sub_data_all)

overfitted_cond_error <- sub_data_all %>%
  select(glob_cellID, cell_line, treatment, time, marker, cellID, fileID, standard) %>%
  add_column(prediction = predictions)  %>%
  group_by(cell_line, treatment, time, marker) %>%
  summarise(error = sqrt(sum((standard - prediction)^2) / n()))
overfitted_cond_error %>% group_by(cell_line) %>% summarise(score = mean(error)) %>% pull(score) %>% mean()

sub_data_err %>%
  group_by(cell_line) %>%
  summarise_at(submissions, mean) %>%
  gather(team, CL_score, submissions) %>%
  group_by(cell_line) %>%
  top_n(-5, CL_score) %>%
  ggplot(aes(cell_line, CL_score, colour = team)) +
  geom_point()
}





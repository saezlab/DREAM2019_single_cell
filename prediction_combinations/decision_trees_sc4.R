## SC4
# Do crossvalidation to find how to combine the predictions, leaving two cell lines out as validation set
# A meta decision tree, arbiter and linear model are built
# Predictions made by selecting value of best submission as predicted by MDT, value taken from submission with lowest
# error predicted by the arbiter, linear combination of all submission basedon error predicted by arbiter and a l
# linear combination of of all submission by linear model that was trained on predicted and true values. 

library(tidyverse)
library(rpart)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

submissions <-  readRDS("./submission_data/intermediate_data/sc4_ranked_teams.rds") %>%
  as.character()
sub_data_values <- readRDS("./submission_data/intermediate_data/sc4_median_and_values.rds")
sub_data_err <- readRDS("./submission_data/intermediate_data/sc4_median_and_error.rds")
cell_lines <- unique(sub_data_err$cell_line)

np_markers <- sub_data_values %>% 
  select(-c(submissions, cell_line, treatment, time, marker, standard)) %>% 
  colnames()

if (FALSE){
# Each row represents a unique combination of two cell lines
# Each row will be used as validation cell lines once
CL_combi <- expand.grid(cell_lines, cell_lines) %>% 
  as_tibble() %>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  filter(! Var1 == Var2) %>%
  group_by(grp = paste(pmax(Var1, Var2), pmin(Var1, Var2), sep = "_")) %>%
  slice(1) %>%
  ungroup() %>%
  select(-grp)
}

# To keep count how often features are selected by MDT and the arbiter
feature_count <- tibble("feature" = c("treatment", "time", np_markers), 
                        "count" = 0)
# Keep track of scores per CV loop
scores <- tibble("CV_loop" = NA,
                 "val_CL" = NA, 
                 "MDT" = NA, 
                 "arbiter" = NA,
                 "lc_arbiter" = NA,
                 "lm" = NA,
                 "single" = NA)
# Save predictions made by each model
all_predictions <- tibble()
for (i in 1:length(cell_lines)) {
  print(paste("Iteration ", i, " out of ", length(cell_lines), ".", sep = ""))
  
  # select training and validation data
  train_data_err <- sub_data_err %>%
    filter(!cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
  train_data_value <- sub_data_values %>%
    filter(!cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
  
  val_data_err <- sub_data_err %>%
    filter(cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
  val_data_value <- sub_data_values %>%
    filter(cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
  
  ## Train models
  # formula for condition + 32 markers
  features <- paste(feature_count$feature, collapse = "+")
  
  # Meta Decision Tree
  
  if (FALSE) {
  MDT_train <- train_data_err %>% #Select 1 best submission per condition
    group_by(treatment, time, marker, best_sub) %>% 
    summarise(n=n()) %>%
    arrange(treatment, time, marker, n) %>% 
    top_n(1) %>%
    group_by(treatment, time, marker)  %>%
    mutate(n=n()) %>%
    nest(best_sub)  %>%
    mutate(sel_sub = map(data, ~.[sample(1:n, 1),])) %>% 
    select(-c(n, data)) %>% 
    unnest() %>% 
    ungroup()}  # If there is no distinguishing cell lines
  MDT <- rpart(paste0("best_sub ~", features), data=train_data_err, method="class") 
  
  #Arbiter
  arbiter <- lapply(submissions, function(x)
  {rpart(paste0(x, "~", features), data=train_data_err, method="anova")})
  names(arbiter) <- submissions
  # Linear model
  lm_formula <- paste0("standard ~", paste(submissions, collapse = "+"))
  linModel <- lm(lm_formula, data=train_data_value)
  
  # Select best single model
  single_model_scores <- train_data_err  %>%
    select(submissions) %>%
    colMeans()
  single_best <- names(single_model_scores)[which.min(single_model_scores)]
  
  ## Assess model performances
  # Select for each condition the value of the MDT estimated best predictor and score the submission
  MDT_pred <-  val_data_value %>% 
    add_column(pred = predict(MDT, val_data_err, type = "class")) %>%
    gather(team, prediction, submissions) %>%
    filter(pred==team) %>%
    select(cell_line, treatment, time, marker, standard, prediction) %>%
    rename(MDT_pred = prediction)
  
  MDT_val_score <- MDT_pred %>%
    group_by(cell_line, treatment, marker) %>%
    summarise(RMSE = sqrt(sum((standard - MDT_pred)^2) / n())) %>%
    pull(RMSE) %>%
    mean()
  
  # Predict error of each submission on the validation data with the arbiter
  arbiter_pred <- lapply(arbiter, FUN = function(a) {predict(a, val_data_err)})   %>% as_tibble()
  
  # Select submission with lowest predicted error and score the combination
  arbiter_pred_value <- val_data_value %>%
    add_column(pred = names(arbiter_pred)[max.col(-arbiter_pred)]) %>%
    gather(team, prediction, submissions) %>%
    filter(pred==team)  %>%
    select(cell_line, treatment, time, marker, standard, prediction) %>%
    rename(arbiter_pred = prediction)
  
  arbiter_val_score <- arbiter_pred_value %>%
    group_by(cell_line, treatment, marker) %>%
    summarise(RMSE = sqrt(sum((standard - arbiter_pred)^2) / n())) %>%
    pull(RMSE) %>%
    mean()
  
  # Linear combination of prediction based on predicted error by the arbiter W is the weights
  W <- lapply(arbiter_pred, function(v) {(1/rowSums(arbiter_pred)/v)}) %>% as_tibble()
  W <- lapply(W, function(x){x/rowSums(W)}) %>% as_tibble()
  lc_arbiter_pred <- val_data_value %>% 
    add_column(lc_arbiter_pred = rowSums(select(val_data_value, submissions) * W)) %>%
    select(cell_line, treatment, time, marker, standard, lc_arbiter_pred)
  
  lc_arbiter_val_score <- lc_arbiter_pred %>%
    group_by(cell_line, treatment, marker) %>%
    summarise(RMSE = sqrt(sum((standard - lc_arbiter_pred)^2) / n())) %>%
    pull(RMSE) %>%
    mean()
  
  # Linear model
  lm_pred <- val_data_value %>% 
    add_column(lm_pred = predict(linModel, val_data_value)) %>%
    select(cell_line, treatment, time, marker, standard, lm_pred)
    
  lm_val_score <- lm_pred %>%
    group_by(cell_line, treatment, marker) %>%
    summarise(RMSE = sqrt(sum((standard - lm_pred)^2) / n())) %>%
    pull(RMSE) %>%
    mean()
  
  # Score predictions of single best model
  SB_pred <-val_data_value  %>%
    select(cell_line, treatment, time, marker, standard, single_best) %>%
    rename(SB_pred  = single_best)
  SB_val_score <- SB_pred %>%
    group_by(cell_line, treatment, marker) %>%
    summarise(RMSE = sqrt(sum((standard - SB_pred)^2) / n())) %>%
    pull(RMSE) %>%
    mean()
  
  ## Update score table
  scores <- scores %>% add_row("CV_loop" = i,
                               "val_CL" = cell_lines[i], 
                               "MDT" = MDT_val_score, 
                               "arbiter" = arbiter_val_score,
                               "lc_arbiter" = lc_arbiter_val_score,
                               "lm" = lm_val_score,
                               "single" = SB_val_score)
  
  ## Count selected features by MDT and arbiter
  MDT_features <- unique(sub("[<>=].*", "", labels(MDT))) %>% 
    as_tibble() %>% 
    count(value) %>% 
    rename(count_MDT=n, feature = value)
  arbiter_features <- lapply(arbiter, function(x){unique(sub("[<>=].*", "", labels(x)))}) %>% 
    unlist() %>% as_tibble() %>% count(value) %>% rename(feature = value)
  feature_count <- feature_count %>%
    left_join(MDT_features) %>%
    left_join(arbiter_features) %>%
    mutate(count_MDT = ifelse(is.na(count_MDT), 0, count_MDT),
           n=  ifelse(is.na(n), 0, n),
           count = count+count_MDT+n) %>%
    select(feature, count)
  
  # Save all predictions
  predictions <- plyr::join_all(list(MDT_pred, arbiter_pred_value, lc_arbiter_pred, lm_pred, SB_pred)) %>%
    as_tibble()
  
  all_predictions <- bind_rows(all_predictions, predictions)
}

scores <- filter(scores, !is.na(CV_loop))

if (FALSE) {
  saveRDS(scores, "./prediction_combinations/SC4/LOO_CV_scores.rds")
  saveRDS(feature_count, "./prediction_combinations/SC4/LOO_CV_feature_counts.rds")
  saveRDS(all_predictions, "./prediction_combinations/SC4/LOO_CV_all_predictions.rds")
  
}

CV_scores <- readRDS("./prediction_combinations/SC4/LOO_CV_scores.rds")
CV_features <- readRDS("./prediction_combinations/SC4/LOO_CV_feature_counts.rds")

colMeans(select(CV_scores, -c(CV_loop, val_CL)))

CV_scores %>%
  gather(model, score, -c(CV_loop, val_CL)) %>%
  ggplot(aes(model, score, fill=model)) +
  geom_boxplot() +
  labs(title = "SC4; leave one out CV") +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

CV_scores %>%
  gather(model, score, -c(CV_loop, val_CL)) %>%
  ggplot(aes(val_CL, score, fill = model)) +
  geom_col(position = "dodge") +
  labs(title = "SC4; leave one out CV") +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

CV_features  %>%
  arrange(desc(count)) %>%
  ggplot(aes(reorder(feature, -count), count/(2*21), fill=feature)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90), 
        axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.position = "none") +
  labs(y="proprtion selected", x="feature", title="SC4; leave one out CV")

## SC2
# Do crossvalidation to find how to combine the predictions, leaving two cell lines out as validation set
# A meta decision tree, arbiter and linear model are built
# Predictions made by selecting value of best submission as predicted by MDT, value taken from submission with lowest
# error predicted by the arbiter, linear combination of all submission basedon error predicted by arbiter and a l
# linear combination of of all submission by linear model that was trained on predicted and true values. 

library(tidyverse)
library(rpart)
library(ranger)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")
source("./scoring_scripts/score_sc2_noFile.R")

# Sample 10.000 cells from multivarite distribution 
sample_MVN <- function(means, covs){
  sample <- MASS::mvrnorm(10000, mu = means, Sigma =  covs, tol=0.005) %>%
    as_tibble()
  return(sample)
}
# From lower triangle matrix to full matrix
makeSymm <- function(x) {
  m <- x %>% separate(variable, into = c("V1", "V2"), sep = "_") %>%
    pivot_wider(names_from = V2, values_from = prediction) %>%
    column_to_rownames("V1") 
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}


submissions <-  readRDS("./submission_data/intermediate_data/sc2_ranked_teams.rds") %>%
  as.character()
# Values are sum of squares error between validaion data and prediction
# to get the score sum up all values of 1 prediction and then take the square root of the sum
sub_data_err <- readRDS("./submission_data/intermediate_data/sc2_err_median_EGF0.rds")
# Statistics per condition of the standard and calculated from the predictions, including median EGF t=0 values per condition
sub_data_values <- readRDS("./submission_data/intermediate_data/sc2_values_median_EGF0.rds") 


# All actual predictions in nested data frame, time 16 already changed to 17
nested_predictions <- readRDS("./submission_data/intermediate_data/sc2_all_predictions_nested.rds")

cell_lines <- unique(sub_data_values$cell_line) %>% 
  sample(12, replace = FALSE) # randomise order

np_markers <- colnames(sub_data_err)[!colnames(sub_data_err) %in% c("cell_line", "treatment",
                                                                    "time", "best_sub", 
                                                                    submissions)]
# To keep count how often features are selected by MDT and the arbiter
feature_count <- tibble("feature" = c("treatment", "time", np_markers), 
                        "count" = 0)
# Keep track of scores per CV loop
scores <- tibble("CV_loop" = NA, 
                 "val_CL" =NA,
                 "MDT" = NA,
                 "arbiter" = NA,
                 "lc_arbiter" = NA,
                 "lm" = NA,
                 "stats_RF" = NA,
                 "stats_lc_arbiter" = NA,
                 "single" = NA)
# Save predictions and SSQ errors made by each model
all_predictions <- tibble()
all_SSQ <- tibble()

# Every loop two cell lines are used as validation, each cell line occurs as validation once
for (i in seq(1, by=2, len= 0.5*length(cell_lines))) {
  print(paste("Iteration", ceiling(i/2), sep = " "))
  
  # select training and validation data
  train_data_err <- sub_data_err %>%
    filter(!cell_line %in% cell_lines[i:(i+1)]) %>%
    arrange(cell_line, treatment, time)
  train_data_value <- sub_data_values %>%
    filter(!cell_line %in% cell_lines[i:(i+1)]) %>%
    arrange(cell_line, treatment, time)
  
  val_data_err <- sub_data_err %>%
    filter(cell_line %in% cell_lines[i:(i+1)]) %>%
    arrange(cell_line, treatment, time)
  val_data_value <- sub_data_values %>%
    filter(cell_line %in% cell_lines[i:(i+1)]) %>%
    arrange(cell_line, treatment, time)
  
  ## Train models
  # Meta Decision Tree
  features <- paste(feature_count$feature, collapse = "+")
  MDT <- rpart(paste0("best_sub ~", features), data=train_data_err, method="class",
               control = rpart.control(minsplit=10))
  
  #Arbiter
  arbiter <- lapply(submissions, function(x)
  {rpart(paste0(x, "~", features), data=train_data_err, method="anova", 
         control = rpart.control(minsplit=10))})
  names(arbiter) <- submissions
  
  # Linear model
  lm_formula <- paste0("standard ~", paste(submissions, collapse = "+"))
  linModel <- lm(lm_formula, data=train_data_value)
  
  # Random forest trained on statistics of the predictions
  stats_RF <- ranger(lm_formula, data=train_data_value, importance = "impurity")
  
  # Select best single model
  single_model_scores <- train_data_err %>%
    select(submissions) %>%
    colSums() %>%
    sqrt()
  single_best <- names(single_model_scores)[which.min(single_model_scores)]
  
  
  ## Assess model performances
  # Select for each condition the value of the MDT estimated best predictor and score the submission
  MDT_pred_best_team <-  val_data_err %>% 
    add_column(pred = predict(MDT, val_data_err, type = "class")) %>%
    gather(team, MDT_SSQ, submissions) %>%
    filter(pred==team)
 
  MDT_pred <- MDT_pred_best_team %>%
    select(cell_line, treatment, time, team) %>%
    left_join(nested_predictions) %>%
    add_column(model="MDT", .before = "cell_line") %>%
    select(-team) %>%
    rename(prediction = data)
  
  MDT_SSQ <- MDT_pred_best_team %>%
    select(cell_line, treatment, time, MDT_SSQ)
    
  MDT_val_score <-MDT_SSQ %>%
    select(MDT_SSQ) %>%
    colSums() %>%
    sqrt()
  
  # Predict error of each submission on the validation data with the arbiter
  arbiter_pred <- lapply(arbiter, FUN = function(a) {predict(a, val_data_err)})   %>% as_tibble()
  
  # Select submission with lowest predicted error and score the combination
  arbiter_pred_best_team <- val_data_err %>%
    add_column(pred = names(arbiter_pred)[max.col(-arbiter_pred)]) %>%
    gather(team, arbiter_SSQ, submissions) %>%
    filter(pred==team) 
  
  arbiter_pred_cells <- arbiter_pred_best_team %>%
    select(cell_line, treatment, time, team) %>%
    left_join(nested_predictions) %>%
    add_column(model = "arbiter", .before =  "cell_line") %>%
    select(-team) %>%
    rename(prediction = data)
  
  arbiter_SSQ <- arbiter_pred_best_team %>%
    select(cell_line, treatment, time, arbiter_SSQ)
  
  arbiter_val_score <- arbiter_SSQ %>%
    select(arbiter_SSQ) %>%
    colSums() %>%
    sqrt()
  
  # Linear combination of prediction based on predicted error by the arbiter W is the weights
  W <- lapply(arbiter_pred, function(v) {(1/rowSums(arbiter_pred)/v)}) %>% as_tibble()
  W <- lapply(W, function(x){x/rowSums(W)}) %>% as_tibble()
  lc_arbiter_pred <- val_data_err %>% 
    select(-c(submissions, best_sub, np_markers)) %>%
    bind_cols(W) %>%
    mutate(time=as.numeric(time)) %>%
    gather(team, weight, submissions) %>%
    left_join(nested_predictions) %>%
    mutate(sample = map2(data, weight, sample_frac)) %>%
    select(cell_line, treatment, time, sample) %>%
    add_column(model  = "lc_arbiter", .before = "cell_line") %>%
    unnest(sample) %>%
    group_by(model, cell_line, treatment, time) %>%
    nest(prediction = -group_cols()) %>%
    ungroup()
  
  lc_arbiter_SSQ <- lc_arbiter_pred  %>% 
    select(-model) %>%
    unnest(prediction) %>%
    data_to_stats() %>%
    rename(pred_stat_value = stat_value) %>%
    arrange(cell_line, treatment, time, stat_variable) %>%
    bind_cols(select(val_data_value, standard)) %>%
    select(cell_line, treatment, time, stat_variable, standard, pred_stat_value) %>%
    group_by(cell_line, treatment, time) %>%
    summarise(lc_arbiter_SSQ = sum((standard - pred_stat_value)^2)) %>%
    ungroup()
  
  lc_arbiter_val_score <- lc_arbiter_SSQ %>%
    select(lc_arbiter_SSQ) %>%
    colSums() %>%
    sqrt()
  
  # Combine predicted error of arbiter on statistics level, make MVN and sample 10,000 cells
  val_data_long <- val_data_value %>%
    select(cell_line, treatment, time, stat_variable, submissions) %>%
    gather(team, pred_stat, submissions) 
  
  stats_lc_arbiter_pred <- val_data_err %>% 
    select(-c(submissions, best_sub, np_markers)) %>%
    bind_cols(W) %>%
    gather(team, weight, submissions) %>%
    right_join(val_data_long) %>%
    group_by(cell_line, treatment, time, stat_variable) %>%
    summarise(prediction = sum(weight*pred_stat))  %>%
    separate(stat_variable, into = c("stat", "variable"),sep="-") %>%
    group_by(cell_line, treatment, time, stat) %>%
    nest(stat_values = c(variable, prediction)) %>%
    pivot_wider(names_from = stat, values_from = stat_values) %>%
    mutate(full_cov = map(cov, makeSymm),
           full_mean = map(mean, ~pull(., prediction))) %>%
    mutate(prediction = map2(full_mean, full_cov, sample_MVN)) %>%
    select(cell_line, treatment, time, prediction)  %>%
    add_column(model="stats_lc_arbiter", .before = "cell_line")
  
  stats_lc_arbiter_SSQ <- stats_lc_arbiter_pred %>%
    select(-model) %>%
    unnest(prediction) %>%
    data_to_stats() %>%
    rename(pred_stat_value = stat_value) %>%
    arrange(cell_line, treatment, time, stat_variable) %>%
    bind_cols(select(val_data_value, standard)) %>%
    select(cell_line, treatment, time, stat_variable, standard, pred_stat_value) %>%
    group_by(cell_line, treatment, time) %>%
    summarise(stats_arbiter_SSQ = sum((standard - pred_stat_value)^2)) %>%
    ungroup()
  
  stats_lc_arbiter_score <- stats_lc_arbiter_SSQ %>%
    select(stats_arbiter_SSQ) %>%
    colSums() %>%
    sqrt()
  
  # Linear model
  lm_pred <- val_data_value %>% 
    add_column(prediction = predict(linModel, val_data_value)) %>%
    select(cell_line, treatment, time, stat_variable, prediction) %>%
    separate(stat_variable, into = c("stat", "variable"),sep="-") %>%
    group_by(cell_line, treatment, time, stat) %>%
    nest(stat_values = c(variable, prediction)) %>%
    pivot_wider(names_from = stat, values_from = stat_values) %>%
    mutate(full_cov = map(cov, makeSymm),
           full_mean = map(mean, ~pull(., prediction))) %>%
    mutate(prediction = map2(full_mean, full_cov, sample_MVN)) %>%
    select(cell_line, treatment, time, prediction) %>%
    add_column(model="lm", .before = "cell_line")
    
  lm_SSQ <- lm_pred %>%
    select(-model) %>%
    unnest(prediction) %>%
    data_to_stats() %>%
    rename(pred_stat_value = stat_value) %>%
    arrange(cell_line, treatment, time, stat_variable) %>%
    bind_cols(select(val_data_value, standard)) %>%
    select(cell_line, treatment, time, stat_variable, standard, pred_stat_value) %>%
    group_by(cell_line, treatment, time) %>%
    summarise(lm_SSQ = sum((standard - pred_stat_value)^2)) %>%
    ungroup() 
  
  lm_val_score <- lm_SSQ %>%
    select(lm_SSQ) %>%
    colSums() %>%
    sqrt()
  
  # Random Forest: Combine the submissions on statistical level, create MVN and sample 10,000 cells w
  stats_RF_pred <- val_data_value %>% 
    add_column(prediction = predict(stats_RF, val_data_value)$predictions) %>%
    select(cell_line, treatment, time, stat_variable, prediction) %>%
    separate(stat_variable, into = c("stat", "variable"),sep="-") %>%
    group_by(cell_line, treatment, time, stat) %>%
    nest(stat_values = c(variable, prediction)) %>%
    pivot_wider(names_from = stat, values_from = stat_values) %>%
    mutate(full_cov = map(cov, makeSymm),
           full_mean = map(mean, ~pull(., prediction))) %>%
    mutate(prediction = map2(full_mean, full_cov, sample_MVN)) %>%
    select(cell_line, treatment, time, prediction) %>%
    add_column(model="stats_RF", .before = "cell_line")
  
  stats_RF_SSQ <- stats_RF_pred %>%
    select(-model) %>%
    unnest(prediction) %>%
    data_to_stats() %>%
    rename(pred_stat_value = stat_value) %>%
    arrange(cell_line, treatment, time, stat_variable) %>%
    bind_cols(select(val_data_value, standard)) %>%
    select(cell_line, treatment, time, stat_variable, standard, pred_stat_value) %>%
    group_by(cell_line, treatment, time) %>%
    summarise(stats_RF_SSQ = sum((standard - pred_stat_value)^2)) %>%
    ungroup()
  
  stats_RF_score <- stats_RF_SSQ %>%
    select(stats_RF_SSQ) %>%
    colSums() %>%
    sqrt()
  
  # Score predictions of single best model 
  SB_pred <- nested_predictions %>%
    filter(team == single_best) %>%
    right_join(select(val_data_err, cell_line, treatment, time)) %>%
    add_column(model = "single best", .before = "team") %>%
    select(-team) %>%
    rename(prediction = data)
  
  SB_SSQ <- val_data_err %>%
    select(cell_line, treatment, time, single_best) %>%
    rename(SB_SSQ = single_best) 
  
  SB_val_score <- SB_SSQ %>%
    select(SB_SSQ) %>%
    colSums() %>%
    sqrt()
  
  ## Update score table 
  scores <- scores %>% add_row("CV_loop" = i, 
                               "val_CL" = paste(cell_lines[i:(i+1)], collapse = "_"),
                               "MDT" = MDT_val_score, 
                               "arbiter" = arbiter_val_score,
                               "lc_arbiter" = lc_arbiter_val_score,
                               "lm" = lm_val_score,
                               "stats_RF" = stats_RF_score,
                               "stats_lc_arbiter" = stats_lc_arbiter_score,
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
           n =  ifelse(is.na(n), 0, n),
           count = count+count_MDT+n) %>%
    select(feature, count)
  
  # Save all predictions
  all_predictions <- bind_rows(all_predictions, list(MDT_pred, arbiter_pred_cells, lc_arbiter_pred, stats_lc_arbiter_pred,
                                                     lm_pred, stats_RF_pred ,SB_pred))
  
  SSQs <- plyr::join_all(list(MDT_SSQ, arbiter_SSQ, lc_arbiter_SSQ, stats_lc_arbiter_SSQ,
                              lm_SSQ, stats_RF_SSQ ,SB_SSQ)) %>%
    as_tibble()
  all_SSQ <- bind_rows(all_SSQ, SSQs)
}

scores <- filter(scores, !is.na(CV_loop))  

if (FALSE) {
  saveRDS(scores, "./prediction_combinations/SC2/L2O_CV_scores.rds")
  saveRDS(feature_count, "./prediction_combinations/SC2/L2O_CV_feature_counts.rds")
  saveRDS(all_predictions, "./prediction_combinations/SC2/L2O_CV_all_predictions.rds")
  saveRDS(all_SSQ, "./prediction_combinations/SC2/L2O_CV_all_SSQ.rds")
  
}

CV_scores <- readRDS("./prediction_combinations/SC2/L2O_CV_scores.rds")
CV_features <- readRDS("./prediction_combinations/SC2/L2O_CV_feature_counts.rds")

colMeans(select(CV_scores, -c(CV_loop, val_CL)))

CV_scores %>%
  gather(model, score, -c(CV_loop, val_CL)) %>%
  ggplot(aes(model, score, fill=model)) +
  geom_boxplot() +
  labs(title = "SC2; leave two out CV") +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

CV_scores %>%
  gather(model, score, -c(CV_loop, val_CL)) %>%
  ggplot(aes(val_CL, score, fill = model)) +
  geom_col(position = "dodge") +
  labs(title = "SC2; leave two out CV") +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

CV_features  %>%
  arrange(desc(count)) %>%
  ggplot(aes(reorder(feature, -count), count/(6*17), fill=feature)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90), 
        axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.position = "none") +
  labs(y="proprtion selected", x="feature", title="SC2; leave two out CV")


# Cross validated combination of models for SC4

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell") 

source("./scoring_scripts/score_sc1_noFile.R")

required_columns <- c("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID", 
                      "p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")
reporters <- c("p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")

# Scoring function wtth input a tibble, not a csv file. Prediction data already in long format
score_sc1_long_format <- function(prediction_data_long, validation_data) {
  
  validation_data_long <- validation_data %>% gather(key = "marker", value = "test", reporters)
  
  combined_data <- validation_data_long %>%
    full_join(prediction_data_long, by = c("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID", "marker"))
  
  ### Calculate score --------------------
  # calculate the RMSE for each condition
  RMSE_cond <- combined_data %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(RMSE = sqrt(sum((test - prediction)^2) / n()))
  final_score <- mean(RMSE_cond$RMSE)
  return(final_score)
}

submission_folder <- "./submission_data/final/SC1"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc1.csv")) 

validation_data <- read_csv("~/Seafile/validation_data/sc1gold.csv") %>% select(required_columns)
cell_lines <- validation_data %>% pull(cell_line) %>% unique

submissions <- list.files(submission_folder)
submissions <- submissions[!startsWith(submissions, "leaderboard")]

# For setting up and checking workflow select only three submissions
n_samples <- 3
selected_submissions <- sample(submissions, n_samples, replace = FALSE)

# Read prediction files
predictions <- bind_rows(lapply(submissions, function(x) {
  read_csv(file.path(submission_folder,x))  %>%
    select(required_columns) %>% 
    add_column(subID = as.numeric(sub(".csv", "", x)))
}))

CV_scores <- tibble("fold" = NA, 
                    "mean" = NA, 
                    "median" = NA, 
                    "best_single" = NA, 
                    "weighted" = NA, 
                    "lm" = NA)

for (c in 1:length(cell_lines)) {
  print(paste("CV fold ", c, "cell line ", cell_lines[c], "as validation"))
  ## Define test and training set, leaving out one cell line at a time
  train_data <- predictions %>%
    filter(cell_line != cell_lines[c])
  val_train <- validation_data %>%
    filter(cell_line != cell_lines[c])
  
  test_data <- predictions %>%
    filter(cell_line == cell_lines[c])
  val_test <- validation_data %>%
    filter(cell_line == cell_lines[c])
  
  ## Train the combination models
  
  # Normailsation of the weights
  model_weights <- train_data %>%
    group_by(subID) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(score = map(data, ~score_sc1(., val_train))) %>% 
    select(subID, score)  %>%
    unnest() %>%
    mutate(min_max_score = max(score) - score,
           weight = min_max_score / sum(min_max_score)) %>%
    select(subID, weight)
  
  # Train linear model
  lm_formula <- as.formula(paste0("val ~", paste(paste0("`", sub(".csv", "", selected_submissions), "`"), collapse =  "+")))
  lm_data <- train_data %>%
    gather(marker, prediction, reporters) %>% 
    spread(subID, prediction) %>% 
    full_join(val_train %>% gather(marker, val, reporters)) %>%
    arrange(glob_cellID, cell_line, treatment, time, cellID, fileID, marker)
  linear_model <- lm(lm_formula, data = lm_data)
  
  
  ## Score the different models and single models on the test set
  
  # Model with normalised weights
  weighted_combi_test_score <- test_data %>%
    left_join(model_weights) %>%
    gather(marker,marker_value, reporters) %>%
    group_by(glob_cellID, cell_line, treatment, time, cellID, fileID, marker) %>%
    summarise(prediction = sum(weight * marker_value)) %>%
    score_sc1_long_format(val_test)
  
  # Predict and score linear model
  lm_test <- test_data %>%
    gather(marker, prediction, reporters) %>% 
    spread(subID, prediction) %>%
    arrange(glob_cellID, cell_line, treatment, time, cellID, fileID, marker)
  lm_test_score <- predict(linear_model, lm_test) %>% 
    as_tibble() %>%
    bind_cols(select(lm_test, glob_cellID, cell_line, treatment, time, cellID, fileID, marker)) %>%
    rename(prediction = value) %>%
    score_sc1_long_format(val_test)
  
  # Single models
  best_single_test_score <- test_data %>% 
    group_by(subID) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(score = map(data, ~score_sc1(., val_test))) %>%
    select(subID, score) %>%
    unnest() %>%
    filter(score == min(score)) %>%
    pull(score)
  
  # Mean of models
  mean_score <- test_data  %>%
    gather(marker,marker_value, reporters) %>%
    group_by(glob_cellID, cell_line, treatment, time, cellID, fileID, marker) %>%
    summarise(prediction = mean(marker_value)) %>%
    score_sc1_long_format(val_test)
  
  # Median of models
  median_score <-  test_data  %>%
    gather(marker,marker_value,reporters)  %>%
    group_by(glob_cellID, cell_line, treatment, time, cellID, fileID, marker) %>%
    summarise(prediction = median(marker_value)) %>%
    score_sc1_long_format(val_test)
  
  CV_scores <- CV_scores %>% add_row("fold" = c, 
                                     "mean" = mean_score, 
                                     "median" = median_score, 
                                     "best_single" = best_single_test_score, 
                                     "weighted" = weighted_combi_test_score,
                                     "lm" = lm_test_score)
  
}

CV_scores <- CV_scores %>% filter(!is.na(fold))
saveRDS(CV_scores, "prediction_combinations/SC1/CV_model_scores.rds")


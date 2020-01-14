# Cross validated combination of models for SC2

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell") 

source("./scoring_scripts/score_sc2_noFile.R")

required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 

submission_folder <- "./submission_data/final/SC2"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc2.csv"))

validation_data <- read_csv("~/Seafile/validation_data/sc2gold.csv") %>% select(required_columns)
cell_lines <- validation_data %>% pull(cell_line) %>% unique

submissions <- list.files(submission_folder)
submissions <- submissions[!startsWith(submissions, "leaderboard")]

# For setting up and checking workflow select only three submissions
n_samples <- 3
selected_submissions <- sample(submissions, n_samples, replace = FALSE)

# Read prediction files
predictions <- tibble()
for (i in 1:length(selected_submissions)) {
  print(paste0(i, " out of ", length(selected_submissions)))
  
  pred <- read_csv(file.path(submission_folder,selected_submissions[i])) %>%
    as_tibble() %>%
    select(required_columns) %>% 
    add_column(subID = as.numeric(sub(".csv", "", selected_submissions[i])))
  predictions <- bind_rows(predictions, pred)
}

CV_scores <- tibble("fold" = NA, 
                    "equal" = NA, 
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
    mutate(score = map(data, ~score_sc2(., val_train)))  %>%
    select(subID, score) %>%
    unnest() %>% 
    mutate(min_max_score = max(score) - score,
           weight = min_max_score / sum(min_max_score)) %>%
    select(subID, weight)
  
  # Train linear model
  lm_formula <- as.formula(paste0("val ~", paste(paste0("`", sub(".csv", "", submissions), "`"), collapse =  "+")))
  lm_data <- train_data %>%
    gather(marker, prediction, -cell_line, -treatment, -time, -subID) %>% 
    spread(subID, prediction) %>% 
    full_join(val_train %>% gather(marker, val, -cell_line, -treatment, -time)) %>%
    arrange(cell_line, treatment, time, marker)
  linear_model <- lm(lm_formula, data = lm_data)
  
  
  ## Score the different models and single models on the test set
  
  # Model with normalised weights
  weighted_combi_test_score <- test_data  %>%
    group_by(subID) %>% 
    nest() %>% 
    ungroup() %>%
    left_join(model_weights) %>%
    mutate(subsample = map2(data, weight, function(x,y) {x %>% 
               group_by(cell_line, treatment, time) %>% 
               sample_frac(y)})) %>%
    select(subsample) %>%
    unnest() %>%
    score_sc2(validation_data)
  
  # Predict and score linear model
  lm_test <- test_data %>%
    gather(marker, prediction, -cell_line, -treatment, -time, -subID) %>% 
    spread(subID, prediction) %>%
    arrange(cell_line, treatment, time, marker)
  lm_test_score <- predict(linear_model, lm_test) %>% 
    as_tibble() %>%
    bind_cols(select(lm_test, cell_line, treatment, time, marker)) %>%
    rename(prediction = value) %>%
    score_sc2_long_format(val_test)
  
  # Single models
  best_single_test_score <- test_data %>% 
    group_by(subID) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(score = map(data, ~score_sc2(., val_test))) %>%
    select(subID, score) %>%
    unnest() %>%
    filter(score == min(score)) %>%
    pull(score)
  
  # Equal sample of models
  equal_sample_size <- test_data %>% 
    group_by(cell_line, treatment, time, subID) %>%
    sample_frac(1/length(unique(test_data$subID))) %>%
    ungroup() %>%
    score_sc2(val_test)
  
  CV_scores <- CV_scores %>% add_row("fold" = c, 
                                     "equal" = equal_sample_size, 
                                     "best_single" = best_single_test_score, 
                                     "weighted" = weighted_combi_test_score,
                                     "lm" = lm_test_score)
  
}

CV_scores <- CV_scores %>% filter(!is.na(fold))
saveRDS(CV_scores, "prediction_combinations/SC2/CV_model_scores.rds")


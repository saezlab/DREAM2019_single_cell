# Cross validated combination of models for SC4

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell") 

source("./scoring_scripts/score_sc4_noFile.R")

required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 
# Scoring function wtth input a tibble, not a csv file
score_sc4_long_format <- function(prediction_data_long, validation_data){
  
  ### Formating -------------------------
  validation_data_long <- validation_data %>% gather(marker, test,-cell_line, -treatment, -time )
  # join the test and validation data
  
  combined_data = full_join(prediction_data_long,
                            validation_data_long,
                            by=c("cell_line", "treatment", "time","marker"))
  
  ### Calculate score --------------------
  # calculate the  distance over all stats
  RMSE_cond = combined_data %>% group_by(cell_line,treatment,marker) %>% 
    summarise(RMSE = sqrt(sum((test - prediction)^2)/n())) 
  
  final_score = mean(RMSE_cond$RMSE)
}

submission_folder <- "./submission_data/final/SC4"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc4.csv")) 

validation_data <- read_csv("~/Seafile/validation_data/sc4gold.csv") %>% select(required_columns)
cell_lines <- validation_data %>% pull(cell_line) %>% unique

submissions <- list.files(submission_folder)
submissions <- submissions[!startsWith(submissions, "leaderboard")]

# For setting up and checking workflow select only three submissions
n_samples <- 3
selected_submissions <- sample(submissions, n_samples, replace = FALSE)

# Read prediction files
predictions <- tibble()
for (i in 1:length(submissions)) {
  print(paste0(i, " out of ", length(submissions)))
  
  pred <- read_csv(file.path(submission_folder,submissions[i])) %>%
    as.tibble() %>%
    select(required_columns) %>% 
    add_column(subID = as.numeric(sub(".csv", "", submissions[i])))
  predictions <- bind_rows(predictions, pred)
}

single_best_sub <- leader_board %>%
  select(objectId, submitterId, score) %>%
  filter(score == min(score)) %>%
  rename(subID = objectId)
 
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
    mutate(score = map(data, ~score_sc4(., val_train)))  %>%
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
  weighted_combi_test_score <- test_data %>%
    left_join(model_weights) %>%
    gather(marker,marker_value,-cell_line, -treatment, -time, -subID, -weight) %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(prediction = sum(weight * marker_value)) %>%
    score_sc4_long_format(val_test)
  
  # Predict and score linear model
  lm_test <- test_data %>%
    gather(marker, prediction, -cell_line, -treatment, -time, -subID) %>% 
    spread(subID, prediction) %>%
    arrange(cell_line, treatment, time, marker)
  lm_test_score <- predict(linear_model, lm_test) %>% 
    as_tibble() %>%
    bind_cols(select(lm_test, cell_line, treatment, time, marker)) %>%
    rename(prediction = value) %>%
    score_sc4_long_format(val_test)
  
  ## Single models
  # Determine single best model on test set
  if (FALSE) {
    best_single_test_score <- test_data %>% 
      group_by(subID) %>% 
      nest() %>% 
      ungroup() %>% 
      mutate(score = map(data, ~score_sc4(., val_test))) %>%
      select(subID, score) %>%
      unnest() %>%
      filter(score == min(score)) %>%
      pull(score)
  }
  
  # Score overall single best model on test set
  if (TRUE) {
    best_single_test_score <- test_data %>% 
      filter(subID == single_best_sub$subID) %>%
      select(required_columns) %>%
      score_sc4(val_test)
  }
  
  # Mean of models
  mean_score <- test_data  %>%
    gather(marker,marker_value,-cell_line, -treatment, -time, -subID) %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(prediction = mean(marker_value)) %>%
    score_sc4_long_format(val_test)
  
  # Median of models
  median_score <-  test_data  %>%
    gather(marker,marker_value,-cell_line, -treatment, -time, -subID)  %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(prediction = median(marker_value)) %>%
    score_sc4_long_format(val_test)
  
 CV_scores <- CV_scores %>% add_row("fold" = c, 
                       "mean" = mean_score, 
                       "median" = median_score, 
                       "best_single" = best_single_test_score, 
                       "weighted" = weighted_combi_test_score,
                       "lm" = lm_test_score)
  
}

CV_scores <- CV_scores %>% filter(!is.na(fold))
if (FALSE) {saveRDS(CV_scores, "prediction_combinations/SC4/CV_model_scores.rds")}

#CV_scores <- readRDS("prediction_combinations/SC4/CV_model_scores.rds")

best_single <- leader_board %>%
  rename(subID = submitterId) %>%
  select(subID, score) %>%
  filter(score == min(score))

CV_scores %>%
  gather(model, score, -fold) %>%
  ggplot(aes(model, score, fill=model)) +
  geom_boxplot()  +
  labs(title = "SC4; 5 fold cross validation")





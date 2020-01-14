# Combine predictions from ramdomly selected submissions
# Combine different numbers of randomly selected submissions

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
source("./scoring_scripts/score_sc4_noFile.R")

submission_folder <- "./submission_data/final/SC4"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc4.csv"))


required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 

# Score predictions if predictions and validation are already in long format.
score_sc4_long_format <- function(prediction_data_long, validation_data) {
  
  validation_data_long <- validation_data %>% gather(marker,test,-cell_line, -treatment, -time )
  combined_data = full_join(prediction_data_long,
                            validation_data_long,
                            by=c("cell_line", "treatment", "time","marker"))
  
  ### Calculate score --------------------
  # calculate the  distance over all stats
  RMSE_cond = combined_data %>% group_by(cell_line,treatment,marker) %>% 
    summarise(RMSE = sqrt(sum((test - prediction)^2)/n())) 
  
  final_score = mean(RMSE_cond$RMSE)
}
validation_data <- read_csv("~/Seafile/validation_data/sc4gold.csv") %>% select(required_columns)

submissions <- list.files(submission_folder) 
submissions <- submissions[!startsWith(submissions, "leaderboard")]

# Read all predictions
all_predictions <- lapply(submissions, function(x) {
  read_csv(file.path(submission_folder,x))  %>%
    select(required_columns) %>% 
    add_column(subID = as.numeric(sub(".csv", "", x)))
})

n_samples <- 10
n_iter<- 10

# Perform sampling n_samples random submissions n_iter number of times
# Score the combinations of the submission and keep track of the scores

repeated_scores <- tibble()
for (j in 1:n_iter) {
  print(j)

  # Select the sampled predictions
  predictions <- bind_rows(sample(all_predictions, n_samples, replace = FALSE))
  
  predictions_long <- predictions %>% 
    gather(marker,marker_value,-cell_line, -treatment, -time, -subID)
  
  # Get scores from the selected predictions 
  scores <- leader_board %>%
    filter(objectId %in% predictions$subID)  %>%
    rename(subID = objectId) %>%
    mutate(min_max_score = max(score) - score,
           weight = min_max_score / sum(min_max_score)) %>%
    select(subID, score, weight)
  
  mean_score <- predictions_long  %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(prediction = mean(marker_value)) %>%
    score_sc4_long_format(validation_data)
  
  median_score <- predictions_long  %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(prediction = median(marker_value)) %>%
    score_sc4_long_format(validation_data)
  
  weighted_by_score <- predictions_long %>%
    left_join(scores) %>% 
    group_by(cell_line, treatment, time, marker) %>%
    summarise(prediction = sum(weight * marker_value)) %>%
    score_sc4_long_format(validation_data)

  all_scores <- scores %>% 
    rename(model = subID) %>%
    add_row(model = c("Mean", "Median", "Weighted"), 
                                   score = c(mean_score, median_score, weighted_by_score)) %>%
    select(model, score) %>%
    arrange(score)
  print(all_scores)
  
  repeated_scores[j, 1] <- mean_score
  repeated_scores[j, 2] <- median_score
  repeated_scores[j, 3] <- weighted_by_score
  
}
colnames(repeated_scores) <- c("Mean", "Median", "Weighted")
#repeated_scores <- repeated_scores %>% as_tibble()
saveRDS(repeated_scores, paste0("prediction_combinations/SC4/SC4_", n_samples, "random_subs_scores.rds"))

# Combine predictions from ramdomly selected submissions
# Combine different numbers of randomly selected submissions

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
source("./scoring_scripts/score_sc1_noFile.R")

submission_folder <- "./submission_data/final/SC1"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_all_rounds_sc1.csv"))

required_columns <- c("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID", 
                      "p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")
reporters <- c("p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")

# Score predictions if predictions and validation are already in long format.
score_sc1_long_format <- function(prediction_data_long, validation_data_long) {
  
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

validation_data <- read_csv("~/Seafile/validation_data/sc1gold.csv") %>% select(required_columns)
validation_data_long <- validation_data %>% gather(marker, test, reporters)

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
  
  predictions_long <- predictions%>% 
    gather(marker, marker_value, reporters)
  
  # Get scores from the selected predictions 
  scores <- leader_board %>%
    filter(objectId %in% predictions$prediction)  %>%
    rename(prediction = objectId) %>%
    mutate(min_max_score = max(score) - score,
           weight = min_max_score / sum(min_max_score)) %>%
    select(prediction, score, weight)
  
    mean_score <- predictions_long  %>%
      group_by(glob_cellID, cell_line, treatment, time, cellID, fileID, marker) %>%
      summarise(prediction = mean(marker_value)) %>%
      score_sc1_long_format(validation_data_long)
    
    median_score <- predictions_long  %>%
      group_by(glob_cellID, cell_line, treatment, time, cellID, fileID, marker) %>%
      summarise(prediction = median(marker_value)) %>%
      score_sc1_long_format(validation_data_long)
    
    weighted_by_score <- predictions_long %>%
      left_join(scores) %>% 
      group_by(glob_cellID, cell_line, treatment, time, cellID, fileID, marker) %>%
      summarise(prediction = sum(weight * marker_value)) %>%
      score_sc1_long_format(validation_data_long)
    
    all_scores <- scores %>% add_row(prediction = c("Mean", "Median", "Weighted"), 
                                     score = c(mean_score, median_score, weighted_by_score)) %>%
      select(prediction, score) %>%
      arrange(score)
    print(all_scores)
    repeated_scores[j, 1] <- mean_score
    repeated_scores[j, 2] <- median_score
    repeated_scores[j, 3] <- weighted_by_score
    
}
colnames(repeated_scores) <- c("Mean", "Median", "Weighted")
#repeated_scores <- repeated_scores %>% as_tibble()
saveRDS(repeated_scores, paste0("prediction_combinations/SC1/SC1_", n_samples, "random_subs_scores.rds"))

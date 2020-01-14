# To combine the predictions in subchallenge 1
library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
source("./scoring_scripts/score_sc1_noFile.R")

submission_folder <- "./submission_data/round3/SC1"

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

# Read prediction files
for (i in 1:length(list.files(submission_folder))) {
  pred <- read_csv(file.path(submission_folder, list.files(submission_folder)[i])) %>%
    as.tibble() %>%
    select(required_columns) %>% 
    add_column(prediction = paste0("pred", i))
  assign(paste0("pred", i), pred)
}

predictions <- pred3 %>%
  bind_rows(pred2, pred1) 

predictions_long <- predictions%>% 
  gather(marker, marker_value, reporters)

# Score the submitted predictions
if (FALSE) {
  scores <- predictions %>% 
    group_by(prediction) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(score = map(data, ~score_sc1(., validation_data))) %>%
    select(prediction, score) %>%
    unnest() %>%
    mutate(min_max_score = max(score) - score,
           weight = min_max_score / sum(min_max_score))
  saveRDS(scores, "prediction_combinations/SC1/SC1_prediction_scores.rds")
} else{scores <- readRDS("prediction_combinations/SC1/SC1_prediction_scores.rds")}


# Perform combing the scores a number of times, for ths SC so far only once, since score is always the same
# Keep track of how well a combing/scoring model does. 
# How often is it best and how ofter better than best single model?

if (FALSE) {
  repeated_scores <- data.frame()
  better_and_best <- data.frame(count = rep(0, 7), row.names = c("Mean_better", "Mean_best",
                                                                 "Median_beter", "Median_best",
                                                                 "Weight_better", "Weight_best",
                                                                 "Single_best"))
for (i in 1:1) {
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
  repeated_scores[i, 1] <- mean_score
  repeated_scores[i, 2] <- median_score
  repeated_scores[i, 3] <- weighted_by_score
  
  if (mean_score < min(scores$score)) {better_and_best[1, 1] <- better_and_best[1, 1]+1 }
  if (mean_score == min(all_scores$score)) {better_and_best[2, 1] <- better_and_best[2, 1]+1 }
  if (median_score < min(scores$score)) {better_and_best[3, 1] <- better_and_best[3, 1]+1 }
  if (median_score == min(all_scores$score)) {better_and_best[4, 1] <- better_and_best[4, 1]+1 }
  if (weighted_by_score < min(scores$score)) {better_and_best[5, 1] <- better_and_best[5, 1]+1 }
  if (weighted_by_score == min(all_scores$score)) {better_and_best[6, 1] <- better_and_best[6, 1]+1 }
  if (min(scores$score) == min(all_scores$score)) {better_and_best[7, 1] <- better_and_best[7, 1]+1 }
}
colnames(repeated_scores) <- c("Mean", "Median", "Weighted")
repeated_scores <- repeated_scores %>% as_tibble()
better_and_best <- better_and_best %>% rownames_to_column() %>% as_tibble() %>% rename(measure = rowname)

saveRDS(repeated_scores, "prediction_combinations/SC1/SC1_scores.rds")
saveRDS(better_and_best, "prediction_combinations/SC1/SC1_performance.rds")
} else {
  repeated_scores <- readRDS("prediction_combinations/SC1/SC1_scores.rds")
  better_and_best <- readRDS("prediction_combinations/SC1/SC1_performance.rds")
}

# PLot results
repeated_scores %>%
  gather(model, score) %>%
  ggplot(aes(model, score, colour = model)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = min(scores$score))) +
  theme(legend.position = "none") +
  scale_y_reverse() +
  labs(title = "SC1")

better_and_best %>%
  ggplot(aes(measure, count, fill = measure)) +
  geom_col() +
  theme(legend.position = "none") +
  labs(title = "SC1")



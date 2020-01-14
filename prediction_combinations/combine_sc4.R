# For combining predictions in SC4

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

source("./scoring_scripts/score_sc4_noFile.R")

# Scoring function wtth input a tibble, not a csv file
score_sc4_long_format <- function(prediction_data_long, validation_data_long){
  
  ### Formating -------------------------
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

subchallenge <-  "prediction_combinations/SC4"
submission_folder <- "./submission_data/round3/SC4"
validation_file <- "~/Seafile/validation_data/sc4gold.csv"

required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 

validation_data <- read_csv(validation_file) %>% select(required_columns)
validation_data_long <- validation_data %>% gather(marker,test,-cell_line, -treatment, -time )

# Read prediction files
for (i in 1:length(list.files(submission_folder))) {
  pred <- read_csv(file.path(submission_folder, list.files(submission_folder)[i])) %>%
    as_tibble() %>%
    select(required_columns) %>% 
    add_column(prediction =paste0("pred", i))
  assign(paste0("pred", i), pred)
}

predictions <- pred3 %>%
  bind_rows(pred2, pred1)

predictions_long <- predictions %>% gather(marker,marker_value,-cell_line, -treatment, -time, -prediction)

# Score the submitted predictions
if (FALSE) {
  scores <- predictions %>% 
    group_by(prediction) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(score = map(data, ~score_sc4(., validation_data))) %>%
    select(prediction, score) %>%
    unnest() %>%
    mutate(min_max_score = max(score) - score,
           weight = min_max_score / sum(min_max_score))
  saveRDS(scores, "prediction_combinations/SC4/SC4_prediction_scores.rds")
} else{scores <- readRDS("prediction_combinations/SC4/SC4_prediction_scores.rds")}

# Perform combing the scores a number of times, for ths SC so far only once, since score is always the same
# Keep track of how well a combing/scoring model does. 
# How ogten is it best and how ofter better than best single model?
if (FALSE) {
  repeated_scores <- data.frame()
  better_and_best <- data.frame(count = rep(0, 7), row.names = c("Mean_better", "Mean_best",
                                                                 "Median_beter", "Median_best",
                                                                 "Weight_better", "Weight_best",
                                                                 "Single_best"))
  for (i in 1:1) {
    mean_score <- predictions_long  %>%
      group_by(cell_line, treatment, time, marker) %>%
      summarise(prediction = mean(marker_value)) %>%
      score_sc4_long_format(validation_data_long)
    
    median_score <- predictions_long  %>%
      group_by(cell_line, treatment, time, marker) %>%
      summarise(prediction = median(marker_value)) %>%
      score_sc4_long_format(validation_data_long)
    
    weighted_by_score <- predictions_long %>%
      left_join(scores) %>% 
      group_by(cell_line, treatment, time, marker) %>%
      summarise(prediction = sum(weight * marker_value)) %>%
      score_sc4_long_format(validation_data_long)
    
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
  
  saveRDS(repeated_scores, "prediction_combinations/SC4/SC4_scores.rds")
  saveRDS(better_and_best, "prediction_combinations/SC4/SC4_performance.rds")

} else {
  repeated_scores <- readRDS("prediction_combinations/SC4/SC4_scores.rds")
  better_and_best <- readRDS("prediction_combinations/SC4/SC4_performance.rds")
}
# PLot results
repeated_scores %>%
  gather(model, score) %>%
  ggplot(aes(model, score, colour = model)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = min(scores$score))) +
  theme(legend.position = "none") +
  scale_y_reverse() +
  labs(title = "SC4")

better_and_best %>%
  ggplot(aes(measure, count, fill = measure)) +
  geom_col() +
  theme(legend.position = "none") +
  labs(title = "SC4")


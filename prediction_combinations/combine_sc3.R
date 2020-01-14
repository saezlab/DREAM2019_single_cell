# To combine predictions of round 3
library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
source("./scoring_scripts/score_sc3_noFile.R")

subchallenge <-  "prediction_combinations/SC3"
submission_folder <- "./submission_data/round3/SC3"
validation_file <- "~/Seafile/validation_data/sc3gold.csv"

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

# Read prediction files
for (i in 1:length(list.files(submission_folder))) {
  pred <- read_csv(file.path(submission_folder, list.files(submission_folder)[i])) %>%
    as.tibble() %>%
    select(required_columns) %>% 
    add_column(prediction =paste0("pred", i))
  assign(paste0("pred", i), pred)
}

predictions <- pred3 %>%
  bind_rows(pred2, pred1)

# Score the submitted predictions
if (FALSE) {
  scores <- predictions %>% 
    group_by(prediction) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(score = map(data, ~score_sc3(., validation_data))) %>%
    select(prediction, score) %>%
    unnest() 
  saveRDS(scores, "prediction_combinations/SC3/SC3_prediction_scores.rds")
} else{scores <- readRDS("prediction_combinations/SC3/SC3_prediction_scores.rds")}

# Perform combing the scores a number of times, because the different samplings lead to different results
# Keep track of how well a combing/scoring model does. 
# How ogten is it best and how ofter better than best single model?
if (FALSE) {
repeated_scores <- data.frame()
better_and_best <- data.frame(count = rep(0, 5), row.names = c("Equal_better", "Equal_best",
                                                               "Weight_better", "Weight_best",
                                                               "Single_best"))

for (i in 1:10) {
  equal_sample_size <- predictions %>% 
    group_by(cell_line, treatment, time, prediction) %>%
    sample_frac(1/length(unique(predictions$prediction))) %>%
    ungroup() %>%
    score_sc3(validation_data)
  
  
  weighted_by_score <- predictions %>%
    group_by(prediction) %>%
    nest() %>%
    ungroup() %>%
    left_join(scores) %>%
    mutate(min_max_score = max(score) - score,
           weight = min_max_score / sum(min_max_score),
           subsample = map2(data, weight, function(x,y) {x %>% 
               group_by(cell_line, treatment, time) %>% 
               sample_frac(y)})) %>%
    select(subsample) %>%
    unnest() %>%
    score_sc3(validation_data)
  all_scores <- scores %>% add_row(prediction = c("Equal", "Weighted"), 
                                   score = c(equal_sample_size, weighted_by_score)) %>%
    arrange(score)
  print(all_scores)
  repeated_scores[i, 1] <- equal_sample_size
  repeated_scores[i, 2] <- weighted_by_score
  if (equal_sample_size < min(scores$score)) {better_and_best[1, 1] <- better_and_best[1, 1]+1 }
  if (equal_sample_size == min(all_scores$score)) {better_and_best[2, 1] <- better_and_best[2, 1]+1 }
  if (weighted_by_score < min(scores$score)) {better_and_best[3, 1] <- better_and_best[3, 1]+1 }
  if (weighted_by_score == min(all_scores$score)) {better_and_best[4, 1] <- better_and_best[4, 1]+1 }
  if (min(scores$score) == min(all_scores$score)) {better_and_best[5, 1] <- better_and_best[5, 1]+1 }
}
colnames(repeated_scores) <- c("Equal", "Weighted")
repeated_scores <- repeated_scores %>% as_tibble()
better_and_best <- better_and_best %>% rownames_to_column() %>% as_tibble() %>% rename(measure = rowname)

saveRDS(repeated_scores, "prediction_combinations/SC2/SC2_scores.rds")
saveRDS(better_and_best, "prediction_combinations/SC2/SC2_performance.rds")

} else {
  repeated_scores <- readRDS("prediction_combinations/SC3/SC3_scores.rds")
  better_and_best <- readRDS("prediction_combinations/SC3/SC3_performance.rds")
}

# PLot results
repeated_scores %>%
  gather(model, score) %>%
  ggplot(aes(model, score, fill = model)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = min(scores$score))) +
  theme(legend.position = "none") +
  scale_y_reverse() +
  labs(title = "SC3")
 
better_and_best %>%
  ggplot(aes(measure, count, fill = measure)) +
  geom_col() +
  theme(legend.position = "none") +
  labs(title = "SC3")

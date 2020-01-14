## SC3; without inner loop
# Combine predictions from ramdomly selected submissions
# Combine different numbers of randomly selected submissions

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function with input a tibble, not a csv file
source("./scoring_scripts/score_sc3_noFile.R")

required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 

validation_data <- read_csv("~/Seafile/validation_data/sc3gold.csv") 
validation_stats <- validation_data  %>% 
  select(required_columns) %>%
  data_to_stats() %>% 
  rename(test_stat_value = stat_value)  %>%
  arrange(cell_line, treatment, time, stat_variable) %>%
  ungroup() %>%
  select(test_stat_value)
submission_folder <- "./submission_data/final/SC3"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc3.csv"))
nested_predictions <- readRDS("./submission_data/intermediate_data/sc3_all_predictions_nested.rds") %>%
  ungroup()

repeated_scores <- tibble("Sample_size" = NA, "Iteration" = NA, "Equal_sample" = NA)
for (n in 1:dim(leader_board)[1]) {
  n_samples <- n
  n_iter<- 100
  
  # Perform sampling n_samples random submissions n_iter number of times
  # Each combination of sumbmissions is combined and scored ten times
  # Because sampling can differ each time. All scores are saves
  
  for (j in 1:n_iter) {
    print(paste0("samples: ", n_samples, " iteration: ", j))
    
    # Select the sampled predictions
    predictions <- sample(unique(nested_predictions$team), n_samples, replace = FALSE)
    
    
    equal_sample_size <- nested_predictions %>%
      filter(team %in% predictions) %>%
      mutate(sample = map(data, ~sample_frac(., 1/n_samples))) %>%
      select(-c(data, team)) %>%
      unnest(sample) %>%
      select(required_columns) %>%
      data_to_stats() %>% 
      rename(pred_stat_value = stat_value) %>%
      arrange(cell_line, treatment, time, stat_variable) %>%
      ungroup() %>%
      bind_cols(validation_stats) %>%
      select(test_stat_value, pred_stat_value) %>%
      as.matrix() %>%
      t() %>%
      dist(method = "euclidean")
    
    print(equal_sample_size)
    
    repeated_scores <- repeated_scores %>%
      add_row(Sample_size = n, 
              Iteration = j, 
              Equal_sample = equal_sample_size)
    
  }
  
}

repeated_scores <- filter(repeated_scores, !is.na(Sample_size))

if (FALSE) {saveRDS(repeated_scores, "prediction_combinations/SC3/SC3_random_subs_scores.rds")}

scores <- readRDS("prediction_combinations/SC3/SC3_random_subs_scores.rds")
ordered_combi <- readRDS("prediction_combinations/SC3/SC3_ordered_combination.rds") %>%
  gather(model, score, Equal_sample)
best_single <- leader_board  %>%
  rename(subID = objectId) %>%
  select(subID, score) %>%
  filter(score == min(score))

best_single <- leader_board  %>%
  rename(subID = objectId) %>%
  select(submitterId, score) %>%
  filter(score == min(score))


scores %>%
  gather(model, score, Equal_sample) %>%
  ggplot(aes(Sample_size, score, colour=model)) +
  geom_point()  +
  labs(title = "SC3", x= "Number of combined predictions")

scores %>%
  gather(model, score, Equal_sample) %>%
  ggplot(aes(as.factor(Sample_size), score, fill = model)) +
  geom_boxplot(position = "dodge") +
  xlim(as.character(1:14)) +
  geom_hline(aes(yintercept = best_single$score, colour = as.factor(best_single$submitterId))) +
  labs(title = "SC3", x= "Number of combined predictions", colour = "best single")  +
  geom_point(data=ordered_combi, aes(x=n, y=score, shape="top n")) +
  scale_shape_manual(values = 23) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

scores %>%
  gather(model, score, Equal_sample) %>%
  group_by(Sample_size, model) %>%
  summarise(score = mean(score)) %>%
  ggplot(aes(Sample_size, score, colour=model)) +
  geom_point() +
  geom_smooth() +
  geom_hline(aes(yintercept = best_single$score, colour = as.factor(best_single$submitterId))) +
  labs(title = "SC3", x= "Number of combined predictions")   +
  geom_point(data=ordered_combi, aes(x=n, y=score, shape="top n"),  stroke =1) +
  scale_shape_manual(values = 3) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

scores %>%
  gather(model, score, Equal_sample) %>%
  arrange(score)

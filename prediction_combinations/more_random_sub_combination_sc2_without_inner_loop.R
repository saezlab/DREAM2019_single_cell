## SC2; without inner loop
# Combine predictions from different number of ramdomly selected submissions
# Combining by sampling cells equally from the selected submissions

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
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

validation_data <- read_csv("./challenge_data/validation_data/sc2gold.csv") %>% select(required_columns)
submission_folder <- "./submission_data/final/SC2"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc2.csv"))

# The predicted cells of the teams, nested so team-cell_line-Tr-Time-dataframe with predictions 
nested_predictions <- readRDS("./submission_data/intermediate_data/sc2_all_predictions_nested.rds") %>%
  ungroup()

# Radnomly sample 1 submission by sampling scores from the leaderboard
repeated_scores <- tibble("Sample_size" = rep(1, 100), "Iteration" = seq(1,100), 
                          "Equal_sample" = sample(leader_board$score, 100, replace = TRUE))

# Perform sampling n submission randomly, n_iter number of times
# Score the combinations of the submission and keep track of the scores
for (n in 2:dim(leader_board)[1]) {
  n_samples <- n
  n_iter<- 100
  
  for (j in 1:n_iter) {
    print(paste0("samples: :", n_samples, " iteration: ", j))
    
    # Select the sampled predictions
    predictions <- sample(unique(nested_predictions$team), n_samples, replace = FALSE)
    
    # Combine by sampling equally from the selected submissions and score it
    equal_sample_size <- nested_predictions %>%
      filter(team %in% predictions) %>%
      mutate(sample = map(data, ~sample_frac(., 1/n_samples))) %>%
      select(-c(data, team)) %>%
      unnest(sample) %>%
      select(required_columns) %>%
      score_sc2(validation_data)
    
    # Update all scores
    repeated_scores <- repeated_scores %>%
      add_row(Sample_size = n, 
              Iteration = j, 
              Equal_sample = equal_sample_size)
  }
  
}

if (FALSE) {saveRDS(repeated_scores, "prediction_combinations/SC2/SC2_random_subs_scores.rds")}

scores <- readRDS("prediction_combinations/SC2/SC2_random_subs_scores.rds")
best_single <- leader_board  %>%
  rename(subID = objectId) %>%
  select(subID, score) %>%
  filter(score == min(score))


scores %>%
  gather(model, score, Equal_sample) %>%
  ggplot(aes(Sample_size, score, colour=model)) +
  geom_point()  +
  labs(title = "SC2", x= "Number of combined predictions")

scores %>%
  gather(model, score, Equal_sample) %>%
  ggplot(aes(as.factor(Sample_size), score, fill = model)) +
  geom_boxplot(position = "dodge") +
  geom_hline(aes(yintercept = best_single$score, colour = as.factor(best_single$subID))) +
  labs(title = "SC2", x= "Number of combined predictions", colour = "best single")

scores %>%
  gather(model, score, Equal_sample) %>%
  group_by(Sample_size, model) %>%
  summarise(score = mean(score)) %>%
  ggplot(aes(Sample_size, score, colour=model)) +
  geom_point() +
  geom_smooth() +
  geom_hline(aes(yintercept = best_single$score, colour = as.factor(best_single$subID))) +
  labs(title = "SC2", x= "Number of combined predictions")  +
  geom_point(data=ordered_combi, aes(x=n, y=score, shape="top n"),  stroke =1) +
  scale_shape_manual(values = 3) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

scores %>%
  gather(model, score, Equal_sample) %>%
  arrange(score) 

scores %>% select(Iteration) %>% unique()


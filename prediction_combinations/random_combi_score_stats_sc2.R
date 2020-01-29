## SC2; without inner loop
# Combine predictions from different number of ramdomly selected submissions
# Combining by sampling cells equally from the selected submissions
# Since we observed that random sampleing 1 or a few submissions performs better than sampling more
# Here we investigate if the covariance get worse when combining submissions causing this phenomenon

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
leader_board <- read_csv("./submission_data/final/SC2/leaderboard_final_sc2.csv")


submissions <-  readRDS("./submission_data/intermediate_data/sc2_ranked_teams.rds") %>%
  as.character()
# The predicted cells of the teams, nested so team-cell_line-Tr-Time-dataframe with predictions 
nested_predictions <- readRDS("./submission_data/intermediate_data/sc2_all_predictions_nested.rds") %>%
  ungroup()
subs_stats <- readRDS("./submission_data/intermediate_data/sc2_values_median_EGF0.rds") %>%
  select(-starts_with("EGF0")) %>% 
  separate(stat_variable, c("stat", "variable"), sep="-")

subs_stats_scores <- subs_stats %>% 
  group_by(stat) %>%
  summarise_at(submissions,~ sum((standard - .)^2))

# Radnomly sample 1 submission
score_per_stat <- sample(subs_stats_scores[,2:17], 100, replace=TRUE) %>% 
  t() %>% 
  as_tibble() %>%
  bind_cols(tibble("Sample_size" = rep(1, 100), "Iteration" = seq(1,100))) %>%
  rename(cov = V1, mean =V2) %>%
  select(Sample_size, Iteration, cov, mean)

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
      data_to_stats() %>%
      arrange(cell_line, treatment,  time, stat_variable) %>%
      select(cell_line, treatment, time, stat_value) %>%
      bind_cols(select(subs_stats, stat, variable, standard)) %>%
      group_by(stat) %>%
      summarise(score = sum((standard - stat_value)^2))
    
    # Update all scores
    score_per_stat <- score_per_stat %>%
      add_row(Sample_size = n, 
              Iteration = j, 
              cov = equal_sample_size[1,2],
              mean = equal_sample_size[2,2])
  }
  
}

score_per_stat <- score_per_stat %>% unnest() %>% mutate(Sample_size = as.factor(Sample_size))

if (FALSE) {saveRDS(score_per_stat, "prediction_combinations/SC2/SC2_random_stat_scores_SSQ.rds")}

score_per_stat %>%
  gather(stat, score, cov, mean) %>%
  ggplot(aes(Sample_size, score, fill = stat)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title = "SC2", x= "Number of combined predictions")

score_per_stat %>%
  mutate(score = sqrt(cov+mean)) %>%
  ggplot(aes(Sample_size, score)) +
  geom_boxplot()

score_per_stat  %>%
  gather(stat, score, cov, mean) %>%
  group_by(Sample_size, Iteration) %>%
  mutate(mean_score = mean(score))  %>%
  ggplot(aes(Sample_size, score, fill = stat)) +
  geom_boxplot() +
  geom_boxplot(inherit.aes = FALSE, aes(Sample_size, mean_score), position = "dodge") +
  labs(title = "SC2", x= "Number of combined predictions")

score_per_stat  %>%
  filter(Sample_size %in% c(1,2)) %>%
  select(-Iteration) %>%
  distinct() %>%
  group_by(Sample_size) %>%
  sample_n(16, replace = FALSE) %>%
  ggplot(aes(cov, mean, colour=Sample_size)) +
  geom_point() +
  geom_abline(slope = 0.5) +
  labs(title = "SC2")

subs_stats_scores  %>%
  gather(team, SSQ, submissions) %>%
  spread(stat, SSQ) %>%
  ggplot(aes(cov, mean)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label=team)) + 
  scale_x_log10() +
  scale_y_log10() +
  geom_point(data = filter(score_per_stat, Sample_size !=1), aes(colour=Sample_size))

  
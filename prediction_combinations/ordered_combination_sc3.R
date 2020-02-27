## SC3
# Select the n best submissions and score the combination of those submissions

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
source("./scoring_scripts/score_sc3_noFile.R")

validation_data <- read_csv("./challenge_data/validation_data/sc3gold.csv") 
leader_board <- read_csv("./submission_data/final/SC3/leaderboard_final_sc3.csv")

# Submissions are ordered from best to worst team
submissions <- readRDS("./submission_data/intermediate_data/sc3_ranked_teams.rds") %>% as.character()

# The predicted cells of the teams, nested so team-cell_line-Tr-Time-dataframe with predictions 
nested_predictions <- readRDS("./submission_data/intermediate_data/sc3_all_predictions_nested.rds") %>%
  ungroup()

# Select the top teams, combine their predictions by taking the mean or median and scoring it
# Keeping track of the scores
all_scores <- tibble("n"=NA, "Equal_sample" = NA)
for (i in 1:length(submissions)) {
  print(paste0(i, " out of ", length(submissions)))
  
  selected_subs <- submissions[1:i]
  
  equal_sample_size <- nested_predictions %>%
    filter(team %in% selected_subs) %>%
    mutate(sample = map(data, ~sample_frac(., 1/i))) %>%
    select(-c(data, team)) %>%
    unnest(sample) %>%
    score_sc3(validation_data)
  
  all_scores <- all_scores %>%
    add_row(n = i,
            Equal_sample = equal_sample_size)
}

all_scores <- all_scores %>% filter(!is.na(n))
saveRDS(all_scores, "prediction_combinations/SC3/SC3_ordered_combination.rds")

equal_sample_size <- nested_predictions %>%
  filter(team %in% submissions[2]) %>%
  unnest(data) %>%
  score_sc3(validation_data)

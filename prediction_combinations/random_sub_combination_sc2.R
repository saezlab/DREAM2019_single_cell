# Combine predictions from ramdomly selected submissions
# Combine different numbers of randomly selected submissions

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
source("./scoring_scripts/score_sc2_noFile.R")

submission_folder <- "./submission_data/final/SC2"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc2.csv"))


required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 

validation_data <- read_csv("~/Seafile/validation_data/sc2gold.csv") %>% select(required_columns)

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
# Each combination of sumbmissions is combined and scored ten times
# Because sampling can differ each time. All scores are saves

repeated_scores <- tibble()
count = 1
for (j in 1:n_iter) {
  print(j)

  # Select the sampled predictions
  predictions <- bind_rows(sample(all_predictions, n_samples, replace = FALSE))
  
  # Get scores from the selected predictions 
  scores <- leader_board %>%
    filter(objectId %in% predictions$subID)  %>%
    rename(subID = objectId) %>%
    mutate(min_max_score = max(score) - score,
           weight = min_max_score / sum(min_max_score)) %>%
    select(subID, score, weight)
  
  for (k in 1:10) {
    print(paste0("Outerloop: ", j, "; inner loop: ", k, "; count : ", count))
    equal_sample_size <- predictions %>% 
      group_by(cell_line, treatment, time, subID) %>%
      sample_frac(1/length(unique(predictions$subID))) %>%
      ungroup() %>%
      score_sc2(validation_data)
    
    weighted_by_score <- predictions %>%
      group_by(subID) %>%
      nest() %>%
      ungroup() %>%
      left_join(scores) %>%
      mutate(subsample = map2(data, weight, function(x,y) {x %>% 
                 group_by(cell_line, treatment, time) %>% 
                 sample_frac(y)})) %>%
      select(subsample) %>%
      unnest() %>%
      score_sc2(validation_data)
    
    all_scores <- scores %>%
      select(subID, score) %>%
      rename(model = subID) %>%
      add_row(model = c("Equal", "Weighted"), 
              score = c(equal_sample_size, weighted_by_score)) %>%
      arrange(score)
    print(all_scores)
    
    repeated_scores[count, 1] <- j
    repeated_scores[count, 2] <- equal_sample_size
    repeated_scores[count, 3] <- weighted_by_score
    
    count = count + 1
  }
}

colnames(repeated_scores) <- c("Subsample_round", "Equal", "Weighted")
#repeated_scores <- repeated_scores %>% as_tibble()
saveRDS(repeated_scores, paste0("prediction_combinations/SC2/SC2_", n_samples, "random_subs_scores.rds"))

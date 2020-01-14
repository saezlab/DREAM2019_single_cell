## SC3
# Select the n best submissions and score the combination of those submissions

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
source("./scoring_scripts/score_sc3_noFile.R")

submission_folder <- "./submission_data/final/SC3"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc3.csv")) %>% 
  arrange(score) %>% 
  rownames_to_column(var ="rank")

required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 

validation_data <- read_csv("~/Seafile/validation_data/sc3gold.csv") %>% 
  select(required_columns)

submissions <- pull(leader_board, objectId)

# Read all predictions
all_predictions <- lapply(submissions, function(x) {
  read_csv(file.path(submission_folder,paste0(x, ".csv")))  %>%
    select(required_columns) %>% 
    add_column(subID =  x)
})

all_scores <- tibble("n"=NA, "Equal_sample" = NA)
for (i in 1:length(submissions)) {
  print(i)
  
  selected_subs <- all_predictions[1:i]
  
  equal_sample_size <- lapply(selected_subs, function(x) {x %>% 
      group_by(cell_line, treatment, time) %>%  sample_frac(1/i)}) %>%
    bind_rows() %>%
    select(required_columns) %>%
    score_sc3(validation_data)
  
  print(equal_sample_size)
  
  all_scores <- all_scores %>%
    add_row(n = i,
            Equal_sample = equal_sample_size)
  
}

all_scores <- all_scores %>% filter(!is.na(n))
saveRDS(all_scores, "prediction_combinations/SC3/SC3_ordered_combination.rds")
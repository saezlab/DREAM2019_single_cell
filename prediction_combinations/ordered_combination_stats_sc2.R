## SC2; combining by statisticss
# Select the n best submissions and score the combination of those submissions

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(Biobase)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

leader_board <- read_csv("./submission_data/final/SC2/leaderboard_final_sc2.csv")
submissions <-  readRDS("./submission_data/intermediate_data/sc2_ranked_teams.rds") %>%
  as.character()
# Statistics per condition of the standard and calculated from the predictions, including median EGF t=0 values per condition
sub_data_values <- readRDS("./submission_data/intermediate_data/sc2_values_median_EGF0.rds") %>%
  select(cell_line, treatment, time, stat_variable, standard, submissions)

# Radnomly sample 1 submission by sampling scores from the leaderboard
all_scores <- tibble("n"=NA, "median_stats" = NA)

for (i in 1:length(submissions)) {
  print(paste0(i, " out of ", length(submissions)))
  
  selected_subs <- submissions[1:i]
  
  equal_sample_size <-  sub_data_values %>%
    mutate(prediction = rowMedians(as.matrix(.[selected_subs]))) %>%
    select(-submissions) %>%
    group_by(cell_line, treatment, time) %>%
    summarise(SSQ = sum((standard - prediction)^2)) %>%
    ungroup() %>%
    pull(SSQ) %>%
    sum() %>%
    sqrt()
    
  all_scores <- all_scores %>%
    add_row(n = i,
            median_stats = equal_sample_size)
  
}
all_scores <- all_scores %>% filter(!is.na(n))
if (FALSE) {saveRDS(all_scores, "prediction_combinations/SC2/SC2_ordered_combination_stats.rds")}



## SC4
# Select the n best submissions and score the combination of those submissions

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
source("./scoring_scripts/score_sc4_noFile.R")
# Score predictions if predictions and validation are already in long format.
score_sc4_long_format <- function(prediction_data_long, validation_data) {
  
  validation_data_long <- validation_data %>% gather(marker,test,-cell_line, -treatment, -time )
  combined_data = full_join(prediction_data_long,
                            validation_data_long,
                            by=c("cell_line", "treatment", "time","marker"))
  
  ### Calculate score --------------------
  # calculate the  distance over all stats
  RMSE_cond = combined_data %>% group_by(cell_line,treatment,marker) %>% 
    summarise(RMSE = sqrt(sum((test - prediction)^2)/n())) 
  
  final_score = mean(RMSE_cond$RMSE)
}

submission_folder <- "./submission_data/final/SC4"
validation_data <- read_csv("~/Seafile/validation_data/sc4gold.csv") %>% select(required_columns)
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc4.csv")) %>% 
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

submissions <- pull(leader_board, objectId)

# Read all predictions
all_predictions <- lapply(submissions, function(x) {
  read_csv(file.path(submission_folder,paste0(x, ".csv")))  %>%
    select(required_columns) %>% 
    add_column(subID =  x)
})

all_scores <- tibble("n" =  NA, "Mean" = NA, "Median" = NA)
for (i in 1:length(submissions)) {
  print(i)
  
  selected_subs_long <- bind_rows(all_predictions[1:i]) %>%
    gather(marker,marker_value,-cell_line, -treatment, -time, -subID)
  
  mean_score <- selected_subs_long  %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(prediction = mean(marker_value)) %>%
    score_sc4_long_format(validation_data)
  
  median_score <- selected_subs_long  %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(prediction = median(marker_value)) %>%
    score_sc4_long_format(validation_data)
  
  all_scores <- all_scores %>% add_row(
    n = i,
    Mean = mean_score,
    Median = median_score)
 
}

all_scores <- all_scores %>% filter(!is.na(n))
saveRDS(all_scores, "prediction_combinations/SC4/SC4_ordered_combination.rds")

 
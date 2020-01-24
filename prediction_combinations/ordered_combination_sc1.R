## SC1
# Select the n best submissions and score the combination of those submissions

library(tidyverse)
library(Biobase)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# The predictions of each team are a column, one column per team
# Thecolumn standard is the goldesn standard/true value
# The columns are ordered, best performing team first and worst team last
all_predictions <- readRDS("./submission_data/intermediate_data/sc1_all_predictions.rds")

leader_board <- read_csv("./submission_data/final/SC1/leaderboard_final_sc1.csv")
submissions <- readRDS("./submission_data/intermediate_data/sc1_ranked_teams.rds") %>% 
  as.character()

# Select the top teams, combine their predictions by taking the mean or median and scoring it
# Keeping track of the scores
all_scores <- tibble("n" =  NA, "Mean" = NA, "Median" = NA)
for (i in 1:dim(leader_board)[1]) {
  print(paste0(i, " out of ", dim(leader_board)[1]))
  
  mean_score <- all_predictions %>%
    mutate(prediction = rowMeans(.[9:(8+i)])) %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(RMSE = sqrt(sum((standard - prediction)^2) / n())) %>%
    pull(RMSE) %>%
    mean()
  
  median_score <- all_predictions %>%
    mutate(prediction = rowMedians(as.matrix(.[9:(8+i)])))  %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(RMSE = sqrt(sum((standard - prediction)^2) / n())) %>%
    pull(RMSE) %>%
    mean()
  
  all_scores <- all_scores %>% add_row(
    n = i,
    Mean = mean_score,
    Median = median_score)

}

all_scores <- all_scores %>% filter(!is.na(n))
saveRDS(all_scores, "prediction_combinations/SC1/SC1_ordered_combination.rds")



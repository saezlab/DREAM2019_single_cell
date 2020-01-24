## SC4
# Select the n best submissions and score the combination of those submissions

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Read all predictions
all_predictions <- lapply(submissions, function(x) {
  read_csv(file.path(submission_folder,paste0(x, ".csv")))  %>%
    select(required_columns) %>% 
    add_column(subID =  x)
})

# The predictions of each team are a column, one column per team
# The column standard is the goldesn standard/true value
all_predictions <- readRDS("./submission_data/intermediate_data/sc4_combined_data.rds")

leader_board <- read_csv("./submission_data/final/SC4/leaderboard_final_sc4.csv")
submissions <- readRDS("./submission_data/intermediate_data/sc4_ranked_teams.rds") %>% 
  as.character()

# Select the top teams, combine their predictions by taking the mean or median and scoring it
# Keeping track of the scores
all_scores <- tibble("n" =  NA, "Median" = NA)
for (i in 1:length(submissions)) {
  print(paste0(i, " out of ", length(submissions)))
  
  median_score <- all_predictions %>%
    mutate(prediction = rowMedians(as.matrix(.[6:(5+i)])))  %>%
    group_by(cell_line,treatment,marker) %>% 
    summarise(RMSE = sqrt(sum((standard - prediction)^2)/n()))  %>%
    pull(RMSE) %>%
    mean()
  
  all_scores <- all_scores %>% add_row(
    n = i,
    Median = median_score)
 
}

all_scores <- all_scores %>% filter(!is.na(n))
saveRDS(all_scores, "prediction_combinations/SC4/SC4_ordered_combination.rds")

 
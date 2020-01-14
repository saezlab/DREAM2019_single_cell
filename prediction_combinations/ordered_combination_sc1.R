## SC1
# Select the n best submissions and score the combination of those submissions

library(tidyverse)
library(Biobase)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")


# Score predictions if predictions and validation are already in long format.
score_sc1_long_format <- function(combined_data) {
  
  #combined_data <- validation_data_long %>%
  #  full_join(prediction_data_long, by = c("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID", "marker"))
  
  ### Calculate score --------------------
  # calculate the RMSE for each condition
  RMSE_cond <- combined_data %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(RMSE = sqrt(sum((standard - prediction)^2) / n()))
  final_score <- mean(RMSE_cond$RMSE)
  return(final_score)
}

#' order_predictions
#'
#' takes a prediction matrix, reformats to long format and keeps only the values.
#' 
#' @param predictions prediction matrix in tibble format
#' @return value vector gathered over reporters and sorted by glob_cellID
order_predictions <- function(predictions){
  
  predictions %>% select(glob_cellID,reporters) %>% 
    arrange(glob_cellID) %>%
    gather("reporter","value",reporters) %>% select(value)
}

# order the predictions the same way for all teams and for the golden standard
# so a simple column bind will merge well.

required_columns <- c("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID", 
                      "p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")
reporters <- c("p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")

submission_folder <- "./submission_data/final/SC1"

leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc1.csv")) %>%
  select(-writeUp, -createdOn) %>% 
  mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
  mutate(submitterId = make.names(submitterId))
leader_board <- leader_board %>% filter(submitterId != "KAUST_RSS") 

validation_data_long <- read_csv("~/Seafile/validation_data/sc1gold.csv") %>%
  select(required_columns) %>%
  arrange(glob_cellID) %>%
  gather(marker, standard, reporters)

#submissions <- pull(leader_board, submitterId)

# read team's predictions from csv files
prediction_data <- leader_board  %>% 
  mutate(predictions = map(file.path(submission_folder,submissions),read_csv))

ordered_predictions <- prediction_data %>%
  mutate(predictions = map(predictions,order_predictions)) %>%
  do(bind_cols(.$predictions))
names(ordered_predictions) <- prediction_data$submitterId

combined_data <- bind_cols(validation_data_long, ordered_predictions)

all_scores <- tibble("n" =  NA, "Mean" = NA, "Median" = NA)
for (i in 1:dim(leader_board)[1]) {
  print(i)
  
  mean_score <- combined_data %>%
    mutate(prediction = rowMeans(.[9:(8+i)])) %>%
    select(glob_cellID, cell_line, treatment, time, cellID, fileID, marker, standard, prediction) %>%
    score_sc1_long_format()
  
  median_score <- combined_data %>%
    mutate(prediction = rowMedians(as.matrix(.[9:(8+i)]))) %>%
    select(glob_cellID, cell_line, treatment, time, cellID, fileID, marker, standard, prediction) %>%
    score_sc1_long_format()
  
  all_scores <- all_scores %>% add_row(
    n = i,
    Mean = mean_score,
    Median = median_score)

}

all_scores <- all_scores %>% filter(!is.na(n))
saveRDS(all_scores, "prediction_combinations/SC1/SC1_ordered_combination.rds")



## SC1
# Combine predictions from ramdomly selected submissions
# Combine different numbers of randomly selected submissions

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

submission_folder <- "./submission_data/final/SC1"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc1.csv"))

required_columns <- c("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID", 
                      "p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")
reporters <- c("p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")

# Score predictions if predictions and validation are already in long format.
score_sc1_long_format <- function(prediction_data_long, validation_data_long) {
  
  combined_data <- validation_data_long %>%
    full_join(prediction_data_long, by = c("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID", "marker"))
  
  ### Calculate score --------------------
  # calculate the RMSE for each condition
  RMSE_cond <- combined_data %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(RMSE = sqrt(sum((test - prediction)^2) / n()))
  final_score <- mean(RMSE_cond$RMSE)
  return(final_score)
}

validation_data <- read_csv("~/Seafile/validation_data/sc1gold.csv") %>% select(required_columns)
validation_data_long <- validation_data %>% gather(marker, test, reporters)

submissions <- pull(leader_board, objectId)

# Read all predictions
all_predictions <- lapply(submissions, function(x) {
  read_csv(file.path(submission_folder,paste0(x, ".csv")))  %>%
    select(required_columns) %>% 
    add_column(subID =  x) %>% 
    gather(marker,marker_value,reporters)
})

# Perform sampling n_samples random submissions n_iter number of times
# Score the combinations of the submission and keep track of the scores

repeated_scores <- tibble("Sample_size" = NA, "Iteration" = NA, "Mean" = NA, "Median" = NA)

for (n in 2:length(submissions)) {
  
  n_samples <- n
  n_iter<- 10
  
  for (j in 1:n_iter) {
    print(j)

    # Select the sampled predictions
    predictions_long <- bind_rows(sample(all_predictions, n_samples, replace = FALSE))
    
    # Get scores from the selected predictions
    mean_score <- predictions_long  %>%
      group_by(glob_cellID, cell_line, treatment, time, cellID, fileID, marker) %>%
      summarise(prediction = mean(marker_value)) %>%
      score_sc1_long_format(validation_data_long)
    
    median_score <- predictions_long  %>%
      group_by(glob_cellID, cell_line, treatment, time, cellID, fileID, marker) %>%
      summarise(prediction = median(marker_value)) %>%
      score_sc1_long_format(validation_data_long)
    
    repeated_scores <- repeated_scores %>% add_row("Sample_size" = n_samples, 
                                                   "Iteration" = j, 
                                                   "Mean" = mean_score, 
                                                   "Median" = median_score)
  }
  
}
repeated_scores <- filter(repeated_scores, !is.na(Sample_size))
#colnames(repeated_scores) <- c("Mean", "Median", "Weighted")
#repeated_scores <- repeated_scores %>% as_tibble()
if (FALSE) {saveRDS(repeated_scores, "prediction_combinations/SC1/SC1_random_subs_scores.rds")}

scores <- readRDS("prediction_combinations/SC1/SC1_random_subs_scores.rds")
ordered_combi <- readRDS("prediction_combinations/SC1/SC1_ordered_combination.rds") %>%
  gather(model, score, c(Mean, Median))

best_single <- leader_board  %>%
  rename(subID = objectId) %>%
  select(submitterId, score) %>%
  filter(score == min(score))

scores %>%
  gather(model, score, Mean, Median) %>%
  ggplot(aes(Sample_size, score, colour=model)) +
  geom_point()  +
  labs(title = "SC1", x= "Number of combined predictions")

scores %>%
  gather(model, score, Mean, Median) %>%
  ggplot(aes(as.factor(Sample_size), score, fill = model)) +
  geom_boxplot(position = "dodge") +
  xlim(as.character(1:22)) +
  geom_hline(aes(yintercept = best_single$score, colour = as.factor(best_single$submitterId))) +
  labs(title = "SC1", x= "Number of combined predictions", colour = "best single") +
  geom_point(data=ordered_combi, aes(x=n, y=score, shape="top n")) + 
  scale_shape_manual(values = 23) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))
  

scores %>%
  gather(model, score, Mean, Median) %>%
  group_by(Sample_size, model) %>%
  summarise(score = mean(score)) %>%
  ggplot(aes(Sample_size, score, colour=model)) +
  geom_point() +
  geom_smooth() +
  geom_hline(aes(yintercept = best_single$score, colour = as.factor(best_single$submitterId))) +
  labs(title = "SC1", x= "Number of combined predictions")  +
  geom_point(data=ordered_combi, aes(x=n, y=score, shape="top n"),  stroke =1) +
  scale_shape_manual(values = 3) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

scores %>%
  gather(model, score, Mean, Median) %>%
  arrange(score) 

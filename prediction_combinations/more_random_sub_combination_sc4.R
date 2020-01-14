## SC4
# Combine predictions from ramdomly selected submissions
# Combine different numbers of randomly selected submissions

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
source("./scoring_scripts/score_sc4_noFile.R")

submission_folder <- "./submission_data/final/SC4"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc4.csv"))

required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 

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
validation_data <- read_csv("~/Seafile/validation_data/sc4gold.csv") %>% select(required_columns)

submissions <- list.files(submission_folder) 
submissions <- submissions[!startsWith(submissions, "leaderboard")]

# Read all predictions
all_predictions <- lapply(submissions, function(x) {
  read_csv(file.path(submission_folder,x))  %>%
    select(required_columns) %>% 
    add_column(subID = as.numeric(sub(".csv", "", x)))
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
    predictions <- bind_rows(sample(all_predictions, n_samples, replace = FALSE))
    
    predictions_long <- predictions %>% 
      gather(marker,marker_value,-cell_line, -treatment, -time, -subID)
    
    # Get scores from the selected predictions 
    scores <- leader_board %>%
      filter(objectId %in% predictions$subID)  %>%
      rename(subID = objectId) %>%
      select(subID, score)
    
    mean_score <- predictions_long  %>%
      group_by(cell_line, treatment, time, marker) %>%
      summarise(prediction = mean(marker_value)) %>%
      score_sc4_long_format(validation_data)
    
    median_score <- predictions_long  %>%
      group_by(cell_line, treatment, time, marker) %>%
      summarise(prediction = median(marker_value)) %>%
      score_sc4_long_format(validation_data)
  
    all_scores <- scores %>% 
      rename(model = subID) %>%
      add_row(model = c("Mean", "Median"), 
              score = c(mean_score, median_score)) %>%
      select(model, score) %>%
      arrange(score)
    print(all_scores)
    
    repeated_scores <- repeated_scores %>% add_row("Sample_size" = n_samples, 
                                                   "Iteration" = j, 
                                                   "Mean" = mean_score, 
                                                   "Median" = median_score)
  }
  
}
repeated_scores <- filter(repeated_scores, !is.na(Sample_size))
#colnames(repeated_scores) <- c("Mean", "Median", "Weighted")
#repeated_scores <- repeated_scores %>% as_tibble()

if (FALSE) {saveRDS(repeated_scores, "prediction_combinations/SC4/SC4_random_subs_scores.rds")}

scores <- readRDS("prediction_combinations/SC4/SC4_random_subs_scores.rds")
ordered_combi <- readRDS("prediction_combinations/SC4/SC4_ordered_combination.rds") %>%
  gather(model, score, c(Mean, Median))

best_single <- leader_board  %>%
  rename(subID = objectId) %>%
  select(submitterId, score) %>%
  filter(score == min(score))
  
scores %>%
  gather(model, score, Mean, Median) %>%
  ggplot(aes(Sample_size, score, colour=model)) +
  geom_point()  +
  labs(title = "SC4", x= "Number of combined predictions")

scores %>%
  gather(model, score, Mean, Median) %>%
  ggplot(aes(as.factor(Sample_size), score, fill = model)) +
  geom_boxplot(position = "dodge") +
  xlim(as.character(1:20)) +
  geom_hline(aes(yintercept = best_single$score, colour = as.factor(best_single$submitterId))) +
  labs(title = "SC4", x= "Number of combined predictions", colour = "best single")  +
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
  labs(title = "SC4", x= "Number of combined predictions")  +
  geom_point(data=ordered_combi, aes(x=n, y=score, shape="top n"),  stroke =1) +
  scale_shape_manual(values = 3) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

scores %>%
  gather(model, score, Mean, Median) %>%
  arrange(score) 


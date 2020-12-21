## SC1
# Combine predictions from ramdomly selected submissions
# Combine different numbers of randomly selected submissions

library(tidyverse)
library(Biobase)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# The predictions of each team are a column, one column per team
# Thecolumn standard is the goldesn standard/true value
all_predictions <- readRDS("./submission_analysis/intermediate_data/sc1_all_NP_predictions.rds")
all_predictions = all_predictions %>% select(-31:-62)

leader_board <- read_csv("./submission_data/final_round/SC1/leaderboard_final_sc1.csv")
submissions <- readRDS("./submission_analysis/intermediate_data/sc1_ranked_teams.rds") %>% as.character()

# Radnomly sample 1 submission by sampling scores from the leaderboard
repeated_scores <- tibble("Sample_size" = rep(1, 100), "Iteration" = seq(1,100), 
                         "Median" = sample(leader_board$score, 100, replace = TRUE)) 

# Perform sampling n_samples random submissions n_iter number of times
# Score the combinations of the submission and keep track of the scores
for (n in 2:length(submissions)) {
  
  n_samples <- n
  n_iter<- 100
  combinations = list()
  
  for (j in 1:n_iter) {
    print(paste0("samples: :", n_samples, " iteration: ", j))

    # Select the sampled predictions
    selected_submissions <- sample(submissions, n_samples)

    # Combine the selected submissions by taking the median and score it
    combined_rmse_per_condition <-  all_predictions %>% 
      mutate(prediction = rowMedians(as.matrix(.[selected_submissions]))) %>%
      group_by(cell_line, treatment, time, marker) %>%
      summarise(RMSE = sqrt(sum((standard - prediction)^2) / n())) 
    
    
    combinations[[j]] = combined_rmse_per_condition %>% mutate(iter = j)
    median_score = combined_rmse_per_condition %>% pull(RMSE) %>%
        mean()
    
    repeated_scores <- repeated_scores %>% add_row("Sample_size" = n_samples, 
                                                   "Iteration" = j, 
                                                   "Median" = median_score)
  }
  
}

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

## SC2
# Combine predictions from ramdomly selected submissions
# Combine different numbers of randomly selected submissions

library(dplyr)
library(readr)
library(tidyr)
library(tibble)

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

repeated_scores <- tibble("Sample_size" = NA, "Iteration" = NA, "Equal_sample" = NA)
for (n in 2:length(submissions)) {
  n_samples <- n
  n_iter<- 25
  
  # Perform sampling n_samples random submissions n_iter number of times
  # Each combination of sumbmissions is combined and scored ten times
  # Because sampling can differ each time. All scores are saves
  
  for (j in 1:n_iter) {
    print(j)
    
    # Select the sampled predictions
    predictions <- sample(all_predictions, n_samples, replace = FALSE)
    
    # Get scores from the selected predictions 
    # scores <- leader_board %>%
    #   filter(objectId %in% predictions$subID)  %>%
    #   rename(subID = objectId) %>%
    #   select(subID, score)
    
    inner_loop_scores <- c()
    
    for (k in 1:3) {
      print(paste0("Iteration: ", j, "; inner loop: ", k, "; sample size: ", n))
      
      equal_sample_size <- lapply(predictions, function(x) {x %>% 
          group_by(cell_line, treatment, time) %>%  sample_frac(1/n_samples)}) %>%
        bind_rows() %>%
        select(required_columns) %>%
        score_sc2(validation_data)
      
      # all_scores <- scores  %>%
      #   rename(model = subID) %>%
      #   add_row(model = "Equal",
      #           score = equal_sample_size) %>%
      #   arrange(score)
      
      print(equal_sample_size)

      inner_loop_scores <- c(inner_loop_scores, equal_sample_size)
    }
    
    repeated_scores <- repeated_scores %>%
      add_row(Sample_size = n, 
              Iteration = j, 
              Equal_sample = mean(inner_loop_scores))
    
  }
  
}

repeated_scores <- filter(repeated_scores, !is.na(Sample_size))
#colnames(repeated_scores) <- c("Subsample_round", "Equal", "Weighted")
#repeated_scores <- repeated_scores %>% as_tibble()

if (FALSE) {saveRDS(repeated_scores, "prediction_combinations/SC2/SC2_random_subs_scores.rds")}
                                
scores <- readRDS("prediction_combinations/SC2/SC2_random_subs_scores.rds")
ordered_combi <- readRDS("prediction_combinations/SC2/SC2_ordered_combination.rds") %>%
  gather(model, score, Equal_sample)
best_single <- leader_board  %>%
  rename(subID = objectId) %>%
  select(submitterId, score) %>%
  filter(score == min(score))


scores %>%
  gather(model, score, Equal_sample) %>%
  ggplot(aes(Sample_size, score, colour=model)) +
  geom_point()  +
  labs(title = "SC2", x= "Number of combined predictions")

scores %>%
  gather(model, score, Equal_sample) %>%
  ggplot(aes(as.factor(Sample_size), score, fill = model)) +
  geom_boxplot(position = "dodge") +
  xlim(as.character(1:16)) +
  geom_hline(aes(yintercept = best_single$score, colour = as.factor(best_single$submitterId))) +
  labs(title = "SC2", x= "Number of combined predictions", colour = "best single") +
  geom_point(data=ordered_combi, aes(x=n, y=score, shape="top n")) + 
  scale_shape_manual(values = 23) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

scores %>%
  gather(model, score, Equal_sample) %>%
  group_by(Sample_size, model) %>%
  summarise(score = mean(score)) %>%
  ggplot(aes(Sample_size, score, colour=model)) +
  geom_point() +
  geom_smooth() +
  geom_hline(aes(yintercept = best_single$score, colour = as.factor(best_single$submitterId))) +
  labs(title = "SC2", x= "Number of combined predictions")  +
  geom_point(data=ordered_combi, aes(x=n, y=score, shape="top n"),  stroke =1) +
  scale_shape_manual(values = 3) +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

scores %>%
  gather(model, score, Equal_sample) %>%
  arrange(score) 



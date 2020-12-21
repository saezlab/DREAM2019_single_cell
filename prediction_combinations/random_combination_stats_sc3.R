## SC3; stats
# Combine predictions from different number of ramdomly selected submissions
# Combining by taking the median of predicted statistics to predict new statistics
# These new statistics are directly scored

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(Biobase)
# 
# setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

leader_board <- read_csv("./submission_data/final_round/SC3/leaderboard_final_sc3.csv")
submissions <-  readRDS("./submission_analysis/intermediate_data/sc3_ranked_teams.rds") %>%
    as.character()
# Statistics per condition of the standard and calculated from the predictions, including median EGF t=0 values per condition
sub_data_values <- readRDS("./submission_analysis/intermediate_data/sc3_stats_conditions.rds") %>%
    select(cell_line, treatment, time, stat_variable, standard, submissions)

# Radnomly sample 1 submission by sampling scores from the leaderboard
repeated_scores <- tibble("Sample_size" = rep(1, 100), "Iteration" = seq(1,100), 
                          "Median" = sample(leader_board$score, 100, replace = TRUE))

# Perform sampling n submission randomly, n_iter number of times
# Score the combinations of the submission and keep track of the scores
for (n in 2:dim(leader_board)[1]) {
    n_samples <- n
    n_iter<- 100
    
    for (j in 1:n_iter) {
        print(paste0("samples: :", n_samples, " iteration: ", j))
        
        # Select the sampled predictions
        selected_submissions <- sample(unique(submissions), n_samples, replace = FALSE)
        
        # Combine by sampling equally from the selected submissions and score it
        median_stats <- sub_data_values %>%
            mutate(prediction = rowMedians(as.matrix(.[selected_submissions]))) %>%
            select(-submissions) %>%
            group_by(cell_line, treatment, time) %>%
            summarise(SSQ = sum((standard - prediction)^2)) %>%
            ungroup() %>%
            pull(SSQ) %>%
            sum() %>%
            sqrt()
        
        # Update all scores
        repeated_scores <- repeated_scores %>%
            add_row(Sample_size = n, 
                    Iteration = j, 
                    Median = median_stats)
    }
    
}

if (TRUE) {
    saveRDS(repeated_scores, "prediction_combinations/SC3/SC3_random_subs_scores_stats.rds")
}

stats_combi <- readRDS("prediction_combinations/SC3/SC3_random_subs_scores_stats.rds") %>%
    rename(Median_stats = Median) %>%
    gather(model, score, Median_stats)

sample_combi <- readRDS("prediction_combinations/SC3/SC3_random_subs_scores.rds") %>%
    gather(model, score, Equal_sample)

boxplot <- stats_combi %>%
    bind_rows(sample_combi) %>%
    mutate(Sample_size = as.factor(Sample_size)) %>%
    ggplot(aes(Sample_size, score, fill = model)) +
    geom_boxplot(position = "dodge") +
    labs(title = "SC3", x= "Number of combined predictions", fill = "combination\nmethod") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15)) +
    theme_bw()
boxplot

violonplot <- stats_combi %>%
    bind_rows(sample_combi) %>%
    mutate(Sample_size = as.factor(Sample_size)) %>%
    ggplot(aes(Sample_size, score, fill = model)) +
    geom_violin(scale = "width") +
    labs(title = "SC3", x= "Number of combined predictions", fill = "combination\nmethod") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15)) +
    theme_bw()
violonplot

# Comparing two different ways of combing SC2 and SC3 predictions
# Firstly sample multipple times a certain fraction of samples from each of the predictions
# compared to a linear combination of the median and covariance in each combination

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(MASS)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell") 

source("./scoring_scripts/score_sc2_noFile.R")


submission_folder <- "./submission_data/final/SC2"

leaderboard <- read_csv(file.path(submission_folder, "leaderboard_final_sc2.csv"))

submissions <- list.files(submission_folder)
submissions <- submissions[!startsWith(submissions, "leaderboard")]
# Sample 10.000 cells from multivarite distribution 
sample_MVN <- function(means, covs){
  sample <- MASS::mvrnorm(10000, mu = means, Sigma =  covs, tol=0.005) %>%
    as_tibble()
  return(sample)
}
# Make incomplete covariance matrix symmetrical again
makeSymm <- function(x) {
  m <- x %>% separate(variable, into = c("V1", "V2"), sep = "-") %>%
    pivot_wider(names_from = V2, values_from = pred_stat_value) %>%
    column_to_rownames("V1") 
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}


required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 

# Read all predictions
all_predictions <- lapply(submissions, function(x) {
  read_csv(file.path(submission_folder,x))  %>%
    dplyr::select(required_columns) %>% 
    add_column(subID = as.numeric(sub(".csv", "", x)))
})

validation_data<- read_csv("~/Seafile/validation_data/sc2gold.csv") %>% dplyr::select(required_columns)
validation_stats <- validation_data %>%
  data_to_stats() %>% 
  dplyr::select(cell_line, treatment, time, stat_variable, stat_value) %>%
  rename(val_stat_value = stat_value)

overall_scores <- tibble("iteration" = NA,
                         "mean_stats" = NA, 
                         "median_stats" = NA, 
                         "weighted_stats" = NA, 
                         "equal_sample" = NA, 
                         "weighted_sample" = NA,
                         "median_sample" = NA,
                         "mean_sample" = NA)
for (k in 1:10) {
  print(k)
  n_samples <- 10
  
  # Select the sampled predictions
  predictions <- bind_rows(sample(all_predictions, n_samples, replace = FALSE))
  
  scores <- leaderboard %>%
    dplyr::filter(objectId %in% predictions$subID)  %>%
    rename(subID = objectId) %>%
    mutate(min_max_score = max(score) - score,
           weight = min_max_score / sum(min_max_score)) %>%
    dplyr::select(subID, score, weight)
  
  predicted_stats <- predictions %>% 
    group_by(subID) %>% 
    nest() %>% 
    mutate(stats = map(data, data_to_stats),
           stats = map(stats, function(x) { x %>% dplyr::select(cell_line, treatment, time, stat_variable, stat_value)})) %>%
    dplyr::select(subID, stats) %>% 
    unnest(cols = c(stats)) %>% 
    ungroup()
  
  mean_stats_score <- predicted_stats %>%
    group_by(cell_line, treatment, time, stat_variable) %>%
    summarise(pred_stat_value = mean(stat_value)) %>%
    left_join(validation_stats) %>%
    ungroup() %>%
    dplyr::select(val_stat_value, pred_stat_value) %>%
    as.matrix() %>%
    t() %>%
    dist(method = "euclidean")
  
  median_stats_score <- predicted_stats %>%
    group_by(cell_line, treatment, time, stat_variable) %>%
    summarise(pred_stat_value = median(stat_value)) %>%
    left_join(validation_stats) %>%
    ungroup() %>%
    dplyr::select(val_stat_value, pred_stat_value) %>%
    as.matrix() %>%
    t() %>%
    dist(method = "euclidean")
  
  weighted_stats_score <- predicted_stats %>%
    left_join(scores) %>%
    mutate(weighted_stat_value = weight*stat_value) %>%
    group_by(cell_line, treatment, time, stat_variable) %>%
    summarise(pred_stat_value = sum(weighted_stat_value)) %>%
    left_join(validation_stats) %>%
    ungroup() %>%
    dplyr::select(val_stat_value, pred_stat_value) %>%
    as.matrix() %>%
    t() %>%
    dist(method = "euclidean")

  equal_sample_size <- predictions %>% 
      group_by(cell_line, treatment, time, subID) %>%
      sample_frac(1/length(unique(predictions$subID))) %>%
      ungroup() %>%
      score_sc2(validation_data)
    
  sampled_from_median_stats <- predicted_stats %>%
      group_by(cell_line, treatment, time, stat_variable) %>%
      summarise(pred_stat_value = median(stat_value)) %>%
      separate(stat_variable, into = c("stat", "variable"),sep="_") %>%
      group_by(cell_line, treatment, time, stat) %>%
      nest(stat_values = c(variable, pred_stat_value)) %>%
      pivot_wider(names_from = stat, values_from = stat_values) %>%
      mutate(full_cov = map(cov, makeSymm),
             full_mean = map(mean, ~pull(., pred_stat_value))) %>%
      mutate(prediction = map2(full_mean, full_cov, sample_MVN)) %>%
      dplyr::select(cell_line, treatment, time, prediction) %>%
      unnest(prediction) %>%
      score_sc2(validation_data)
    
  sampled_from_mean_stats <- predicted_stats %>%
      group_by(cell_line, treatment, time, stat_variable) %>%
      summarise(pred_stat_value = mean(stat_value)) %>%
      separate(stat_variable, into = c("stat", "variable"),sep="_") %>%
      group_by(cell_line, treatment, time, stat) %>%
      nest(stat_values = c(variable, pred_stat_value)) %>%
      pivot_wider(names_from = stat, values_from = stat_values) %>%
      mutate(full_cov = map(cov, makeSymm),
             full_mean = map(mean, ~pull(., pred_stat_value))) %>%
      mutate(prediction = map2(full_mean, full_cov, sample_MVN)) %>%
      dplyr::select(cell_line, treatment, time, prediction) %>%
      unnest(prediction) %>%
      score_sc2(validation_data)
    
  weighted_by_score <- predictions %>%
      group_by(subID) %>%
      nest() %>%
      ungroup() %>%
      left_join(scores) %>%
      mutate(min_max_score = max(score) - score,
             weight = min_max_score / sum(min_max_score),
             subsample = map2(data, weight, function(x,y) {x %>% 
                 group_by(cell_line, treatment, time) %>% 
                 sample_frac(y)})) %>%
      dplyr::select(subsample) %>%
      unnest() %>%
      score_sc2(validation_data)

  overall_scores <- overall_scores %>% add_row("iteration" = k,
                                               "mean_stats" = mean_stats_score, 
                                               "median_stats" = median_stats_score, 
                                               "weighted_stats" = weighted_stats_score, 
                                               "equal_sample" = equal_sample_size, 
                                               "weighted_sample" = weighted_by_score,
                                               "median_sample" = sampled_from_median_stats,
                                               "mean_sample" = sampled_from_mean_stats)
  
}
overall_scores <- overall_scores %>% dplyr::filter(!is.na(iteration))

if (FALSE) {saveRDS(overall_scores, "./prediction_combinations/SC2/combi_methods_N10_scores.rds")}

overall_scores <- readRDS("./prediction_combinations/SC2/combi_methods_N10_scores.rds")

overall_scores %>%
  gather(method, score, c(mean_stats, median_stats,weighted_stats, equal_sample, weighted_sample)) %>%
  ggplot(aes(method, score, fill = method)) +
  theme_bw() +
  geom_boxplot(position = "dodge")
  
overall_scores %>%
  gather(method, score, c(mean_stats, median_stats,weighted_stats, equal_sample, weighted_sample)) %>%
  ggplot(aes(iteration, score, fill = method)) +
  theme_bw() +
  geom_bar(stat="identity", position = "dodge") +
  labs(fill = "combination\nmethod")


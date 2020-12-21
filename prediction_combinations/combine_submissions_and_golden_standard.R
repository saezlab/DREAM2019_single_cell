# SC1-SC4
# Combine all the prediction and the golden standard into 1 dataframe
# For Sc1 and SC4 the predicted values of each team and the standard will be one column
# For SC2 and SC3 the predicted values are nested and the dataframe will be team-cell_line-Tr-Time-dataframe with predictions
# and the statistics of the predictions will by as in SC1/4 with standard and a column per team
# By from Atilla Gabor on 10-12-2019, edited by Alice Driessen on 27-01-2020


library(tidyverse)
setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# SC1 --------------------------------------------------------------------------


# make sure to download the data of each team to this folder: 
# naming convention for files: `submissionID`.csv
submission_folder = "./submission_data/final/SC1/"

# read leaderboard
SC_leaderboard = read_csv(
  file.path(submission_folder,"leaderboard_final_sc1.csv")) %>%
  select(-writeUp, -createdOn) %>% 
  mutate(submissions = paste0(objectId,".csv")) %>% 
  arrange(score) %>% 
  mutate(submitterId = make.names(submitterId))

# Team KAUST_RSS didnt follow the template but managed to pass the validation 
SC_leaderboard <- SC_leaderboard %>% filter(submitterId != "KAUST_RSS") 

# get the team names ranked by performance
ranked_teams <- factor(SC_leaderboard$submitterId,
                       levels = SC_leaderboard$submitterId)
saveRDS(ranked_teams, "./submission_data/intermediate_data/sc1_ranked_teams.rds")

# read team's predictions from csv files
prediction_data <- SC_leaderboard  %>% 
  mutate(predictions = map(file.path(submission_folder,submissions),read_csv))

reporters <- c("p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")

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
ordered_predictions <- prediction_data %>%
  mutate(predictions = map(predictions,order_predictions)) %>%
  do(bind_cols(.$predictions))
names(ordered_predictions) <- prediction_data$submitterId

# read golden standard
gs <- read_csv("./challenge_data/validation_data/sc1gold.csv") %>%
  arrange(glob_cellID) %>%
  gather(key = "marker", value = "standard", reporters)

## combined data: 
# columns describe the conditions, 
# 	standard: golden standard measurement
combined_data <- bind_cols(gs,ordered_predictions)
saveRDS(combined_data, "./submission_data/intermediate_data/sc1_all_predictions.rds")

combined_median_data <- combined_data %>%
  group_by(cell_line, treatment, time, marker) %>%
  summarise_at(c("standard", as.character(ranked_teams)), median)
saveRDS(combined_median_data, "./submission_data/intermediate_data/sc1_median_conditions.rds")

# calculate the RMSE error for each conditions for each team. 
RMSE_conditions <- combined_data %>% 
  group_by(cell_line, treatment, time, marker) %>%
  summarise_at(as.character(ranked_teams),~ sqrt(sum((standard - .)^2) / n())) %>%
  ungroup()
saveRDS(RMSE_conditions, "./submission_data/intermediate_data/sc1_rmse_conditions.rds")



# --------------------------------- SC2 ---------------------------------------------


source("./scoring_scripts/score_sc2_noFile.R")
# make sure to download the data of each team to this folder: 
# naming convention for files: `submissionID`.csv
submission_folder = "./submission_data/final_round/SC2/"

# read leaderboard
SC_leaderboard = read_csv(file.path(submission_folder,"leaderboard_final_sc2.csv")) %>%
  select(-writeUp, -createdOn) %>% 
  mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
  mutate(submitterId = make.names(submitterId))

# get the team names ranked by performance
ranked_teams <- factor(SC_leaderboard$submitterId,levels = SC_leaderboard$submitterId)
#saveRDS(ranked_teams, "./submission_data/intermediate_data/sc2_ranked_teams.rds")

# read team's predictions from csv files and compute the stats
required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 
#' read_to_stats
#' 
#' reads data and compute stats --> we only work with the stats afterwards
read_to_stats <- function(file_name){
  read_csv(file_name) %>% select(required_columns) %>% data_to_stats(.) %>% ungroup()
}

prediction_data <- SC_leaderboard %>% 
  mutate(predictions = map(file.path(submission_folder,submissions),read_to_stats))


#' order_predictions
#'
#' takes a prediction matrix, reformats to long format and keeps only the values.
#' 
#' @param stat_matrix  matrix with the statistics in tibble format
#' @return value vector gathered over reporters and sorted by glob_cellID
order_stats <- function(stat_matrix){
  
  stat_matrix %>% arrange(cell_line, treatment,  time, stat_variable) %>%
    select(stat_value)
}

# order the predictions the same way for all teams and for the golden standard
# so a simple column bind will merge well.

ordered_predictions <- prediction_data %>% 
  mutate(predictions = map(predictions,order_stats)) %>%
  do(bind_cols(.$predictions))
names(ordered_predictions) <- prediction_data$submitterId

# read golden standard
gs <- read_to_stats("./challenge_data/validation_data/sc2gold.csv") %>%
  arrange(cell_line, treatment,  time, stat_variable) %>% 
  rename(standard=stat_value)

## combined data: 
# columns describe the conditions, 
# 	standard:  stats computed from the golden standard measurement

combined_statistics <- bind_cols(gs,ordered_predictions)  %>%
  unite("cond_id", cell_line, treatment, time,sep="_", remove = FALSE) %>% 
  select(cond_id, cell_line, treatment, time, stat_variable, standard, as.character(ranked_teams))
saveRDS(combined_statistics, "./submission_data/intermediate_data/sc2_stats_conditions.rds")


# calculate ahead the squared diff between the stats for each condition
condition_stats_SumSquared <- combined_statistics %>% 
  group_by(cond_id) %>%
  summarise_at(as.character(ranked_teams),~ sum((standard - .)^2))
saveRDS(condition_stats_SumSquared,"./submission_data/intermediate_data/sc2_stats_sumSquared_conditions.rds")

# get all predictions and save them as nested dataframce
all_predictions <- lapply(SC_leaderboard$objectId, function(x) {
  read_csv(file.path(submission_folder,paste0(x, ".csv")))  %>%
    select(required_columns)
})
names(all_predictions) <- prediction_data$submitterId

nested_predictions <- bind_rows(all_predictions, .id = "team") %>%
  group_by(team, cell_line, treatment, time) %>%
  arrange(team, cell_line, treatment, time) %>%
  nest() %>%
  ungroup() 
write_rds(nested_predictions, "./submission_analysis/intermediate_data/sc2_all_predictions_nested.rds")


# SC3 --------------------------------------------------------------------------
source("./scoring_scripts/score_sc3_noFile.R")

# make sure to download the data of each team to this folder: 
# naming convention for files: `submissionID`.csv
submission_folder = "./submission_data/final/SC3/"

# read leaderboard
SC_leaderboard = read_csv(file.path(submission_folder,"leaderboard_final_sc3.csv")) %>%
  select(-writeUp, -createdOn) %>% 
  mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
  mutate(submitterId = make.names(submitterId))


# get the team names ranked by performance
ranked_teams <- factor(SC_leaderboard$submitterId,levels = SC_leaderboard$submitterId)
saveRDS(ranked_teams, "./submission_data/intermediate_data/sc3_ranked_teams.rds")

# read team's predictions from csv files and compute the stats
required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 
#' read_to_stats
#' 
#' reads data and compute stats --> we only work with the stats afterwards
read_to_stats <- function(file_name){
  
  read_csv(file_name) %>% select(required_columns) %>% data_to_stats(.) %>% ungroup()
}

prediction_data <- SC_leaderboard %>% 
  mutate(predictions = map(file.path(submission_folder,submissions),read_to_stats))


#' order_predictions
#'
#' takes a prediction matrix, reformats to long format and keeps only the values.
#' 
#' @param stat_matrix  matrix with the statistics in tibble format
#' @return value vector gathered over reporters and sorted by glob_cellID
order_stats <- function(stat_matrix){
  
  stat_matrix %>% arrange(cell_line, treatment,  time, stat_variable) %>%
    select(stat_value)
}

# order the predictions the same way for all teams and for the golden standard
# so a simple column bind will merge well.

ordered_predictions <- prediction_data %>% mutate(predictions = map(predictions,order_stats)) %>%
  do(bind_cols(.$predictions))
names(ordered_predictions) <- prediction_data$submitterId

# read golden standard
gs <- read_to_stats("./challenge_data/validation_data/sc3gold.csv") %>%
  arrange(cell_line, treatment,  time, stat_variable) %>% 
  rename(standard=stat_value)

## combined data: 
# columns describe the conditions, 
# 	standard:  stats computed from the golden standard measurement

combined_statistics <- bind_cols(gs,ordered_predictions)  %>%
  unite("cond_id", cell_line, treatment, time,sep="_", remove = FALSE) %>% 
  select(cond_id, cell_line, treatment, time, stat_variable, standard, as.character(ranked_teams))
saveRDS(combined_statistics, "./submission_data/intermediate_data/sc3_stats_conditions.rds")

# calculate ahead the squared diff between the stats for each condition
condition_stats_SumSquared <- combined_statistics %>% 
  group_by(cond_id) %>% 
  summarise_at(as.character(ranked_teams),~ sum((standard - .)^2))
saveRDS(condition_stats_SumSquared, "./submission_data/intermediate_data/sc3_stats_sumSquared_conditions.rds")

# get all predictions and save them as nested dataframce
all_predictions <- lapply(SC_leaderboard$objectId, function(x) {
  read_csv(file.path(submission_folder,paste0(x, ".csv")))  %>%
    select(required_columns)  
})
names(all_predictions) <-prediction_data$submitterId

nested_predictions <- lapply(all_predictions, function(x){x  %>%
    group_by(cell_line, treatment, time) %>%
    arrange(cell_line, treatment, time) %>%
    nest() %>%
    ungroup()}) %>%
  bind_rows(.id = "team")

saveRDS(nested_predictions, "./submission_data/intermediate_data/sc3_all_predictions_nested.rds")



# SC4 --------------------------------------------------------------------------


# make sure to download the data of each team to this folder: 
# naming convention for files: `submissionID`.csv
submission_folder = "./submission_data/final/SC4/"

# read leaderboard
SC_leaderboard = read_csv(file.path(submission_folder,"leaderboard_final_sc4.csv")) %>%
  select(-writeUp, -createdOn) %>% 
  mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
  mutate(submitterId = make.names(submitterId))

# get the team names ranked by performance
ranked_teams <- factor(SC_leaderboard$submitterId,levels = SC_leaderboard$submitterId)
saveRDS(ranked_teams, "./submission_data/intermediate_data/sc4_ranked_teams.rds")

# read team's predictions from csv files
prediction_data <- SC_leaderboard  %>% 
  mutate(predictions = map(file.path(submission_folder,submissions),read_csv))

reporters <- c('b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
               'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
               'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
               'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
               'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
               'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
               'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
               'p.STAT3', 'p.STAT5') 

#' order_predictions
#'
#' takes a prediction matrix, reformats to long format and keeps only the values.
#' 
#' @param predictions prediction matrix in tibble format
#' @return value vector gathered over reporters and sorted by glob_cellID
order_predictions <- function(predictions){
  
  predictions %>%	arrange(cell_line, treatment,  time) %>%
    gather("reporter","value",reporters) %>% select(value)
}

# order the predictions the same way for all teams and for the golden standard
# so a simple column bind will merge well.

ordered_predictions <- prediction_data %>% 
  mutate(predictions = map(predictions,order_predictions)) %>%
  do(bind_cols(.$predictions))
names(ordered_predictions) <- prediction_data$submitterId

# read golden standard
gs <- read_csv("./challenge_data/validation_data/sc4gold.csv") %>%
  arrange(cell_line, treatment,  time) %>%
  gather(key = "marker", value = "standard", reporters)

## combined data: 
# columns describe the conditions, 
# 	standard: golden standard measurement

combined_data <- bind_cols(gs,ordered_predictions)
saveRDS(combined_data, "./submission_data/intermediate_data/sc4_combined_data.rds")

# calculate the RMSE error for each conditions for each team. 
RMSE_conditions <- combined_data %>% 
  group_by(cell_line, treatment, marker) %>%
  summarise_at(as.character(ranked_teams),~ sqrt(sum((standard - .)^2) / n())) %>%
  ungroup()
saveRDS(RMSE_conditions, "./submission_data/intermediate_data/sc4_rmse_conditions.rds")

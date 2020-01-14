# Combine all submssions to one file
library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

## -------------------------------- SC1 ---------------------------------
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
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc1.csv") )%>%
  select(-writeUp, -createdOn) %>% 
  mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
  mutate(submitterId = make.names(submitterId))
leader_board <- leader_board %>% filter(submitterId != "KAUST_RSS") 

validation_data_long <- read_csv("~/Seafile/validation_data/sc1gold.csv")  %>%
  select(required_columns) %>%
  arrange(glob_cellID) %>%
  gather(marker, standard, reporters)

submissions <- pull(leader_board, submitterId)

# read team's predictions from csv files
prediction_data <- leader_board  %>% 
  mutate(predictions = map(file.path(submission_folder,submissions),read_csv))

ordered_predictions <- prediction_data %>%
  mutate(predictions = map(predictions,order_predictions)) %>%
  do(bind_cols(.$predictions))
names(ordered_predictions) <- prediction_data$submitterId

combined_data <- bind_cols(validation_data_long, ordered_predictions)
saveRDS(combined_data, "./submission_data/intermediate_data/sc1_all_predictions.rds")

##-------------------------------------------------------- SC2
submission_folder <- "./submission_data/final/SC2"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc2.csv")) %>%
  arrange(score) %>%
  mutate(submitterId = sub("@", "X.", submitterId),
         submitterId = sub(" ", ".", submitterId))
required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 
all_predictions <- lapply(leader_board$objectId, function(x) {
  read_csv(file.path(submission_folder,paste0(x, ".csv")))  %>%
    select(required_columns)
})
names(all_predictions) <-leader_board$submitterId

nested_predictions <- bind_rows(all_predictions, .id = "team") %>%
  group_by(team, cell_line, treatment, time) %>%
  arrange(team, cell_line, treatment, time) %>%
  nest() %>%
  ungroup() %>%
  mutate(time=ifelse(time==16, 17, time))

saveRDS(nested_predictions, "./submission_data/intermediate_data/sc2_all_predictions_nested.rds")


## SC3
submission_folder <- "./submission_data/final/SC3"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc3.csv")) %>%
  arrange(score) %>%
  mutate(submitterId = sub("@", "X.", submitterId),
         submitterId = sub(" ", ".", submitterId))
required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 
all_predictions <- lapply(leader_board$objectId, function(x) {
  read_csv(file.path(submission_folder,paste0(x, ".csv")))  %>%
    select(required_columns)  
})
names(all_predictions) <-leader_board$submitterId

nested_predictions <- lapply(all_predictions, function(x){x  %>%
    group_by(cell_line, treatment, time) %>%
    arrange(cell_line, treatment, time) %>%
    nest() %>%
    ungroup() %>%
    mutate(time=ifelse(time==16, 17, time),
           time = ifelse(time==14, 13, time))}) %>%
  bind_rows(.id = "team")

saveRDS(nested_predictions, "./submission_data/intermediate_data/sc3_all_predictions_nested.rds")


## SC4
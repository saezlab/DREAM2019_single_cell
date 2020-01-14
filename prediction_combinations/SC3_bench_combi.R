# To compare combined prediction scores to bechmark, random, replicates for SC3
library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
source("./scoring_scripts/score_sc3_noFile.R")

validation_file <- "~/Seafile/validation_data/sc3gold.csv"
fileID <- read.csv("~/Desktop/BQ internship/BQ/DREAM/FileID_table.csv") %>% as_tibble()

required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 

validation_data <- read_csv(validation_file)

valA <- validation_data %>% left_join(select(fileID, fileID, time_course)) %>%
  filter(time_course == "A") %>%
  select(-time_course)
valB <- validation_data %>% left_join(select(fileID, fileID, time_course)) %>%
  filter(time_course == "B") %>%
  select(-time_course)
score_replicates <- score_sc3(valA, valB)

#better_and_best <- readRDS("prediction_combinations/SC3/SC23performance.rds")

single_model_scores <- readRDS("prediction_combinations/SC3/SC3_prediction_scores.rds") %>% select(prediction, score) %>% rename(model = prediction)
combi_model_scores <- readRDS("prediction_combinations/SC3/SC3_scores.rds") %>% gather(model, score)
random_scores <- readRDS("prediction_combinations/SC3/SC3_random_scores.rds")
#attila_score <- read_csv( "~/Seafile/dry_run/aim1_2_1_predictions.csv") %>% score_sc3(validation_data)


single_scores <- tibble(model = c("best_single", "replicatesT0"),
                        score = c(min(single_model_scores$score), score_replicates))

repeated_scores <- combi_model_scores %>% bind_rows(random_scores)

all_scores <- repeated_scores  %>%
  bind_rows(single_scores, single_model_scores) %>%
  filter(model != "best_single")

# PLot results
repeated_scores  %>%
  ggplot(aes(model, score, fill = model)) +
  geom_boxplot() +
  geom_hline(data = single_scores, mapping = aes(yintercept = score, color=model))  +
  scale_y_reverse() +
  labs(title = "SC3")

repeated_scores %>%
  ggplot(aes(score, fill = model)) +
  geom_density(kernel = "gaussian", alpha = 0.5)  +
  geom_segment(single_scores, mapping = aes(x = score, y = 0, xend = score, yend = Inf, color=model))  +
  labs(title = "SC3")


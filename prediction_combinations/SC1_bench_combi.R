# To compare combined prediction scores to bechmark, random, replicates for SC1
library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
source("./scoring_scripts/score_sc1_noFile.R")

required_columns <- c("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID", 
                      "p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")
reporters <- c("p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")

validation_data <- read_csv("~/Seafile/validation_data/sc1gold.csv") %>% select(required_columns)

#better_and_best <- readRDS("prediction_combinations/SC1/SC1_performance.rds")

single_model_scores <- readRDS("prediction_combinations/SC1/SC1_prediction_scores.rds") %>% select(prediction, score) %>% rename(model = prediction)
combi_model_scores <- readRDS("prediction_combinations/SC1/SC1_scores.rds") %>% gather(model, score)
random_scores <- readRDS("prediction_combinations/SC1/SC1_random_scores.rds")
attila_score <- read_csv( "prediction_combinations/SC1/SC1_predictions.csv") %>% score_sc1(validation_data)

single_scores <- tibble(model = c("Attila", "best_single"),
                       score = c(attila_score,
                                 min(single_model_scores$score))) %>%
  bind_rows(combi_model_scores)

repeated_scores <- random_scores

all_scores <- repeated_scores  %>%
  bind_rows(single_scores, single_model_scores) %>%
  filter(model != "best_single")

# PLot results
combi_model_scores %>%
  bind_rows(repeated_scores) %>%
  ggplot(aes(model, score, colour = model)) +
  geom_boxplot()  +
  geom_hline(data = filter(single_scores, !model %in% combi_model_scores$model), mapping = aes(yintercept = score, color=model))  +
  scale_y_reverse() +
  labs(title = "SC1")

repeated_scores %>%
  ggplot(aes(score, fill = model)) +
  geom_density(kernel = "gaussian", alpha = 0.5) +
  geom_segment(single_scores, mapping = aes(x = score, y = 0, xend = score, yend = Inf, color=model)) +
  labs(title = "SC1")



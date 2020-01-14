# To compare combined prediction scores to bechmark, random, replicates for SC1
library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Scoring function wtth input a tibble, not a csv file
source("./scoring_scripts/score_sc4_noFile.R")

required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 

validation_data <- read_csv("~/Seafile/validation_data/sc4gold.csv") %>% select(required_columns)

#better_and_best <- readRDS("prediction_combinations/SC1/SC1_performance.rds")

single_model_scores <- readRDS("prediction_combinations/SC4/SC4_prediction_scores.rds") %>% select(prediction, score) %>% rename(model = prediction)
combi_model_scores <- readRDS("prediction_combinations/SC4/SC4_scores.rds") %>% gather(model, score)
random_scores <- readRDS("prediction_combinations/SC4/SC4_random_scores.rds")
attila_score <- read_csv( "~/Seafile/dry_run/aim2_model_predictions_test.csv") %>% score_sc4(validation_data)

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
  ggplot(aes(model, score, fill = model, colour=model)) +
  geom_boxplot()  +
  geom_hline(data = filter(single_scores, !model %in% combi_model_scores$model), mapping = aes(yintercept = score, color=model))  +
  scale_y_reverse() +
  labs(title = "SC4")

repeated_scores %>%
  ggplot(aes(score, fill = model)) +
  geom_density(kernel = "gaussian", alpha = 0.5) +
  geom_segment(single_scores, mapping = aes(x = score, y = 0, xend = score, yend = Inf, color=model))  +
  labs(title = "SC4")



# Plot the scores of the random combinations, ordered combinations and cross validated results

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

library(tidyverse)

subchallenge <- "SC1"
challenge_folder <- file.path("prediction_combinations", subchallenge)

if (subchallenge == "SC1") {
# ----------------------------------- SC1 ---------------------------------------------------
CV_results <- readRDS(file.path(challenge_folder, "LOO_CV_RF_scores.rds"))
CV_scores <- CV_results %>% 
  select(-CV_loop, -val_CL) %>%
  colMeans() %>% 
  enframe() %>%
  rename(model = name, score=value) 
random_methods <- c("Mean", "Median")
labels_random <- c("random mean", "random median")
} else if (subchallenge == "SC2") {# ------------------ SC2 ----------------
CV_results <- readRDS(file.path(challenge_folder, "L2O_CV_scores.rds"))
CV_scores <- CV_results %>% 
  select(-CV_loop, -val_CL) %>% 
  apply(c(1,2), function(x){x^2}) %>%
  colSums() %>%
  sqrt() %>%
  enframe() %>%
  rename(model = name, score=value)
random_methods <- "Equal_sample"
labels_random <- "random"
 } else if(subchallenge == "SC3") { # --------------- SC3 --------------
   CV_results <- readRDS(file.path(challenge_folder, "L4O_CV_scores.rds"))
   CV_scores <- CV_results %>% 
     select(-CV_loop, -val_CL) %>% 
     apply(c(1,2), function(x){x^2}) %>%
     colSums() %>%
     sqrt() %>%
     enframe() %>%
     rename(model = name, score=value) 
   random_methods <- "Equal_sample"
 } else if (subchallenge == "SC4") {# --------------- SC4 -------------------
  CV_results <- readRDS(file.path(challenge_folder, "LOO_CV_scores.rds"))
  CV_scores <- CV_results %>% 
    select(-CV_loop, -val_CL) %>%
    colMeans() %>% 
    enframe() %>%
    rename(model = name, score=value)
  random_methods <- c("Mean", "Median")
 }

# Find single best performer of subchallenge based on the leaderboard
# Combine with CV scores into single scores tibble
leader_board <- read_csv(paste0("./submission_data/final/", subchallenge, "/leaderboard_final_", subchallenge, ".csv"))
single_scores <- leader_board  %>%
  select(submitterId, score) %>%
  filter(score == min(score)) %>%
  rename(model = submitterId) %>%
  bind_rows(CV_scores)

# The scores when combining the top n performing teams, biased/gold standard
ordered_scores <- readRDS(file.path(challenge_folder, paste0(subchallenge, "_ordered_combination.rds"))) %>%
  gather(model, score, random_methods) %>%
  rename(topn_model = model)
# Scores after combing random teams 100 times
random_scores <- readRDS(file.path(challenge_folder, paste0(subchallenge, "_random_subs_scores.rds"))) %>%
  gather(model, score, random_methods)

# Plot random scores as boxplots, top n scores as diamonds and all other scores as horizontal lines
if (subchallenge %in% c("SC2", "SC3")) {
  box_plot <- random_scores  %>%
    ggplot(aes(as.factor(Sample_size), score, fill = model)) +
    geom_boxplot(position = "dodge", aes(fill ="random")) +
    xlim(as.character(1:dim(leader_board)[1])) +
    geom_hline(data = single_scores, mapping = aes(yintercept = score, colour = model))   +
    geom_point(inherit.aes = FALSE, data=ordered_scores, aes(x=n, y=score, shape = "top n", fill = "top n")) +
    scale_shape_manual(values = 23) +
    labs(title = subchallenge, x= "Number of combined predictions") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15),
          legend.text = element_text(size=13)) +
    guides(fill = guide_legend(title= "combinations", override.aes = list(linetype = 0, shape = "",size = 4)),
           shape = guide_legend(title = "top n",override.aes = list(size = 4, fill = "#00BFC4")),
           colour = guide_legend(title = "model"))
  
  
  #  Plot random scores as line (geom_smooth), top n scores as plusses and all other scores as horizontal lines
  line_plot <- random_scores  %>%
    group_by(Sample_size, model) %>%
    summarise(score = mean(score)) %>%
    ggplot(aes(Sample_size, score, colour=model)) +
    geom_point(aes(colour = "random")) +
    geom_smooth(aes(colour = "random"))  +
    geom_hline(data = single_scores, mapping = aes(yintercept = score, colour = model))   +
    geom_point(inherit.aes = FALSE, data=ordered_scores, aes(x=n, y=score, shape="top n"), colour = "#00BFC4", stroke =1) +
    scale_shape_manual(values = 3) +
    xlim(as.character(1:dim(leader_board)[1])) +
    labs(title = subchallenge, x= "Number of combined predictions", shape = " ") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15),
          legend.text = element_text(size=13))
} else if (subchallenge %in% c("SC1", "SC4")) {
  box_plot <- random_scores  %>%
    ggplot(aes(as.factor(Sample_size), score, fill = model)) +
    geom_boxplot(position = "dodge") +
    xlim(as.character(1:dim(leader_board)[1])) +
    geom_hline(data = single_scores, mapping = aes(yintercept = score, colour = model))   +
    geom_point(inherit.aes = FALSE, data=ordered_scores, aes(x=n, y=score, shape = topn_model, fill = topn_model)) +
    scale_shape_manual(values = c(21, 24)) +
    labs(title = subchallenge, x= "Number of combined predictions") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15),
          legend.text = element_text(size=13)) +
    guides(fill = guide_legend(title= "random"),
           shape = guide_legend(title = "top n",override.aes = list(size = 4, fill = c("#F8766D", "#00BFC4"))),
           colour = guide_legend(title = "model"))
  
  
  #  Plot random scores as line (geom_smooth), top n scores as plusses and all other scores as horizontal lines
  line_plot <- random_scores  %>%
    mutate(model = paste0("random ", tolower(model))) %>%
    group_by(Sample_size, model) %>%
    summarise(score = mean(score)) %>%
    ggplot(aes(Sample_size, score, colour=model)) +
    geom_point() +
    geom_smooth()  +
    geom_hline(data = single_scores, mapping = aes(yintercept = score, colour = model))   +
    geom_point(inherit.aes = FALSE, data=ordered_scores, aes(x=n, y=score, shape=topn_model, fill = topn_model)) +
    scale_shape_manual(values = c(21, 24)) +
    xlim(as.character(1:dim(leader_board)[1])) +
    labs(title = subchallenge, x= "Number of combined predictions") +
    guides(fill = guide_legend(title = "top n"), shape = guide_legend(title = "top n")) + 
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15),
          legend.text = element_text(size=13))
}
box_plot


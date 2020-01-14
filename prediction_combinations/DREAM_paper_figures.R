# Make plots for the DREAM challenge paper for SC1-SC4

# Plot the scores of the random combinations, ordered combinations and cross validated results

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

library(tidyverse)

subchallenge <- "SC3"
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
  mutate(model = "challenge winner") %>%
  bind_rows(CV_scores) %>%
  filter(model !="SC_MDT") %>%
  mutate(model = case_when(model == "MDT" ~ "cond. best team",
                           model == "arbiter" ~ "cond. err. predictor (best)",
                           model == "lm" ~ "S.C. linear model",
                           model == "lc_arbiter" ~ "cond. err. predictor (weighted)",
                           model == "single" ~ "predicted winner",
                           model == "SC_MDT" ~ "S.C. best team",
                           TRUE ~ model)) %>%
  arrange(score)
single_scores$model <- factor(single_scores$model, levels = c(single_scores$model))

# The scores when combining the top n performing teams, biased/gold standard
ordered_scores <- readRDS(file.path(challenge_folder, paste0(subchallenge, "_ordered_combination.rds"))) %>%
  gather(model, score, random_methods) %>%
  rename(topn_model = model) %>%
  filter(topn_model != "Mean")
# Scores after combing random teams 100 times
random_scores <- readRDS(file.path(challenge_folder, paste0(subchallenge, "_random_subs_scores.rds"))) %>%
  gather(model, score, random_methods) %>%
  filter(model != "Mean")

# Plot random scores as boxplots, top n scores as diamonds and all other scores as horizontal lines
if (subchallenge == "SC2") {
  box_plot <- random_scores  %>%
    ggplot(aes(as.factor(Sample_size), score, fill = model)) +
    geom_boxplot(position = "dodge") +
    xlim(as.character(1:dim(leader_board)[1])) +
    ylim(22, 100) +
    geom_hline(data = single_scores, mapping = aes(yintercept = score, colour = model))   +
    geom_point(data=ordered_scores, aes(x=n, y=score, shape = "equal sample", fill = topn_model), size = 2.5) +
    scale_shape_manual(values = 23) +
    scale_fill_discrete(labels = "equal sample") +
    labs(title = subchallenge, x= "Number of combined predictions") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15),
          legend.text = element_text(size=13)) +
    guides(fill = guide_legend(title= "random", override.aes = list(shape = "", labels ="equal sample")),
           shape = guide_legend(title = "top n",override.aes = list(size = 4, fill = "#F8766D")),
           colour = guide_legend(title = "model"))

} else if (subchallenge == "SC3") {  
  box_plot <- random_scores  %>%
    ggplot(aes(as.factor(Sample_size), score, fill = model)) +
    geom_boxplot(position = "dodge") +
    xlim(as.character(1:dim(leader_board)[1])) +
    geom_hline(data = single_scores, mapping = aes(yintercept = score, colour = model))   +
    geom_point(data=ordered_scores, aes(x=n, y=score, shape = "equal sample", fill = topn_model), size = 2.5) +
    scale_shape_manual(values = 23) +
    scale_fill_discrete(labels = "equal sample") +
    labs(title = subchallenge, x= "Number of combined predictions") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15),
          legend.text = element_text(size=13)) +
    guides(fill = guide_legend(title= "combinations", override.aes = list(shape = "")),
           shape = guide_legend(title = "top n",override.aes = list(size = 4, fill = "#F8766D")),
           colour = guide_legend(title = "model"))

} else if (subchallenge  == "SC1") {
  box_plot <- random_scores  %>%
    ggplot(aes(as.factor(Sample_size), score, fill = model)) +
    geom_boxplot(position = "dodge") +
    xlim(as.character(1:dim(leader_board)[1])) +
    ylim(NA,1)+
    geom_hline(data = single_scores, scores, mapping = aes(yintercept = score, colour = model))   +
    geom_point(inherit.aes = FALSE, data=ordered_scores, aes(x=n, y=score, shape = topn_model, fill = topn_model)) +
    scale_shape_manual(values = 23) +
    labs(title = subchallenge, x= "Number of combined predictions") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15),
          legend.text = element_text(size=13)) +
    guides(fill = guide_legend(title= "random"),
           shape = guide_legend(title = "top n",override.aes = list(size = 4, fill = "#F8766D")),
           colour = guide_legend(title = "model"))
} else if (subchallenge  == "SC4") {
  box_plot <- random_scores  %>%
    ggplot(aes(as.factor(Sample_size), score, fill = model)) +
    geom_boxplot(position = "dodge") +
    xlim(as.character(1:dim(leader_board)[1])) +
    ylim(NA, 0.47) +
    geom_hline(data = single_scores, scores, mapping = aes(yintercept = score, colour = model))   +
    geom_point(inherit.aes = FALSE, data=ordered_scores, aes(x=n, y=score, shape = topn_model, fill = topn_model)) +
    scale_shape_manual(values = 23) +
    labs(title = subchallenge, x= "Number of combined predictions") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15),
          legend.text = element_text(size=13)) +
    guides(fill = guide_legend(title= "random"),
           shape = guide_legend(title = "top n",override.aes = list(size = 4, fill = "#F8766D")),
           colour = guide_legend(title = "model"))
}
box_plot


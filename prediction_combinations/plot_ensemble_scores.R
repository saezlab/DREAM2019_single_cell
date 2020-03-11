# Make plots for the DREAM challenge paper for SC1-SC4

# Plot the scores of the random combinations, ordered combinations and cross validated results

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

library(tidyverse)
library(wesanderson)

subchallenge <- "SC1"
challenge_folder <- file.path("prediction_combinations", subchallenge)

if (subchallenge == "SC1") {
  # ----------------------------------- SC1 ---------------------------------------------------
  CV_results1 <- readRDS(file.path(challenge_folder, "LOO_CV_RF_scores.rds"))
  CV_results2 <- readRDS(file.path(challenge_folder, "LOO_CV_DREAM_scores.rds"))
  CV_results3 <- readRDS(file.path(challenge_folder,"LOO_CV_asvRF_scores.rds"))
  CV_results4 <- readRDS(file.path(challenge_folder,"LOO_CV_asvRF_incl32_scores.rds"))
  CV_results5 <- readRDS(file.path(challenge_folder,"CV_asvRF_incl32_scores.rds")) 
  CV_results6 <- readRDS(file.path(challenge_folder,"Tr_CV_asvRF_incl32_scores.rds"))
  CV_results <- left_join(CV_results1, CV_results2) %>% 
    left_join(CV_results3) %>%
    left_join(CV_results4) %>%
    left_join(CV_results5) %>%
    left_join(CV_results6)
  CV_scores <- CV_results %>% 
    select(-CV_loop, -val_CL) %>%
    colMeans() %>% 
    enframe() %>%
    rename(model = name, score=value) %>%
    filter(!model %in% c("SC_MDT", "RF", "asv_RF", "ttm_MDT"))
  random_methods <- "Median"
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
  CV_results <- readRDS(file.path(challenge_folder, "LOO_CV_RF_scores.rds"))
  CV_scores <- CV_results %>% 
    select(-CV_loop, -val_CL) %>%
    colMeans() %>% 
    enframe() %>%
    rename(model = name, score=value)
  random_methods <- "Median"
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
  arrange(score) %>%
  mutate(model = case_when(model == "single" ~ "predicted winner",
                           TRUE ~ model)) %>%
  arrange(score)
single_scores$model <- factor(single_scores$model, levels = c("arbiter", "MDT", "lc_arbiter", "lm", "challenge winner",
                                                              "predicted winner","stats_lc_arbiter", "stats_RF", 
                                                              "cond_lm", "cond_RF", "asv_RF_32", "asv_RF_32_Tr", "CV_asv_RF_32"))
col <- wes_palette(11, name = "Darjeeling1", type = "continuous")
names(col) <-  c("arbiter", "MDT", "lc_arbiter", "lm", "challenge winner",
                 "predicted winner","stats_lc_arbiter", "stats_RF", 
                 "cond_lm", "cond_RF", "asv_RF_32", "asv_RF_32_Tr", "CV_asv_RF_32")
# The scores when combining the top n performing teams, biased/gold standard
ordered_scores <- readRDS(file.path(challenge_folder, paste0(subchallenge, "_ordered_combination.rds"))) %>%
  gather(model, score, random_methods) %>%
  rename(topn_model = model)

# Scores after combing random teams 100 times
random_scores <- readRDS(file.path(challenge_folder, paste0(subchallenge, "_random_subs_scores.rds"))) %>%
  gather(model, score, random_methods)

# Plot random scores as boxplots, top n scores as diamonds and all other scores as horizontal lines
if (subchallenge == "SC2") {
  box_plot <- random_scores  %>%
    ggplot(aes(as.factor(Sample_size), score, fill = model)) +
    theme_bw() +
    coord_cartesian(clip = "off") +
    geom_boxplot(position = "dodge", show.legend = FALSE) +
    geom_point(inherit.aes=FALSE, data = leader_board, mapping=aes(x=1, y=score), colour = "#00BFC4", show.legend = FALSE)+
    geom_hline(data = single_scores, mapping = aes(yintercept = score, colour = model), show.legend = FALSE) + 
    ggrepel::geom_text_repel(inherit.aes = FALSE, data = single_scores, aes(label = model, x = 17, y =score, colour = model), 
                             direction = "y", 
                             hjust = 0, 
                             segment.size = 0.2,
                             xlim = c(18, 22),
                             ylim = c(20, 80))  +
    scale_color_manual(values= col, guide = "none") +
    geom_point(data=ordered_scores, aes(x=n, y=score, shape = "equal sample", fill = topn_model), size = 2.5, show.legend = FALSE) +
    scale_shape_manual(values = 23) +
    labs(title = subchallenge, x= "Number of combined predictions") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15),
          plot.margin = unit(c(0.1,4.1,0.1,0.1), "cm")) 

} else if (subchallenge == "SC3") {  
  box_plot <- random_scores  %>%
    ggplot(aes(as.factor(Sample_size), score, fill = model)) +
    theme_bw() +
    coord_cartesian(clip = "off") +
    geom_boxplot(position = "dodge", show.legend = FALSE)+
    geom_point(inherit.aes=FALSE, data = leader_board, mapping=aes(x=1, y=score), colour = "#00BFC4", show.legend = FALSE) +
    geom_hline(data = single_scores, mapping = aes(yintercept = score, colour = model), show.legend = FALSE) + 
    ggrepel::geom_text_repel(inherit.aes = FALSE, data = single_scores, aes(label = model, x = 15, y =score, colour = model), 
                             direction = "y", 
                             hjust = 0, 
                             segment.size = 0.2,
                             xlim = c(16, 20),
                             ylim = c(50, 150))  +
    scale_color_manual(values= col, guide = "none") +
    geom_point(data=ordered_scores, aes(x=n, y=score, shape = "equal sample", fill = topn_model), size = 2.5, show.legend = FALSE) +
    scale_shape_manual(values = 23) +
    scale_fill_discrete(labels = "equal sample") +
    labs(title = subchallenge, x= "Number of combined predictions") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15),
          plot.margin = unit(c(0.1,4.2,0.1,0.1), "cm"))

} else if (subchallenge  == "SC1") {
  box_plot <- random_scores  %>%
    ggplot(aes(as.factor(Sample_size), score, fill = model)) +
    theme_bw() +
    coord_cartesian(clip = "off") +
    geom_boxplot(position = "dodge", show.legend = FALSE) +
    geom_point(inherit.aes=FALSE, data = leader_board, mapping=aes(x=1, y=score), colour = "#00BFC4", show.legend = FALSE)+
    ylim(NA,1) +
    geom_hline(data = single_scores, mapping = aes(yintercept = score, colour = model), show.legend = FALSE)  +
    ggrepel::geom_text_repel(inherit.aes = FALSE, data = single_scores, aes(label = model, x = 23, y =score, colour = model), 
                             direction = "y", 
                             hjust = 0, 
                             segment.size = 0.2,
                             na.rm = TRUE,
                             xlim = c(24, 28),
                             ylim = c(0, 1))  +
    scale_color_manual(values= col, guide = "none") +
    geom_point(inherit.aes = FALSE, data=ordered_scores, aes(x=n, y=score, shape = topn_model, fill = topn_model), show.legend = FALSE) +
    scale_shape_manual(values = 23) +
    labs(title = subchallenge, x= "Number of combined predictions") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15),
          plot.margin = unit(c(0.1,4,0.1,0.1), "cm"))
} else if (subchallenge  == "SC4") {
  box_plot <- random_scores  %>%
    ggplot(aes(as.factor(Sample_size), score, fill = model)) +
    theme_bw() +
    coord_cartesian(clip = "off") +
    geom_boxplot(position = "dodge", show.legend = FALSE)+
    geom_point(inherit.aes=FALSE, data = leader_board, mapping=aes(x=1, y=score), colour = "#00BFC4", show.legend = FALSE) +
    ylim(NA, 0.47) +
    geom_hline(data = single_scores, mapping = aes(yintercept = score, colour = model), show.legend = FALSE) + 
    ggrepel::geom_text_repel(inherit.aes = FALSE, data = single_scores, aes(label = model, x = 21, y =score, colour = model), 
                             direction = "y", 
                             hjust = 0, 
                             segment.size = 0.2,
                             xlim = c(22, 27),
                             ylim = c(0.25, 0.45)) +
    scale_color_manual(values=col, guide = "none") +
    geom_point(inherit.aes = FALSE, data=ordered_scores, aes(x=n, y=score, shape = topn_model, fill = topn_model), show.legend = FALSE) +
    scale_shape_manual(values = 23) +
    labs(title = subchallenge, x= "Number of combined predictions") +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          title = element_text(size=15),
          plot.margin = unit(c(0.1,4,0.1,0.1), "cm"))
}
box_plot
#ggsave(paste0(subchallenge, "_scores.pdf"), box_plot, device="pdf", path = "~/Desktop/", width =8.04 , height = 7.67, units="in")

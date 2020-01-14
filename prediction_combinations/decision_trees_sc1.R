## SC1
# Do crossvalidation to find how to combine the predictions, leaving two cell lines out as validation set
# A meta decision tree, arbiter and linear model are built
# Predictions made by selecting value of best submission as predicted by MDT, value taken from submission with lowest
# error predicted by the arbiter, linear combination of all submission basedon error predicted by arbiter and a l
# linear combination of of all submission by linear model that was trained on predicted and true values. 

library(tidyverse)
library(rpart)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# Participating team names
submissions <-  readRDS("./submission_data/intermediate_data/sc1_ranked_teams.rds") %>%
  as.character()
# Median prediction per team per condition (cell line, treatment, time, marker) and the median values of 32 non predicted markers\)
sub_data_median <- readRDS("./submission_data/intermediate_data/sc1_median_conditions_np.rds")
# Error per team per condition (cell line, treatment, time, marker) and the median values of 32 non predicted markers
sub_data_err <- readRDS( "./submission_data/intermediate_data/sc1_condErr_medianVals.rds") 
# All single cell predictions of all teams and true values (standard) as well as 32 marker values of non predicted markers
sub_data_all <- readRDS("./submission_data/intermediate_data/sc1_all_NP_predictions.rds")

np_markers <- c("b.CATENIN", "cleavedCas", "CyclinB", "GAPDH", "IdU", "Ki.67", "p.4EBP1", 
                "p.AKT.Thr308.", "p.AMPK", "p.BTK", "p.CREB", "p.FAK", "p.GSK3b", "p.H3", "p.JNK",
                "p.MAP2K3", "p.MAPKAPK2", "p.MEK", "p.MKK3.MKK6", "p.MKK4", "p.NFkB", "p.p38", "p.p53",
                "p.p90RSK", "p.PDPK1", "p.RB", "p.S6K", "p.SMAD23", "p.SRC", "p.STAT1", "p.STAT3", 
                "p.STAT5")
if (FALSE) {
# Each row represents a unique combination of two cell lines
# Each row will be used as validation cell lines once
CL_combi <- expand.grid(cell_lines, cell_lines) %>% 
  as_tibble() %>%
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  filter(! Var1 == Var2) %>%
  group_by(grp = paste(pmax(Var1, Var2), pmin(Var1, Var2), sep = "_")) %>%
  slice(1) %>%
  ungroup() %>%
  select(-grp)}

cell_lines <- unique(sub_data_err$cell_line)

# To keep count how often features are selected by MDT and the arbiter
feature_count <- tibble("feature" = c("treatment", "time", "marker", np_markers), 
                        "count" = 0)
# Keep track of scores per CV loop
scores <- tibble("CV_loop" = NA,
                 "val_CL" = NA,
                 "MDT" = NA, 
                 "arbiter" = NA,
                 "lc_arbiter" = NA,
                 "lm" = NA,
                 "single" = NA,
                 "SC_MDT" = NA)
# Save predictions made by each model
all_predictions <- tibble()
for (i in 1:length(cell_lines)) {
  
  print(paste("Iteration ", i, " out of ", length(cell_lines), ".", sep = ""))
  
  # select training and validation data
  train_data_err <- sub_data_err %>% # Used for arbiter, MDT and best single
    filter(!cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
  train_data_median <- sub_data_median %>% # Used for lm
    filter(!cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)

  val_data_err <- sub_data_err %>%
    filter(cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
  val_data_median <- sub_data_median %>%
    filter(cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
  val_data_all <- sub_data_all %>%
    filter(cell_line %in% cell_lines[i]) %>%
    arrange(cell_line, treatment, time, marker)
    
  ## Train models
  # formula for condition + 32 markers
  features <- paste(feature_count$feature, collapse = "+")
  # Meta Decision Tree
  if (FALSE) {
  MDT_train <- train_data_err %>% 
    group_by(treatment, time, marker, best_sub) %>% 
    summarise(n=n()) %>%
    arrange(treatment, time, marker, n) %>% 
    top_n(1) %>%
    group_by(treatment, time, marker)  %>%
    mutate(n=n()) %>%
    nest(best_sub)  %>%
    mutate(sel_sub = map(data, ~.[sample(1:n, 1),])) %>% 
    select(-c(n, data)) %>% 
    unnest() } # If there is no distinguishing cell lines
  MDT <- rpart(paste0("best_sub ~", features), data=train_data_err, method="class")
  
  #Arbiter
  arbiter <- lapply(submissions, function(x)
    {rpart(paste0(x, "~", features), data=train_data_err, method="anova")})
  names(arbiter) <- submissions
  
  # Linear model
  lm_formula <- paste0("standard ~", paste(submissions, collapse = "+"))
  linModel <- lm(lm_formula, data=train_data_median)
  
  # Select best single model
  single_model_scores <- train_data_err  %>%
    select(submissions) %>%
    colMeans()
  single_best <- names(single_model_scores)[which.min(single_model_scores)]
  
    
  ## ------------- Assess model performanceson condition level -------------------------------
  # Select for each condition the MDT estimated best predictor
  # Use error of the best predictors over the entire condition to score
  MDT_pred <-  val_data_err %>% 
    add_column(pred = predict(MDT, val_data_err, type = "class")) %>%
    select(-c(np_markers, best_sub)) %>%
    gather(team, MDT_error, submissions) %>%
    filter(pred==team) %>%
    select(cell_line, treatment, time, marker, MDT_error)
  
  MDT_val_score <- MDT_pred %>%
    pull(MDT_error) %>%
    mean()

  # Predict error of each submission on the validation data with the arbiter
  # Error per condition, not per single cell
  arbiter_pred <- lapply(arbiter, FUN = function(a) {predict(a, val_data_err)})   %>% as_tibble()

  # Select submission with lowest predicted error and score the combination
  arbiter_pred_value <- val_data_err %>%
    add_column(pred = names(arbiter_pred)[max.col(-arbiter_pred)]) %>%
    select(-c(np_markers, best_sub)) %>%
    gather(team, arbiter_error, submissions) %>%
    filter(pred==team)  %>%
    select(cell_line, treatment, time, marker, arbiter_error)
  
  arbiter_val_score <- arbiter_pred_value  %>%
    pull(arbiter_error) %>%
    mean()
  
  ## ------------------ Assess models on single cell level ---------------------------------
  # Linear combination of prediction based on predicted error by the arbiter W is the weights
  SC_arbiter_pred <- lapply(arbiter, FUN = function(a) {predict(a, val_data_all)})   %>% as_tibble()
  W <- lapply(SC_arbiter_pred, function(v) {(1/rowSums(SC_arbiter_pred)/v)}) %>% as_tibble()
  W <- lapply(W, function(x){x/rowSums(W)}) %>% as_tibble()
  lc_arbiter_pred <- val_data_all %>% 
    add_column(prediction = rowSums(select(val_data_all, submissions) * W)) %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(lc_arbiter_error = sqrt(sum((standard - prediction)^2) / n())) %>%
    ungroup() %>%
    select(cell_line, treatment, time, marker, lc_arbiter_error)
  
  lc_arbiter_val_score <- lc_arbiter_pred %>%
    pull(lc_arbiter_error) %>%
    mean()
  
  # Linear model
  lm_pred <- val_data_all %>% 
    add_column(prediction = predict(linModel, val_data_all)) %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(lm_error = sqrt(sum((standard - prediction)^2) / n())) %>%
    select(cell_line, treatment, time, marker, lm_error)
    
  lm_val_score <- lm_pred %>%
    pull(lm_error) %>%
    mean()
  
  # Score predictions of single best model
  SB_pred <- val_data_err %>% 
    select(cell_line, treatment, time, marker, single_best) %>%
    rename(SB_error = single_best)
  
  SB_val_score <- SB_pred %>%
    pull(SB_error) %>%
    mean()
  
  # MDT on single cells
  SC_MDT_pred <-  val_data_all %>% 
    add_column(pred = predict(MDT, val_data_all, type = "class")) %>%
    select(-np_markers) %>%
    gather(team, prediction, submissions) %>%
    filter(pred==team)  %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(SC_MDT_error = sqrt(sum((standard - prediction)^2) / n())) %>%
    select(cell_line, treatment, time, marker, SC_MDT_error)
  
  SC_MDT_val_score <- SC_MDT_pred %>%
    pull(SC_MDT_error) %>%
    mean()
  
  ## Update score table 
  scores <- scores %>% add_row("CV_loop" = i,
                               "val_CL" = cell_lines[i],
                               "MDT" = MDT_val_score, 
                               "arbiter" = arbiter_val_score,
                               "lc_arbiter" = lc_arbiter_val_score,
                               "lm" = lm_val_score,
                               "single" = SB_val_score,
                               "SC_MDT" = SC_MDT_val_score)
  
  ## Count selected features by MDT and arbiter
  MDT_features <- unique(sub("[<>=].*", "", labels(MDT))) %>% 
    as_tibble() %>% 
    count(value) %>% 
    rename(count_MDT=n, feature = value)
  arbiter_features <- lapply(arbiter, function(x){unique(sub("[<>=].*", "", labels(x)))}) %>% 
    unlist() %>% as_tibble() %>% count(value) %>% rename(feature = value)
  feature_count <- feature_count %>%
    left_join(MDT_features) %>%
    left_join(arbiter_features) %>%
    mutate(count_MDT = ifelse(is.na(count_MDT), 0, count_MDT),
           n=  ifelse(is.na(n), 0, n),
           count = count+count_MDT+n) %>%
    select(feature, count)
  
  # Save all predictions
  predictions <- plyr::join_all(list(MDT_pred, arbiter_pred_value, lc_arbiter_pred, lm_pred, SB_pred, SC_MDT_pred)) %>%
    as_tibble()
  
  all_predictions <- bind_rows(all_predictions, predictions)
}
scores <- filter(scores, !is.na(CV_loop))  
if (FALSE) {
  saveRDS(scores, "./prediction_combinations/SC1/LOO_CV_scores.rds")
  saveRDS(feature_count, "./prediction_combinations/SC1/LOO_CV_feature_counts.rds")
  saveRDS(all_predictions, "./prediction_combinations/SC1/LOO_CV_all_predictions.rds")
  
}

CV_scores <- readRDS("./prediction_combinations/SC1/LOO_CV_scores.rds")
CV_features <- readRDS("./prediction_combinations/SC1/LOO_CV_feature_counts.rds")

colMeans(select(CV_scores, -c(CV_loop, val_CL)))

CV_scores %>%
  gather(model, score, -c(CV_loop, val_CL)) %>%
  ggplot(aes(model, score, fill=model)) +
  geom_boxplot() +
  labs(title = "SC1; leave one out CV") +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

CV_scores %>%
  gather(model, score, -c(CV_loop, val_CL)) %>%
  ggplot(aes(val_CL, score, fill = model)) +
  geom_col(position = "dodge") +
  labs(title = "SC1; leave one out CV") +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.text = element_text(size=13))

CV_features  %>%
  arrange(desc(count)) %>%
  ggplot(aes(reorder(feature, -count), count/(6*23), fill=feature)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90), 
        axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        title = element_text(size=15),
        legend.position = "none") +
  labs(y="proprtion selected", x="feature", title="SC1; leave one out CV")

wilcox.test(CV_scores$lc_arbiter, CV_scores$single, paired = TRUE, alternative = "greater")
t.test(CV_scores$lc_arbiter, CV_scores$single, paired = TRUE)

# To investigate the decision trees
if (FALSE) {
val <- val_data_err %>% select(cell_line, treatment, time, marker, submissions) %>%
  gather(team, true_error, submissions)
arbiter_pred %>% gather(team, pred_error, submissions) %>%
  bind_cols(val) %>%
  filter(team == "icx_bxai") %>%
  ggplot(aes(pred_error, true_error, colour = marker)) + 
  geom_point() +
  facet_wrap(~time)

select(val_data_err, cell_line, treatment, time, marker) %>% bind_cols(arbiter_pred) %>%
  arrange(icx_bxai)

arbiter_pred %>% gather(team, pred_error, submissions)%>%
  bind_cols(val) %>%
  select(-team1) %>%
  group_by(cell_line, treatment, time, marker) %>%
  top_n(n = -1, wt=pred_error) %>%
  arrange(team, cell_line, treatment, time, marker) %>%
  select(cell_line, treatment, time, marker, team, true_error, pred_error) %>%
  pull(pred_error) %>%
  mean()
val_data_err %>%
  add_column(pred = names(arbiter_pred)[max.col(-arbiter_pred)]) %>%
  select(cell_line, treatment, time, marker, pred) %>%
  arrange(pred, cell_line, treatment, time, marker)

plot(MDT)
text(MDT)
var <- val_data_err %>%
  select(-(np_markers))%>% 
  gather(team, true_error, submissions) %>%
  mutate(condID = paste(cell_line, treatment, time, marker, sep = "_")) 
temp <- var %>% filter(team == "icx_bxai")%>% arrange(true_error) %>% pull(condID)
var %>% mutate(condID = factor(condID, levels = temp)) %>%
  ggplot(aes(condID, true_error, colour = team)) +
  geom_point()
    
arbiter_pred %>% gather(team, pred_error, submissions)%>%
  bind_cols(val) %>%
  select(-team1) %>%
  group_by(cell_line, treatment, time, marker) %>%
  top_n(n = -1, wt=pred_error) %>%
  arrange(team, cell_line, treatment, time, marker) %>%
  select(cell_line, treatment, time, marker, team, true_error, pred_error) %>%
  pull(team) %>%
  table()

sub_data_err %>% pull(best_sub) %>% table() %>% as_tibble() %>% arrange(-n)
arbiter_pred %>% gather(team, pred_error, submissions)%>%
  bind_cols(val) %>%
  select(-team1) %>%
  group_by(cell_line, treatment, time, marker) %>%
  top_n(n = -1, wt=pred_error) %>%
  pull(team) %>%
  table() %>%
  sort()

train_data_err %>% pull(best_sub) %>%table() %>% sort()
plot(arbiter[[1]])
text(arbiter[[1]])
}



# Get per subchallenge data on non-predicted conditions
library(tidyverse)
library(RSQLite)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

#To load new cell line
con <- dbConnect(RSQLite::SQLite(), "~/Desktop/BQ internship/BQ/DREAM/single_cell_dream_cls.sqlite")

#------------------------------------------ SC1 ------------------------------------------------------------------------ 
## We want the 32 markers not in the validation data per cell and median per condition
sub_data_value <- readRDS("./submission_data/intermediate_data/sc1_all_predictions.rds")
cell_lines <- unique(sub_data_value$cell_line)

all_cell_line_data <- tibble()
for (i in 1:length(cell_lines)) {
sc_data <- dbReadTable(con,  dbQuoteIdentifier(con, cell_lines[i])) %>%
  as_tibble()
all_cell_line_data <- sub_data_value %>% 
  filter(cell_line == cell_lines[i]) %>%
  select(glob_cellID, cell_line, treatment, time, cellID, fileID) %>%
  distinct() %>%
  left_join(sc_data) %>%
  bind_rows(all_cell_line_data)
}

dbDisconnect(con)
reporters <- c("p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")

# Single cells
all_cell_line_data <- all_cell_line_data %>% select(-reporters)
#saveRDS(all_cell_line_data, "./submission_data/intermediate_data/sc1_all_non_predicted.rds")

# Median per condition
all_median <- all_cell_line_data  %>%
  select(-c(glob_cellID, cellID, fileID)) %>% 
  group_by(cell_line, treatment, time) %>%
  summarise_all(median)
#saveRDS(all_median, "./submission_data/intermediate_data/sc1_median_non_predicted.rds")

submissions <-  readRDS("./submission_data/intermediate_data/sc1_ranked_teams.rds") %>%
  as.character()
sub_data_err <- readRDS("./submission_data/intermediate_data/sc1_rmse_conditions.rds")
sub_data_err <- sub_data_err %>%
  add_column(best_sub = names(sub_data_err[submissions])[max.col(-sub_data_err[submissions])])
sub_data_median <- readRDS("./submission_data/intermediate_data/sc1_median_conditions.rds")
sub_data_median_np <- readRDS("./submission_data/intermediate_data/sc1_median_non_predicted.rds")
sub_data_all <- readRDS("./submission_data/intermediate_data/sc1_all_predictions.rds")
sub_data_all_np <- readRDS("./submission_data/intermediate_data/sc1_all_non_predicted.rds")


np_markers <- sub_data_median_np %>% ungroup() %>% select(-c(cell_line, treatment, time)) %>% colnames()

sub_data_err <- sub_data_err %>% full_join(sub_data_median_np) %>%
  select(cell_line, treatment, time, marker, np_markers, submissions, best_sub)
sub_data_median <- sub_data_median %>% full_join(sub_data_median_np) %>%
  select(cell_line, treatment, time, marker, np_markers, standard, submissions)
sub_data_all <- sub_data_all %>% full_join(sub_data_all_np)

saveRDS(sub_data_err, "./submission_data/intermediate_data/sc1_condErr_medianVals.rds")
saveRDS(sub_data_median, "./submission_data/intermediate_data/sc1_median_conditions_np.rds")
saveRDS(sub_data_all, "./submission_data/intermediate_data/sc1_all_NP_predictions.rds")




#------------------------------------------ SC2 ------------------------------------------------------------------------ 
required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 
submissions <-  readRDS("./submission_data/intermediate_data/sc2_ranked_teams.rds") %>%
  as.character()
sub_data_err <- readRDS("./submission_data/intermediate_data/sc2_stats_sumSquared_conditions.rds")
sub_data_err <- sub_data_err %>%
  add_column(best_sub = names(sub_data_err[submissions])[max.col(-sub_data_err[submissions])]) %>%
  separate(cond_id, into = c("cell_line", "treatment", "time"), sep = "_") %>%
  mutate(time=ifelse(time==16, 17, time))
sub_data_values <- readRDS("./submission_data/intermediate_data/sc2_stats_conditions.rds") %>%
  mutate(stat_variable = sub("v_","v-", stat_variable),
         stat_variable = sub("n_","n-", stat_variable)) %>%
  mutate(time=ifelse(time==16, 17, time))

cell_lines <- unique(sub_data_values$cell_line)

all_cell_line_data <- tibble()
for (i in 1:length(cell_lines)) {
  sc_data <- dbReadTable(con,  dbQuoteIdentifier(con, cell_lines[i])) %>% 
    as_tibble() %>%
    filter(treatment == "EGF" & time == 0) %>% 
    select(required_columns)  %>%
    group_by(cell_line) %>%
    summarise_at(vars(-time, -treatment, -group_cols()), median) %>%
    ungroup()
  all_cell_line_data <- bind_rows(all_cell_line_data, sc_data)
}

colnames(all_cell_line_data) <- c("cell_line", paste0("EGF0_", colnames(all_cell_line_data)[2:36]))

err_median_EGF <- sub_data_err %>% 
  mutate(time=as.numeric(time)) %>%
  left_join(all_cell_line_data)

values_median_EGF0 <- sub_data_values %>% 
  mutate(time=as.numeric(time)) %>%
  left_join(all_cell_line_data)

saveRDS(err_median_EGF, "./submission_data/intermediate_data/sc2_err_median_EGF0.rds")
saveRDS(values_median_EGF0, "./submission_data/intermediate_data/sc2_values_median_EGF0.rds")
  
## ------------------------------------- SC3 --------------------------------------

required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 
submissions <-  readRDS("./submission_data/intermediate_data/sc3_ranked_teams.rds") %>%
  as.character()
sub_data_err <- readRDS("./submission_data/intermediate_data/sc3_stats_sumSquared_conditions.rds")
sub_data_err <- sub_data_err %>%
  add_column(best_sub = names(sub_data_err[submissions])[max.col(-sub_data_err[submissions])]) %>%
  separate(cond_id, into = c("cell_line", "treatment", "time"), sep = "_") %>%
  mutate(time=ifelse(time==16, 17, time))
sub_data_values <- readRDS("./submission_data/intermediate_data/sc3_stats_conditions.rds") %>%
  mutate(stat_variable = sub("v_","v-", stat_variable),
         stat_variable = sub("n_","n-", stat_variable)) %>%
  mutate(time=ifelse(time==16, 17, time))

cell_lines <- unique(sub_data_values$cell_line)

all_cell_line_data <- tibble()
for (i in 1:length(cell_lines)) {
  sc_data <- dbReadTable(con,  dbQuoteIdentifier(con, cell_lines[i])) %>% 
    as_tibble() %>%
    filter(treatment == "EGF" & time == 0) %>% 
    select(required_columns)  %>%
    group_by(cell_line) %>%
    summarise_at(vars(-time, -treatment, -group_cols()), median) %>%
    ungroup()
  all_cell_line_data <- bind_rows(all_cell_line_data, sc_data)
}

colnames(all_cell_line_data) <- c("cell_line", paste0("EGF0_", colnames(all_cell_line_data)[2:36]))

err_median_EGF <- sub_data_err %>% 
  mutate(time=as.numeric(time)) %>%
  left_join(all_cell_line_data) %>%
  mutate(time =ifelse(time==14, 13, time))

values_median_EGF0 <- sub_data_values %>% 
  mutate(time=as.numeric(time)) %>%
  left_join(all_cell_line_data) %>%
  mutate(time =ifelse(time==14, 13, time))

saveRDS(err_median_EGF, "./submission_data/intermediate_data/sc3_err_median_EGF0.rds")
saveRDS(values_median_EGF0, "./submission_data/intermediate_data/sc3_values_median_EGF0.rds")

# -------------------------------------------- SC4 ------------------------------------------------
## Median values of all markers at full t=0 per cell line
submissions <-  readRDS("./submission_data/intermediate_data/sc4_ranked_teams.rds") %>%
  as.character()
sub_data_values <- readRDS("./submission_data/intermediate_data/sc4_combined_data.rds")
cell_lines <- unique(sub_data_values$cell_line)

all_cell_line_data <- tibble()
for (i in 1:length(cell_lines)) {
  sc_data <- dbReadTable(con,  dbQuoteIdentifier(con, cell_lines[i])) %>% 
    as_tibble() %>%
    filter(treatment == "full" & time == 0) %>% 
    select(-c(fileID, cellID, time, treatment))  %>%
    group_by(cell_line) %>%
    summarise_at(vars(-group_cols()), median) %>%
    ungroup()
  all_cell_line_data <- bind_rows(all_cell_line_data, sc_data)
}
dbDisconnect(con)
colnames(all_cell_line_data) <- c("cell_line", paste0("full0_", colnames(all_cell_line_data)[2:38]))

# Caldulate error per team and condition (CL, treatment, time, marker)
# Then find best performing team per condition
sub_data_err <- sub_data_values %>% mutate_at(submissions, ~sqrt((standard-.)^2))
sub_data_err <- sub_data_err %>%
  add_column(best_sub = names(sub_data_err[submissions])[max.col(-sub_data_err[submissions])])

# Save both the predicted values and error data frames including the median marker values at full0
median_and_values <- sub_data_values %>% left_join(all_cell_line_data)
median_and_error <- sub_data_err %>% left_join(all_cell_line_data)
saveRDS(median_and_values, "./submission_data/intermediate_data/sc4_median_and_values.rds")
saveRDS(median_and_error, "./submission_data/intermediate_data/sc4_median_and_error.rds")

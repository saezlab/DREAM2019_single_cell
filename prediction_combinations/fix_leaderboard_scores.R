## Fix leaderboard

library(tidyverse)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell")

# ------------------------ SC1 -----------------------------
submission_folder <- "./submission_data/final/SC1"
leader_board <- read_csv(file.path(submission_folder, "leaderboard_final_sc1.csv"))
new_LB <-leader_board %>% mutate(score = case_when(score > 1 ~ score*1e-6,
                                          TRUE ~ score))
new_LB$score  
leader_board$score
write_csv(new_LB, file.path(submission_folder, "leaderboard_final_sc1.csv"))

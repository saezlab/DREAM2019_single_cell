


library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(ggrepel)

source("./submission_analysis/utilities.R")

colors  = brewer.pal(n = 8,"Dark2")

LB_sc1 <- read_csv("./submission_data/final_round/leaderboard_all_rounds_sc1.csv") %>%
	fix_leaderboard_raw() %>% filter(!is.na(round)) %>%
	mutate(createdOn = fix_submission_dates(createdOn)) %>%
	mutate(round = paste0("round_",round)) %>% mutate(sc=1)
LB_sc2 <- read_csv("./submission_data/final_round/leaderboard_all_rounds_sc2.csv") %>%
	fix_leaderboard_raw() %>%  filter(!is.na(round)) %>%
	mutate(createdOn = fix_submission_dates(createdOn)) %>%
	mutate(round = paste0("round_",round)) %>% mutate(sc=2)
LB_sc3 <- read_csv("./submission_data/final_round/leaderboard_all_rounds_sc3.csv") %>%
	fix_leaderboard_raw() %>%  filter(!is.na(round)) %>%
	mutate(createdOn = fix_submission_dates(createdOn)) %>%
	mutate(round = paste0("round_",round)) %>% mutate(sc=3)
LB_sc4 <- read_csv("./submission_data/final_round/leaderboard_all_rounds_sc4.csv") %>%
	fix_leaderboard_raw() %>%  filter(!is.na(round)) %>%
	mutate(createdOn = fix_submission_dates(createdOn)) %>%
	mutate(round = paste0("round_",round)) %>% mutate(sc=4)



LB_all <- LB_sc1 %>% bind_rows(LB_sc2, LB_sc3,LB_sc4)

## RANK them based on their best prediction, NOT official

LB_all %>% filter(status == "SCORED") %>% 
	group_by(submitterId, sc) %>%
	top_n(-1, score) %>%
	filter(!duplicated(submitterId,sc)) %>%  # for duplicates
	ungroup() %>%
	arrange(sc,score) %>%
	select(-objectId, -status, -FAILURE_REASON) %>%
	group_by(sc) %>% mutate(sc_rank = rank(score)) %>%
	select(submitterId, sc,sc_rank) %>% 
	group_by(submitterId) %>% 
	summarise(subchallenges = paste(sc,collapse = ", "),
			  ranks = paste(sc_rank,collapse = ", "), 
			  min_rank = min(sc_rank)) %>%
	arrange(min_rank) %>%
	select(-min_rank) %>%
	rowwise() %>%
	do(cat(paste(.,"\t"),"\n"))


### 
LB_all %>% filter(sc == 1 ) %>% select(submitterId, round) %>% table()
LB_all %>% filter(sc == 2 ) %>% select(submitterId, round) %>% table()
LB_all %>% filter(sc == 3 ) %>% select(submitterId, round) %>% table()
LB_all %>% filter(sc == 4 ) %>% select(submitterId, round) %>% table()


LB_all %>% select(sc, round, submitterId) %>% table()
	


## Find the best submission for teams that didnt submit to round 3 (for forwarding):

LB_all %>% filter(submitterId =="Raghava_India_SCS") %>% group_by(sc) %>% top_n(-1,score)

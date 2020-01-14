# evaluating the rankings of the top performing teams in SC2
#
# author: Attila Gabor
# 
# Short summary:
# In this subchallenge the score is based on squared sum of the modelling 
# errors (distance in median and covariance) across all conditions.
# To assess if the small differences between the top performing teams are 
# significant we utilise bootstrap sampling of the conditions. We resample the
# experimental conditions (defined by the triplet: cell-line, treatment and time)
# with replacement and calculate the score for each team. 
# 
# Finally we compute the Bayes factor over the bootstraps for
# consecutively ranked teams
# B:= [sum(score(teamA) < score(teamB))] / [sum(score(teamA) > score(teamB))]
# we consider a Bayes factor larger than 3 significant 


library(tidyverse)
library(purrr)
source("./scoring_scripts/score_sc2.R")
# SC2 --------------------------------------------------------------------------

# make sure to download the data of each team to this folder: 
# naming convention for files: `submissionID`.csv
submission_folder = "./submission_data/final_round/SC2/"

# read leaderboard
SC_leaderboard = read_csv(file.path(submission_folder,"leaderboard_final_sc2.csv")) %>%
	select(-writeUp, -createdOn) %>% 
	mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
	mutate(submitterId = make.names(submitterId))

# we take a subset of teams from the top of the rankings. 
N_top <- 16
ranked_teams <- factor(SC_leaderboard$submitterId[1:N_top],levels = SC_leaderboard$submitterId[1:N_top])

SC_leaderboard %>% 
	ggplot(aes(factor(submitterId,levels = SC_leaderboard$submitterId),score)) +
	geom_point() + 
	theme_bw() +
	xlab("Teams") + 
	ylab("RMSE error") + 
	ggtitle("Final leaderboard for subchallenge II") + 
	theme(axis.text.x = element_text(hjust=1,angle = 45))
ggsave("./submission_analysis/figures/SC2_leaderboard.pdf",width = 6,height = 5)


# read team's predictions from csv files and compute the stats
required_columns <- c('cell_line','treatment', 'time',
					  'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
					  'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
					  'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
					  'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
					  'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
					  'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
					  'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
					  'p.STAT3', 'p.STAT5') 
#' read_to_stats
#' 
#' reads data and compute stats --> we only work with the stats afterwards
read_to_stats <- function(file_name){
	read_csv(file_name) %>% select(required_columns) %>% data_to_stats(.)
}

prediction_data <- SC_leaderboard  %>%
 	slice(1:N_top) %>% 
	mutate(predictions = map(file.path(submission_folder,submissions),read_to_stats))


#' order_predictions
#'
#' takes a prediction matrix, reformats to long format and keeps only the values.
#' 
#' @param stat_matrix  matrix with the statistics in tibble format
#' @return value vector gathered over reporters and sorted by glob_cellID
order_stats <- function(stat_matrix){
	
	stat_matrix %>% arrange(cell_line, treatment,  time, stat_variable) %>%
		select(stat_value)
}

# order the predictions the same way for all teams and for the golden standard
# so a simple column bind will merge well.

ordered_predictions <- prediction_data %>% mutate(predictions = map(predictions,order_stats)) %>%
	do(bind_cols(.$predictions))
names(ordered_predictions) <- prediction_data$submitterId

# read golden standard
gs <- read_to_stats("./challenge_data/validation_data/sc2gold.csv") %>%
	arrange(cell_line, treatment,  time, stat_variable) %>% 
	rename(standard=stat_value)

## combined data: 
# columns describe the conditions, 
# 	standard:  stats computed from the golden standard measurement

combined_statistics <- bind_cols(gs,ordered_predictions)

# calculate the RMSE error for each conditions for each team. 
RMSE_conditions <- combined_statistics %>% 
	summarise_at(as.character(ranked_teams),~ sqrt(sum((standard - .)^2)))


## Bootstrapping: 
# we repeat N=1000 times the boostrapping from the conditions and compare the 
# the averaged RMSE over the samples.
N_bootstrap <- 1000
set.seed(123)

# add ID for conditions
combined_statistics <- combined_statistics %>%
	mutate(cond_id = paste(cell_line,treatment,time,sep = "_")) %>%
	select(cond_id,everything())


# for each bootstrap:
# 1. we calculate ahead the squared diff between the stats for each condition
# 2. for each BS sample select conditions with replacement
# 3. sum up the squared differences and take the square root of the sum 

condition_stats_SumSquared <- combined_statistics %>% 
	group_by(cond_id) %>% 
	summarise_at(as.character(ranked_teams),~ sum((standard - .)^2))


bootstrap_stats <- tibble(BS_sample = seq(N_bootstrap)) %>%
	mutate( map(BS_sample, .f = ~ condition_stats_SumSquared %>% 
					sample_frac(size = 1, replace = TRUE) %>%
					summarise_at(as.character(ranked_teams),~ sqrt(sum(.))))
	) %>% unnest()



# Plot the performance of the team on the bootstrap samples
bootstrap_stats %>% 
	gather(teams,score,-BS_sample) %>%
	mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
	ggplot() + geom_boxplot(aes(teams, score,color=teams)) + 
	ggtitle("Subchallenge II: performance over bootstrap samples") + 
	geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) + 
	theme_bw() +
	theme(axis.text.x = element_text(hjust=1,angle = 45))
ggsave("./submission_analysis/figures/SCII_score_bootstrap.pdf",width = 10,height = 6)

# show the ranks across Bootstraps
bootstrap_stats %>% ungroup() %>%
	gather(team, score, as.character(ranked_teams), -BS_sample) %>%
	group_by(BS_sample) %>% mutate(BS_ranks = rank(score)) %>%
	ggplot(aes(factor(team,levels = levels(ranked_teams)), BS_ranks)) +
	geom_violin(aes(fill=team),scale = "width",draw_quantiles = 0.5) + 
	ggtitle("Subchallenge II: ranks over bootstrap samples") + 
	geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,rank(score))) + 
	theme_bw() +
	xlab("teams") + 
	ylab("ranks") + 
	theme(axis.text.x = element_text(hjust=1,angle = 45))
ggsave("./submission_analysis/figures/SCII_score_bootstrap_ranks.pdf",width = 10,height = 6)


# Compute the Bayes score based on the score: how many times teamA is better than 
# teamB over the number of time the opposit happened: 
compute_bayes_rank <- function(team_a,team_b,bootstrap_table){
	team_a <- as.character(team_a)
	team_b <- as.character(team_b)
	sum(bootstrap_table[,team_a]<bootstrap_table[,team_b])/sum(bootstrap_table[,team_a]>bootstrap_table[,team_b])
}


Bayes_factors <- tibble(team_A = ranked_teams[-16], team_B = ranked_teams[-1]) %>%
	rowwise() %>%
	mutate(Bayes_factor = map2_dbl(team_A,team_B,compute_bayes_rank,bootstrap_stats))


Bayes_factors %>% mutate(Bayes_factor = round(Bayes_factor,1),
						 Bayes_factor = ifelse(
						 	is.infinite(Bayes_factor),
						 	">500",as.character(Bayes_factor))) %>%
	write_tsv("./submission_analysis/figures/SCII_bayes_factors.tsv",col_names = T)

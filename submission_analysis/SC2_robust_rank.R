# evaluating the rankings of the top performing teams in SC1
#
# author: Attila Gabor
# 
# Short summary:
# In each of the subchallenges the scores are based on averaging the modelling 
# errors across all conditions.
# To assess if the small differences between the top performing teams are 
# significant we utilise bootstrap sampling of the conditions. Then we calculate
# the modelling error on each sample. Finally we compute the Bayes factor


library(tidyverse)
library(purrr)
source("./scoring_scripts/score_sc2.R")
# SC2 --------------------------------------------------------------------------

# we take a subset of teams from the top of the rankings. 
N_top <- 10

# make sure to download the data of each team to this folder: 
# naming convention for files: `submissionID`.csv
submission_folder = "./submission_data/final_round/SC2/"

# read leaderboard
SC_leaderboard = read_csv(file.path(submission_folder,"leaderboard_final_sc2.csv")) %>%
	select(-writeUp, -createdOn) %>% 
	mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
	mutate(submitterId = make.names(submitterId))

ranked_teams <- factor(SC_leaderboard$submitterId[1:N_top],levels = SC_leaderboard$submitterId[1:N_top])

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
	#group_by(cell_line, treatment, time, marker) %>%
	summarise_at(as.character(ranked_teams),~ sqrt(sum((standard - .)^2)))


## Bootstrapping: 
# we repeat N=1000 times the boostrapping from the conditions and compare the 
# the averaged RMSE over the samples.

N_bootstrap <- 1000
set.seed(123)

bootstrap_stats <- tibble(BS_sample = seq(N_bootstrap)) %>%
	mutate( map(BS_sample, .f = ~ combined_statistics %>% 
						 	sample_frac(size = 1, replace = TRUE) %>%
						 	summarise_at(as.character(ranked_teams),~ sqrt(sum((standard - .)^2)))
				)) %>% unnest()

# order the boostraps by the best team and plot RMSE vs Bootstraps
bootstrap_stats %>% arrange(icx_bxai) %>% 
	mutate(BS_sample = 1:N_bootstrap) %>%
	gather(teams,RMSE,-BS_sample) %>%
	mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
	ggplot() + geom_point(aes(BS_sample, RMSE,color=teams)) + 
	theme_bw() + 
	ggtitle("Subchallenge II: performance for each bootstrap sample")

# plotting by boxplot
bootstrap_stats %>% 
	gather(teams,RMSE,-BS_sample) %>%
	mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
	ggplot() + geom_boxplot(aes(teams, RMSE,color=teams)) + 
	ggtitle("Subchallenge II: performance over bootstrap samples") + 
	geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) + 
	theme_bw() +
	theme(axis.text.x = element_text(hjust=1,angle = 45))
ggsave("./submission_analysis/figures/SCII_score_bootstrap.pdf",width = 10,height = 6)

Bayes_factors <- bootstrap_stats %>% 
	summarise(icx_bxai_vs_pqiu = sum(icx_bxai<X.pqiu)/sum(icx_bxai>X.pqiu), 
			  pqiu_vs_orangeballs = sum(X.pqiu<orangeballs)/sum(X.pqiu>orangeballs), 
			  orangeballs_vs_PaL = sum(orangeballs<PaL)/sum(orangeballs>PaL), 
			  PaL_vs_NAD = sum(PaL<NAD)/sum(PaL>NAD) )

# Bayes_factors <- bootstrap_stats %>% 
# 	summarise(icx_bxai_vs_pqiu = sum(icx_bxai<X.pqiu)/N_bootstrap, 
# 			  pqiu_vs_orangeballs = sum(X.pqiu<orangeballs)/N_bootstrap, 
# 			  orangeballs_vs_PaL = sum(orangeballs<PaL)/N_bootstrap, 
# 			  PaL_vs_NAD = sum(PaL<NAD)/N_bootstrap )

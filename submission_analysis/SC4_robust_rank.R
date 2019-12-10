# evaluating the rankings of the top performing teams in SC4
#
# author: Attila Gabor
# 
# Short summary:
# In the subchallenge the score is based on the RMSE averaged  
# across all conditions.
# To assess if the small differences between the top performing teams are 
# significant we utilise bootstrap sampling of the conditions. We resample the
# experimental conditions (defined by cell-line and treatment)
# with replacement and calculate the score for each team. 
# 
# Finally we compute the Bayes factor over the bootstraps for
# consecutively ranked teams
# B:= [sum(score(teamA) < score(teamB))] / [sum(score(teamA) > score(teamB))]
# we consider a Bayes factor larger than 3 significant 



library(tidyverse)
library(purrr)

# SC4 --------------------------------------------------------------------------


# make sure to download the data of each team to this folder: 
# naming convention for files: `submissionID`.csv
submission_folder = "./submission_data/final_round/SC4/"

# read leaderboard
SC_leaderboard = read_csv(file.path(submission_folder,"leaderboard_final_sc4.csv")) %>%
	select(-writeUp, -createdOn) %>% 
	mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
	mutate(submitterId = make.names(submitterId))

# we take a subset of teams from the top of the rankings. 
N_top <- 20

ranked_teams <- factor(SC_leaderboard$submitterId[1:N_top],levels = SC_leaderboard$submitterId[1:N_top])

SC_leaderboard %>%  
	ggplot(aes(factor(submitterId,levels = SC_leaderboard$submitterId),score)) +
	geom_point() + 
	theme_bw() +
	xlab("Teams") + 
	ylab("RMSE error") + 
	ggtitle("Final leaderboard for subchallenge IV") + 
	theme(axis.text.x = element_text(hjust=1,angle = 45))
ggsave("./submission_analysis/figures/SC4_leaderboard.pdf",width = 10,height = 6)


# read team's predictions from csv files
prediction_data <- SC_leaderboard  %>%
 	slice(1:N_top) %>% 
	mutate(predictions = map(file.path(submission_folder,submissions),read_csv))
	
reporters <- c('b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
					  'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
					  'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
					  'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
					  'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
					  'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
					  'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
					  'p.STAT3', 'p.STAT5') 

#' order_predictions
#'
#' takes a prediction matrix, reformats to long format and keeps only the values.
#' 
#' @param predictions prediction matrix in tibble format
#' @return value vector gathered over reporters and sorted by glob_cellID
order_predictions <- function(predictions){
	
	predictions %>%	arrange(cell_line, treatment,  time) %>%
		gather("reporter","value",reporters) %>% select(value)
}

# order the predictions the same way for all teams and for the golden standard
# so a simple column bind will merge well.

ordered_predictions <- prediction_data %>% mutate(predictions = map(predictions,order_predictions)) %>%
	do(bind_cols(.$predictions))
names(ordered_predictions) <- prediction_data$submitterId

# read golden standard
gs <- read_csv("./challenge_data/validation_data/sc4gold.csv") %>%
	arrange(cell_line, treatment,  time) %>%
	gather(key = "marker", value = "standard", reporters)

## combined data: 
# columns describe the conditions, 
# 	standard: golden standard measurement

combined_data <- bind_cols(gs,ordered_predictions)

# calculate the RMSE error for each conditions for each team. 
RMSE_conditions <- combined_data %>% 
	group_by(cell_line, treatment, marker) %>%
	summarise_at(as.character(ranked_teams),~ sqrt(sum((standard - .)^2) / n())) %>%
	ungroup()
	
# check against the leaderboard, should be the same: 
RMSE_conditions %>% summarise_at(as.character(ranked_teams),mean) 

## Bootstrapping: 
# we repeat N=1000 times the boostrapping from the conditions and compare the 
# the averaged RMSE over the samples.

N_bootstrap <- 1000
set.seed(123)

bootstrap_RMSE <- tibble(BS_sample = seq(N_bootstrap)) %>%
	mutate( map(BS_sample, .f = ~ RMSE_conditions %>% 
						 	sample_frac(size = 1, replace = TRUE) %>%
						 	summarise_at(as.character(ranked_teams),mean))) %>% unnest()


# plotting by boxplot
bootstrap_RMSE %>% 
	gather(teams,RMSE,-BS_sample) %>%
	mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
	ggplot() + geom_boxplot(aes(teams, RMSE,color=teams)) + 
	ggtitle("Subchallenge IV: performance over bootstrap samples") + 
	geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) +
	theme_bw() +
	theme(axis.text.x = element_text(hjust=1,angle = 45))

ggsave("./submission_analysis/figures/SC4_score_bootstrap.pdf",width = 10,height = 6)



# show the ranks across Bootstraps
bootstrap_RMSE %>% ungroup() %>%
	gather(team, score, as.character(ranked_teams), -BS_sample) %>%
	group_by(BS_sample) %>% mutate(BS_ranks = rank(score)) %>%
	ggplot(aes(factor(team,levels = levels(ranked_teams)), BS_ranks)) +
	# geom_boxplot(aes(color=team)) +
	geom_violin(aes(fill=team),scale = "width",draw_quantiles = 0.5) + 
	# geom_jitter(aes(color=team)) + 
	ggtitle("Subchallenge IV: ranks over bootstrap samples") + 
	geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,rank(score))) + 
	theme_bw() +
	xlab("teams") + 
	ylab("ranks") + 
	theme(axis.text.x = element_text(hjust=1,angle = 45))

ggsave("./submission_analysis/figures/SCIV_score_bootstrap_ranks.pdf",width = 10,height = 6)


# Compute the Bayes score based on the score: how many times teamA is better than 
# teamB over the number of time the opposit happened: 
compute_bayes_rank <- function(team_a,team_b,bootstrap_table){
	team_a <- as.character(team_a)
	team_b <- as.character(team_b)
	sum(bootstrap_table[,team_a]<bootstrap_table[,team_b])/sum(bootstrap_table[,team_a]>bootstrap_table[,team_b])
}


Bayes_factors <- tibble(team_A = ranked_teams[-20], team_B = ranked_teams[-1]) %>%
	rowwise() %>%
	mutate(Bayes_factor = map2_dbl(team_A,team_B,compute_bayes_rank,bootstrap_RMSE))


Bayes_factors %>% mutate(Bayes_factor = round(Bayes_factor,1),
						 Bayes_factor = ifelse(
						 	is.infinite(Bayes_factor),
						 	">500",as.character(Bayes_factor))) %>%
	write_tsv("./submission_analysis/figures/SCIV_bayes_factors.tsv",col_names = T)


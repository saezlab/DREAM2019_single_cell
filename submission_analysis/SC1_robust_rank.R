# evaluating the rankings of the top performing teams in SC1
#
# author: Attila Gabor
# 
# Short summary: Bootstrap ranking
# In the subchallenge the score is based on averaging the modelling 
# errors across all conditions.
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

# SC1 --------------------------------------------------------------------------


# make sure to download the data of each team to this folder: 
# naming convention for files: `submissionID`.csv
submission_folder = "./submission_data/final_round/SC1/"

# read leaderboard
SC_leaderboard = read_csv(
		file.path(submission_folder,"leaderboard_final_sc1.csv")) %>%
	select(-writeUp, -createdOn) %>% 
	mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
	mutate(submitterId = make.names(submitterId))

# Team KAUST_RSS didnt follow the template but managed to pass the validation 
SC_leaderboard <- SC_leaderboard %>% filter(submitterId != "KAUST_RSS") 

# we take a subset of teams from the top of the rankings. 
N_top <- nrow(SC_leaderboard)

ranked_teams <- factor(SC_leaderboard$submitterId[1:N_top],
					   levels = SC_leaderboard$submitterId[1:N_top])

# submission statistics
SC_leaderboard %>% summarise(mean = mean(score),
                             std = sd(score),
                             median = median(score),
                             range_L = range(score)[[1]],
                             range_U = range(score)[[2]])
                            


SC_leaderboard %>%  
	ggplot(aes(factor(submitterId,levels = SC_leaderboard$submitterId),score)) +
	geom_point() + 
	theme_bw() +
	xlab("Teams") + 
	ylab("RMSE error") + 
	ggtitle("Final leaderboard for subchallenge I") + 
	theme(axis.text.x = element_text(hjust=1,angle = 45))
ggsave("./submission_analysis/figures/SC_leaderboard.pdf",width = 6,height = 5)

# read team's predictions from csv files
prediction_data <- SC_leaderboard  %>%
 	slice(1:N_top) %>% 
	mutate(predictions = map(file.path(submission_folder,submissions),read_csv))
	
reporters <- c("p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")

#' order_predictions
#'
#' takes a prediction matrix, reformats to long format and keeps only the values.
#' 
#' @param predictions prediction matrix in tibble format
#' @return value vector gathered over reporters and sorted by glob_cellID
order_predictions <- function(predictions){
	
	predictions %>% select(glob_cellID,reporters) %>% 
		arrange(glob_cellID) %>%
		gather("reporter","value",reporters) %>% select(value)
}

# order the predictions the same way for all teams and for the golden standard
# so a simple column bind will merge well.


ordered_predictions <- prediction_data %>%
	mutate(predictions = map(predictions,order_predictions)) %>%
	do(bind_cols(.$predictions))
names(ordered_predictions) <- prediction_data$submitterId

# read golden standard
gs <- read_csv("./challenge_data/validation_data/sc1gold.csv") %>%
	arrange(glob_cellID) %>%
	gather(key = "marker", value = "standard", reporters)

## combined data: 
# columns describe the conditions, 
# 	standard: golden standard measurement

combined_data <- bind_cols(gs,ordered_predictions)

# calculate the RMSE error for each conditions for each team. 
RMSE_conditions <- combined_data %>% 
	group_by(cell_line, treatment, time, marker) %>%
	summarise_at(as.character(ranked_teams),~ sqrt(sum((standard - .)^2) / n())) %>%
	ungroup()

median_conditions <- combined_data %>% 
	group_by(cell_line, treatment, time, marker) %>%
	summarise_at(c("standard",as.character(ranked_teams)),~ median(.)) %>%
	ungroup()

## Bootstrapping: 
# we repeat N=1000 times the boostrapping from the conditions and compare the 
# the averaged RMSE over the samples.

N_bootstrap <- 1000
set.seed(123)

bootstrap_RMSE <- tibble(BS_sample = seq(N_bootstrap)) %>%
	mutate( map(BS_sample, .f = ~ RMSE_conditions %>% 
						 	sample_frac(size = 1, replace = TRUE) %>%
						 	summarise_at(as.character(ranked_teams),mean))) %>% unnest()



# save intermediate summary results for postchallange analyis
if(FALSE){
	write_rds(RMSE_conditions,"./submission_analysis/intermediate_data/sc1_rmse_conditions.rds")
	write_rds(median_conditions,"./submission_analysis/intermediate_data/sc1_median_conditions.rds")
	write_rds(bootstrap_RMSE,"./submission_analysis/intermediate_data/sc1_bootstrap_rmse.rds")
	write_rds(ranked_teams, "./submission_analysis/intermediate_data/sc1_ranked_teams.rds")
}

bootstrap_RMSE <- read_rds("./submission_analysis/intermediate_data/sc1_bootstrap_rmse.rds")
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc1_ranked_teams.rds")


bootstrap_RMSE %>% 
    gather(teams,RMSE,-BS_sample) %>% group_by(teams) %>%
    summarise(mean_RMSE_BS = mean(RMSE),
              std_RMSE_BS = sd(RMSE))

# Plot the performance of the team on the bootstrap samples
bootstrap_RMSE %>% 
	gather(teams,RMSE,-BS_sample) %>%
	mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
	ggplot() + geom_boxplot(aes(teams, RMSE,color=teams)) + 
	ggtitle("Subchallenge I: performance over bootstrap samples") + 
	geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) +
	theme_bw() +
	theme(axis.text.x = element_text(hjust=1,angle = 45))

ggsave("./submission_analysis/figures/SCI_score_bootstrap.pdf",width = 10,height = 6)


# Rank the teams based on the bootstraps and plot the rank distribution 
bootstrap_RMSE %>% ungroup() %>%
	gather(team, score, as.character(ranked_teams), -BS_sample) %>%
	group_by(BS_sample) %>% mutate(BS_ranks = rank(score)) %>%
	ggplot(aes(factor(team,levels = levels(ranked_teams)), BS_ranks)) +
	#geom_boxplot(aes(color=team)) +
	geom_violin(aes(fill=team),scale = "width",draw_quantiles = 0.5) + 
	ggtitle("Subchallenge I: ranks over bootstrap samples") + 
	geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,rank(score))) + 
	theme_bw() +
	xlab("teams") + 
	ylab("ranks") + 
	theme(axis.text.x = element_text(hjust=1,angle = 45))

ggsave("./submission_analysis/figures/SCI_score_bootstrap_ranking.pdf",width = 10,height = 6)


# Compute the Bayes score based on the score: how many times teamA is better than 
# teamB over the number of time the opposit happened: 
compute_bayes_rank <- function(team_a,team_b,bootstrap_table){
	team_a <- as.character(team_a)
	team_b <- as.character(team_b)
	sum(bootstrap_table[,team_a]<bootstrap_table[,team_b])/sum(bootstrap_table[,team_a]>bootstrap_table[,team_b])
}


Bayes_factors <- tibble(team_A = ranked_teams[-22], team_B = ranked_teams[-1]) %>%
	rowwise() %>%
	mutate(Bayes_factor = map2_dbl(team_A,team_B,compute_bayes_rank,bootstrap_RMSE))


Bayes_factors %>% mutate(Bayes_factor = round(Bayes_factor,1),
						 Bayes_factor = ifelse(
						 	is.infinite(Bayes_factor),
						 	">500",as.character(Bayes_factor))) %>%
	write_tsv("./submission_analysis/figures/SCI_bayes_factors.tsv",col_names = T)

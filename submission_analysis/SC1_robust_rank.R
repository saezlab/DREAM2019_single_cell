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

# SC1 --------------------------------------------------------------------------

# we take a subset of teams from the top of the rankings. 
N_top <- 10

# make sure to download the data of each team to this folder: 
# naming convention for files: `submissionID`.csv
submission_folder = "./submission_data/final_round/SC1/"

# read leaderboard
SC1_leaderboard = read_csv(file.path(submission_folder,"leaderboard_final_sc1.csv")) %>%
	select(-writeUp, -createdOn) %>% 
	mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
	mutate(submitterId = make.names(submitterId))

ranked_teams <- factor(SC1_leaderboard$submitterId[1:N_top],levels = SC1_leaderboard$submitterId[1:N_top])

# read team's predictions from csv files
prediction_data <- SC1_leaderboard  %>%
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

ordered_predictions <- prediction_data %>% mutate(predictions = map(predictions,order_predictions)) %>%
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


## Bootstrapping: 
# we repeat N=1000 times the boostrapping from the conditions and compare the 
# the averaged RMSE over the samples.

N_bootstrap <- 1000
set.seed(123)

bootstrap_RMSE <- tibble(BS_sample = seq(N_bootstrap)) %>%
	mutate( map(BS_sample, .f = ~ RMSE_conditions %>% 
						 	sample_frac(size = 1, replace = TRUE) %>%
						 	summarise_at(as.character(ranked_teams),mean))) %>% unnest()

# order the boostraps by the best team and plot RMSE vs Bootstraps
bootstrap_RMSE %>% arrange(icx_bxai) %>% 
	mutate(BS_sample = 1:N_bootstrap) %>%
	gather(teams,RMSE,-BS_sample) %>%
	mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
	ggplot() + geom_point(aes(BS_sample, RMSE,color=teams)) + 
	theme_bw() + 
	ggtitle("Subchallenge I: performance for each bootstrap sample")

# plotting by boxplot
bootstrap_RMSE %>% 
	gather(teams,RMSE,-BS_sample) %>%
	mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
	ggplot() + geom_boxplot(aes(teams, RMSE,color=teams)) + 
	ggtitle("Subchallenge I: performance over bootstrap samples") + 
	geom_point(data = SC1_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) +
	theme_bw() +
	theme(axis.text.x = element_text(hjust=1,angle = 45))

ggsave("./submission_analysis/figures/SCI_score_bootstrap.pdf",width = 10,height = 6)

Bayes_factors <- bootstrap_RMSE %>% 
	summarise(icx_bxai_vs_NAD = sum(icx_bxai<NAD)/sum(icx_bxai>NAD), 
			  NAD_vs_SingleCellLand = sum(NAD<SingleCellLand)/sum(NAD>SingleCellLand), 
			  SingleCellLand_vs_PaL = sum(SingleCellLand<PaL)/sum(SingleCellLand>PaL), 
			  PaL_vs_Raghava_India_SCS = sum(PaL<Raghava_India_SCS)/sum(PaL>Raghava_India_SCS) )

Bayes_factors <- bootstrap_RMSE %>% 
	summarise(icx_bxai_vs_NAD = sum(icx_bxai<NAD)/1000, 
			  NAD_vs_SingleCellLand = sum(NAD<SingleCellLand)/1000, 
			  SingleCellLand_vs_PaL = sum(SingleCellLand<PaL)/1000, 
			  PaL_vs_Raghava_India_SCS = sum(PaL<Raghava_India_SCS)/1000 )

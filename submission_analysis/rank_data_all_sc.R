# Plot the Bootstrap results together with the data usage heatmap 

library(tidyverse)
library(cowplot)

### SC1 -----------------------------------------------------------------------
submission_folder <- "./submission_data/final_round/SC1"
bootstrap_RMSE <- read_rds("./submission_analysis/intermediate_data/sc1_bootstrap_rmse.rds")
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc1_ranked_teams.rds")

data_usage <- read_csv("./submission_data/final_round/data_usage.csv") %>% 
	filter(challenge==1) %>% 
	mutate(submitterId = make.names(submitterId))

# read leaderboard
SC_leaderboard = read_csv(
	file.path(submission_folder,"leaderboard_final_sc1.csv")) %>%
	select(-writeUp, -createdOn) %>% 
	mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
	mutate(submitterId = make.names(submitterId))

# Team KAUST_RSS didnt follow the template but managed to pass the validation 
SC_leaderboard <- SC_leaderboard %>% filter(submitterId != "KAUST_RSS") 

# Plot the performance of the team on the bootstrap samples
# the magic numbers come from script sc_prediction_anchors.R
# 
rank_plot <- bootstrap_RMSE %>% 
	gather(teams,RMSE,-BS_sample) %>%
	mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
	ggplot() + geom_boxplot(aes(teams, RMSE,color=teams)) + 
	geom_hline(yintercept = c(0.74, 0.91,1.44,2.21)) +
	geom_text(data = tibble(x_coord=c(1,1,1,1),
							y_coord=c(0.74,1,1.44,2.21),
							data_type=c("best in condition prediction","reference","shuffle by condition","shuffle all")),
			  aes(x_coord,y_coord,label=data_type), hjust = 0,vjust = 1 )+ 
	ggtitle("Subchallenge I: performance over bootstrap samples") + 
	geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) +
	theme_bw() +
	guides(color="none")+ 
	theme(axis.text.x = element_blank()) + labs(x=NULL)

data_usage_plot <- data_usage %>% select(-score,-challenge) %>%
	mutate_each(~ifelse(is.na(.),0,.)) %>% 
	gather(data_type,in_use, -submitterId) %>%
	mutate(submitterId = factor(submitterId,levels = levels(ranked_teams))) %>%
	ggplot() + 
	geom_tile(aes(submitterId,data_type,fill=as.factor(in_use)),color="black") + 
	scale_fill_manual(values = c("1"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],"0"="white")) +
	theme_bw() +
	xlab("teams") + 
	ylab("Data type used") + 
	theme(axis.text.x = element_text(hjust=1,angle = 45),legend.position = "none") 

plot_grid(rank_plot,data_usage_plot, axis = "lr", align = "v", rel_heights =  c(1.5, 1), ncol = 1)
ggsave("./submission_analysis/figures/sc1_rank_data_anchors.pdf", width = 6,height = 5)

### SC2 -----------------------------------------------------------------------
subchallange <- 2
submission_folder <- "./submission_data/final_round/SC2"
bootstrap_RMSE <- read_rds("./submission_analysis/intermediate_data/sc2_bootstrap_stats.rds")
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc2_ranked_teams.rds")

data_usage <- read_csv("./submission_data/final_round/data_usage.csv") %>% 
	filter(challenge==subchallange) %>% 
	mutate(submitterId = make.names(submitterId))

# read leaderboard
SC_leaderboard = read_csv(
	file.path(submission_folder,"leaderboard_final_sc2.csv")) %>%
	select(-writeUp, -createdOn) %>% 
	mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
	mutate(submitterId = make.names(submitterId))

# Plot the performance of the team on the bootstrap samples
rank_plot <- bootstrap_RMSE %>% 
	gather(teams,RMSE,-BS_sample) %>%
	mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
	ggplot() + geom_boxplot(aes(teams, RMSE,color=teams)) + 
	ggtitle("Subchallenge II: performance over bootstrap samples") + 
	geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) +
	theme_bw() +
	guides(color="none")+ 
	theme(axis.text.x = element_blank()) + labs(x=NULL)

data_usage_plot <- data_usage %>% select(-score,-challenge) %>%
	mutate_each(~ifelse(is.na(.),0,.)) %>% 
	gather(data_type,in_use, -submitterId) %>%
	mutate(submitterId = factor(submitterId,levels = levels(ranked_teams))) %>%
	ggplot() + 
	geom_tile(aes(submitterId,data_type,fill=as.factor(in_use)),color="black") + 
	scale_fill_manual(values = c("1"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],"0"="white")) +
	theme_bw() +
	xlab("teams") + 
	ylab("Data type used") + 
	theme(axis.text.x = element_text(hjust=1,angle = 45),legend.position = "none") 

plot_grid(rank_plot,data_usage_plot, axis = "lr", align = "v", rel_heights =  c(1.5, 1), ncol = 1)
ggsave("./submission_analysis/figures/sc2_rank_data.pdf", width = 6,height = 5)

### SC3 -----------------------------------------------------------------------

subchallange <- 3
submission_folder <- "./submission_data/final_round/SC3"
bootstrap_RMSE <- read_rds("./submission_analysis/intermediate_data/sc3_bootstrap_stats.rds")
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc3_ranked_teams.rds")

data_usage <- read_csv("./submission_data/final_round/data_usage.csv") %>% 
	filter(challenge==subchallange) %>% 
	mutate(submitterId = make.names(submitterId))

# read leaderboard
SC_leaderboard = read_csv(
	file.path(submission_folder,"leaderboard_final_sc3.csv")) %>%
	select(-writeUp, -createdOn) %>% 
	mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
	mutate(submitterId = make.names(submitterId))

# Plot the performance of the team on the bootstrap samples
rank_plot <- bootstrap_RMSE %>% 
	gather(teams,RMSE,-BS_sample) %>%
	mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
	ggplot() + geom_boxplot(aes(teams, RMSE,color=teams)) + 
	ggtitle("Subchallenge III: performance over bootstrap samples") + 
	geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) +
	theme_bw() +
	guides(color="none")+
	theme(axis.text.x = element_blank()) + labs(x=NULL)

data_usage_plot <- data_usage %>% select(-score,-challenge) %>%
	mutate_each(~ifelse(is.na(.),0,.)) %>% 
	gather(data_type,in_use, -submitterId) %>%
	mutate(submitterId = factor(submitterId,levels = levels(ranked_teams))) %>%
	ggplot() + 
	geom_tile(aes(submitterId,data_type,fill=as.factor(in_use)),color="black") + 
	scale_fill_manual(values = c("1"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],"0"="white")) +
	theme_bw() +
	xlab("teams") + 
	ylab("Data type used") + 
	theme(axis.text.x = element_text(hjust=1,angle = 45),legend.position = "none") 

plot_grid(rank_plot,data_usage_plot, axis = "lr", align = "v", rel_heights =  c(1.5, 1), ncol = 1)
ggsave("./submission_analysis/figures/sc3_rank_data.pdf", width = 6,height = 5)
### SC4 -----------------------------------------------------------------------

subchallange <- 4
submission_folder <- "./submission_data/final_round/SC4"
bootstrap_RMSE <- read_rds("./submission_analysis/intermediate_data/sc4_bootstrap_RMSE.rds")
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc4_ranked_teams.rds")

data_usage <- read_csv("./submission_data/final_round/data_usage.csv") %>% 
	filter(challenge==subchallange) %>% 
	mutate(submitterId = make.names(submitterId))

# read leaderboard
SC_leaderboard = read_csv(
	file.path(submission_folder,"leaderboard_final_sc4.csv")) %>%
	select(-writeUp, -createdOn) %>% 
	mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
	mutate(submitterId = make.names(submitterId))

# Plot the performance of the team on the bootstrap samples
rank_plot <- bootstrap_RMSE %>% 
	gather(teams,RMSE,-BS_sample) %>%
	mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
	ggplot() + geom_boxplot(aes(teams, RMSE,color=teams)) + 
	ggtitle("Subchallenge IV: performance over bootstrap samples") + 
	geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) +
	theme_bw() +
	guides(color="none")+
	theme(axis.text.x = element_blank()) + labs(x=NULL)

data_usage_plot <- data_usage %>% select(-score,-challenge) %>%
	mutate_each(~ifelse(is.na(.),0,.)) %>% 
	gather(data_type,in_use, -submitterId) %>%
	mutate(submitterId = factor(submitterId,levels = levels(ranked_teams))) %>%
	ggplot() + 
	geom_tile(aes(submitterId,data_type,fill=as.factor(in_use)),color="black") + 
	scale_fill_manual(values = c("1"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],"0"="white")) +
	theme_bw() +
	xlab("teams") + 
	ylab("Data type used") + 
	theme(axis.text.x = element_text(hjust=1,angle = 45),legend.position = "none") 

plot_grid(rank_plot,data_usage_plot, axis = "lr", align = "v", rel_heights =  c(1.5, 1), ncol = 1)
ggsave("./submission_analysis/figures/sc4_rank_data.pdf", width = 6,height = 5)

# Plot the Bootstrap results together with the data usage heatmap 

library(tidyverse)
library(cowplot)

### SC1 ------------------------------------------------------------------------
submission_folder <- "./submission_data/final_round/SC1"
bootstrap_RMSE <- read_rds("./submission_analysis/intermediate_data/sc1_bootstrap_rmse.rds")
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc1_ranked_teams.rds")

data_usage <- read_csv("./submission_data/final_round/data_usage.csv") %>% 
	filter(challenge==1) %>% 
	mutate(submitterId = make.names(submitterId))

method_usage <- readxl::read_xlsx("./methods/sc1_methods_processed.xlsx",sheet = 2)
descriptor_order <- method_usage$Question

method_usage <-  method_usage %>%
	gather(submitterId,answer,-Question) %>% spread(Question,answer)

# we remove the Anand_1812 from this analysis: 
method_usage <- method_usage[-which(method_usage$submitterId =="anand_1812/CSBL"),]


# read leaderboard
SC_leaderboard = read_csv(
	file.path(submission_folder,"leaderboard_final_sc1.csv")) %>%
	select(-writeUp, -createdOn) %>% 
	mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
	mutate(submitterId = make.names(submitterId))

# Team KAUST_RSS didnt follow the template but managed to pass the validation 
SC_leaderboard <- SC_leaderboard %>% filter(submitterId != "KAUST_RSS") 


# fix naming
ranked_teams <- factor(gsub("X.","",ranked_teams,fixed = T),levels = gsub("X.","",ranked_teams,fixed = T))
names(bootstrap_RMSE) <- gsub("X.","",names(bootstrap_RMSE),fixed = T)
data_usage$submitterId <- gsub("X.","",data_usage$submitterId,fixed = T)


# Plot the performance of the team on the bootstrap samples
# the magic numbers come from script sc_prediction_anchors.R
# 

scores_ordered_teams <- read_rds("./submission_data/combination_scoring/SC1/SC1_ordered_combination.rds")
scores_random_teams <- read_rds("./submission_data/combination_scoring/SC1/SC1_random_subs_scores.rds")
RF_combination <- read_rds("./submission_data/combination_scoring/SC1/LOO_CV_asvRF_incl32_scores.rds")


sum_scores_random_teams <- scores_random_teams %>%  group_by(Sample_size) %>%
	summarise(score_min = quantile(Median,0,25),
			  score_max = quantile(Median,0.75),
			 score_med = median(Median)) 

ref_numbers <- c(ordered_teams = scores_ordered_teams,
				 random_teams = 0.857,
				 asv_RF_32 = 0.877,
				 reference_model= 0.91,
				 random = 1.44) %>% enframe(value = "score") %>%
	mutate(x = 1)


my_colors  = RColorBrewer::brewer.pal(8,"Dark2")


rank_plot <- bootstrap_RMSE %>% 
	gather(teams,RMSE,-BS_sample) %>%
	group_by(teams) %>% summarise(score_min = min(RMSE),
								  score_max = max(RMSE),
								  score_med = median(RMSE)) %>% 
	ungroup() %>%
	mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
	ggplot() +
	# individual performance:
	#geom_errorbar(aes(teams,ymin = score_min,ymax=score_max) ,color = my_colors[[1]]) + 
	geom_point(aes(teams,score_med), color = my_colors[[1]]) + 
	#geom_line(aes(teams,score_med,group=1), color = my_colors[[1]]) + 
	
    # combination of ordered predictions:
	geom_point(data =scores_ordered_teams, aes(n,Median),color = my_colors[[2]] ) +
	geom_line(data =scores_ordered_teams, aes(n,Median),color = my_colors[[2]] ) +
	
    # combination of random teams
	geom_errorbar(data=sum_scores_random_teams, aes(Sample_size,ymin = score_min,ymax=score_max) ,color = my_colors[[3]]) + 
	geom_point(data=sum_scores_random_teams, aes(Sample_size,score_med), color = my_colors[[3]]) + 
	geom_line(data=sum_scores_random_teams, aes(Sample_size,score_med), color = my_colors[[3]]) + 
	
    #geom_boxplot(aes(teams, RMSE, group=teams),size=0.2, outlier.size = 0.5) + 
	#geom_hline(yintercept = c(0.74,0.772, 0.91,1.44,2.21)) +
	#geom_text(data = tibble(x_coord=c(1,1,1,1,1),
#							y_coord=c(0.74,0.772,1,1.44,2.21),
#							data_type=c("best in condition prediction","RF all CL","reference","shuffle by condition","shuffle all")),
#			  aes(x_coord,y_coord,label=data_type), hjust = 0,vjust = 1 )+
	geom_hline(yintercept = mean(RF_combination$asv_RF_32 ))+
	#ggtitle("Subchallenge I: performance over bootstrap samples") + 
	#geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) +
	theme_bw() +
	guides(color="none")+ 
	theme(axis.text.x = element_blank(),
		  axis.ticks = element_blank()) + labs(x=NULL) + 
	coord_cartesian(ylim = c(0.84,1)) + ylab("Score (RMSE)") 
print(rank_plot)

# 
# data_usage_plot <- data_usage %>% select(-score,-challenge) %>%
# 	mutate_each(~ifelse(is.na(.),0,.)) %>% 
# 	gather(data_type,in_use, -submitterId) %>%
# 	mutate(submitterId = factor(submitterId,levels = levels(ranked_teams))) %>%
# 	ggplot() + 
# 	geom_tile(aes(submitterId,data_type,fill=as.factor(in_use)),color="black") + 
# 	scale_fill_manual(values = c("1"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],"0"="white")) +
# 	theme_bw() +
# 	xlab("teams") + 
# 	ylab("Data type used") + 
# 	theme(axis.text.x = element_text(hjust=1,angle = 45),legend.position = "none")  + 
# 	theme(axis.text.x = element_blank()) + labs(x=NULL)
# 

method_usage_plot <- method_usage %>% select(-5,-11,-12) %>%
	mutate(submitterId = factor(submitterId,levels = levels(ranked_teams))) %>%
	gather(method_info,method_value,-submitterId) %>%
	mutate(method_info = factor(method_info,levels = rev(descriptor_order))) %>%
	ggplot() + 
	geom_tile(aes(submitterId,method_info,fill=as.factor(method_value)),color="black") + 
	scale_fill_manual(values = c("y"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],"n"="white")) +
	theme_bw() +
	xlab("Teams") + 
	ylab("Methods") + 
	theme(axis.text.x = element_blank(),
		  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		  panel.border = element_blank(),
		  axis.ticks = element_blank(),legend.position = "none") + labs(x=NULL) 

data_usage_plot <- method_usage %>% select(1,5,11,12) %>% 
	mutate(submitterId = factor(submitterId,levels = levels(ranked_teams))) %>%
	gather(method_info,method_value,-submitterId) %>%
	mutate(method_info = factor(method_info,levels = rev(descriptor_order))) %>%
	ggplot() + 
	geom_tile(aes(submitterId,method_info,fill=as.factor(method_value)),color="black") + 
	scale_fill_manual(values = c("y"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],"n"="white")) +
	theme_bw() +
	xlab("Teams") + 
	ylab("Data") + 
	theme(axis.text.x = element_text(hjust=1,angle = 45),
		  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		  axis.ticks = element_blank(),
		  panel.border = element_blank(),legend.position = "none") 


plot_grid(rank_plot,
		  method_usage_plot, 
		  data_usage_plot,
		  axis = "lr", align = "v", rel_heights =  c(1.5, 1,1.2), ncol = 1)
ggsave("./submission_analysis/figures/sc1_rank_data_anchors.pdf", width = 6,height = 5)


ref_numbers <- c(ordered_teams = scores_ordered_teams,
					  random_teams = 0.857,
					  asv_RF_32 = 0.877,
					  reference_model= 0.91,
					  random = 1.44) %>% enframe(value = "score") %>%
	mutate(x = 1)
					  

library(ggrepel)
ref_models <- ref_numbers  %>% ggplot() + geom_point(aes(x,score)) + 
	geom_hline(aes(yintercept = score)) +
	#geom_text_repel(aes(x,score,label=data_type),xlim = c(1.025,NA)) +
	theme_bw() +
	#ylim(0.7,2.25,)+
	theme(axis.text.x = element_blank()) + labs(x=NULL) +
	scale_y_continuous(position = "right",limits =c(0.7,2.25) )

ggsave("./submission_analysis/figures/sc1_ref_models_score.pdf", width = 1,height = 2)

#plot_grid(rank_plot,ref_models,method_usage_plot, axis = "lr", align = "v", rel_heights =  c(1.5, 1), ncol = 2)


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

# DREAM challenge post analysis
# SC4: compare the teams and combination of predictions, data-usage etc
# A. Gabor

library(tidyverse)
library(cowplot)



### SC4 -----------------------------------------------------------------------
subchallange <- 4
submission_folder <- "./submission_data/final_round/SC4"
bootstrap_RMSE <- read_rds("./submission_analysis/intermediate_data/sc4_bootstrap_rmse.rds")
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc4_ranked_teams.rds")



method_usage <- readxl::read_xlsx("./methods/sc4_methods_processed.xlsx",sheet = 2)
descriptor_order <- method_usage$Question

method_usage <-  method_usage %>%
    gather(submitterId,answer,-Question) %>% spread(Question,answer)

# we remove the Anand_1812 from this analysis: 
#method_usage <- method_usage[-which(method_usage$submitterId =="anand_1812/CSBL"),]


# read leaderboard
SC_leaderboard = read_csv(
    file.path(submission_folder,"leaderboard_final_sc4.csv")) %>%
    select(-writeUp, -createdOn) %>% 
    mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
    mutate(submitterId = make.names(submitterId))

# Team KAUST_RSS didnt follow the template but managed to pass the validation 
#SC_leaderboard <- SC_leaderboard %>% filter(submitterId != "KAUST_RSS") 


# fix naming
ranked_teams <- factor(gsub("X.","",ranked_teams,fixed = T),levels = gsub("X.","",ranked_teams,fixed = T))
names(bootstrap_RMSE) <- gsub("X.","",names(bootstrap_RMSE),fixed = T)



# Plot the performance of the team on the bootstrap samples

scores_ordered_teams <- read_rds("./submission_data/combination_scoring/SC4/SC4_ordered_combination.rds")
scores_random_teams <- read_rds("./submission_data/combination_scoring/SC4/SC4_random_subs_scores.rds")
RF_combination <- read_rds("./submission_data/combination_scoring/SC4/LOO_CV_RF_scores.rds")


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
    geom_point(aes(teams,score_med, color = "team prediction")) + 
    #geom_line(aes(teams,score_med,group=1), color = my_colors[[1]]) + 
    
    # combination of ordered predictions:
    geom_point(data =scores_ordered_teams, aes(n,Median,color = "combination top N") ) +
    geom_line(data =scores_ordered_teams, aes(n,Median,color = "combination top N")) +
    
    # combination of random teams
    geom_errorbar(data=sum_scores_random_teams, aes(Sample_size,ymin = score_min,ymax=score_max, color ="random combination" ) ) + 
    geom_point(data=sum_scores_random_teams, aes(Sample_size,score_med), color = my_colors[[3]]) + 
    geom_line(data=sum_scores_random_teams, aes(Sample_size,score_med,color = "random combination")) + 
    #geom_boxplot(aes(teams, RMSE, group=teams),size=0.2, outlier.size = 0.5) + 
    #geom_hline(yintercept = c(0.74,0.772, 0.91,1.44,2.21)) +
    #geom_text(data = tibble(x_coord=c(1,1,1,1,1),
    #							y_coord=c(0.74,0.772,1,1.44,2.21),
    #							data_type=c("best in condition prediction","RF all CL","reference","shuffle by condition","shuffle all")),
    #			  aes(x_coord,y_coord,label=data_type), hjust = 0,vjust = 1 )+
    geom_hline(data = NULL, aes(yintercept = mean(RF_combination$lm), color = "RF combination"))+
    geom_hline(data = NULL,  aes(yintercept = 0.428, color = "average CL"))+
    #ggtitle("Subchallenge I: performance over bootstrap samples") + 
    #geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "left",
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0),
          #legend.background = element_blank(),
          legend.box.margin =  margin(0, 0, 0, 0),
          legend.spacing.y = unit(0.05, 'cm')) +
    labs(x=NULL) + 
    #coord_cartesian(ylim = c(0,0.7)) +
    ylab("Score (RMSE)") +
    scale_color_manual(values = c("team prediction"= my_colors[[1]],
                                  "combination top N" = my_colors[[2]],
                                  "random combination" = my_colors[[3]],
                                  "RF combination" = my_colors[[5]],
                                  "average CL" = my_colors[[4]]))
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

method_usage_plot <- method_usage %>% 
    select(-2,-6,-8,-11,-12,-13,-14) %>%  # remove data columns
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


hide_names <- tibble(submitterId = ranked_teams, alt_name = as.character(1:length(ranked_teams)))
hide_names$alt_name = factor(hide_names$alt_name,levels =hide_names$alt_name)



data_usage_plot <- method_usage %>% 
    select(1,2,6,8,11,12,13,14) %>%  # remove data columns
    mutate(submitterId = factor(submitterId,levels = levels(ranked_teams))) %>%
    gather(method_info,method_value,-submitterId) %>%
    mutate(method_info = factor(method_info,levels = rev(descriptor_order))) %>%
    left_join(hide_names,by="submitterId") %>% 
    select(-submitterId) %>%
    ggplot() + 
    geom_tile(aes(alt_name,method_info,fill=as.factor(method_value)),color="black") + 
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
          axis = "lr", align = "v", rel_heights =  c(1.5, 0.65, 0.9), ncol = 1)
ggsave("./publication/figures/figure5/sc4_leaderboard_anonym.pdf", width = 6,height = 5)





### Simplified leaderboard


hide_names <- tibble(teams = ranked_teams, alt_name = as.character(1:length(ranked_teams)))
hide_names$alt_name = factor(hide_names$alt_name,levels =hide_names$alt_name)



rank_plot <- bootstrap_RMSE %>% 
    gather(teams,RMSE,-BS_sample) %>%
    mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
    left_join(hide_names,by="teams") %>% 
    select(-teams) %>%
    ggplot() +
    geom_violin(aes(alt_name,RMSE,fill="team prediction"),draw_quantiles = 0.5) + 
    geom_hline(data = NULL,  aes(yintercept = 0.428, color = "average cell line"))+
    geom_hline(data = NULL,  aes(yintercept = 0.3408274, color = "reference model"))+
    theme_bw() +
    theme(# axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.2,0.6),
        legend.title = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        #legend.background = element_blank(),
        legend.box.margin =  margin(0, 0, 0, 0),
        legend.spacing.y = unit(0.05, 'cm')) +
    xlab("Teams") + 
    ylab("Score (RMSE)") +
    scale_fill_manual(values = c("team prediction"= my_colors[[1]])) +
    scale_color_manual(values = c("average cell line" = my_colors[[4]],
                                  "reference model"=my_colors[[3]]))
print(rank_plot)
ggsave("./publication/figures/figure5/sc4_leaderboard_anonym_rank_only_nocombination.pdf", width = 6,height = 3)

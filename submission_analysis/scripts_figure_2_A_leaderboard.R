# scirpts for figure 2:

library(tidyverse)
library(cowplot)

### A) Leaderboard overview ------------------------------------------------------

submission_folder <- "./submission_data/final_round/SC1"
bootstrap_RMSE <- read_rds("./submission_analysis/intermediate_data/sc1_bootstrap_rmse.rds")
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc1_ranked_teams.rds")
# ranked_teams2 <- read_rds("./submission_analysis/intermediate_data/sc2_ranked_teams.rds")
# ranked_teams3 <- read_rds("./submission_analysis/intermediate_data/sc3_ranked_teams.rds")
# ranked_teams4 <- read_rds("./submission_analysis/intermediate_data/sc4_ranked_teams.rds")

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



# Plot the performance of the team on the bootstrap samples

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



rank_plot <- 
    # uncomment if for points: 
    # bootstrap_RMSE %>% 
    # gather(teams,RMSE,-BS_sample) %>%
    # group_by(teams) %>% summarise(score_min = min(RMSE),
    #                               score_max = max(RMSE),
    #                               score_med = median(RMSE)) %>% 
    # ungroup() %>%
    # mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
    # ggplot() +
    # # individual performance:
    # geom_point(aes(teams,score_med, color = "team prediction")) + 

# violin plot version: 
bootstrap_RMSE %>% 
    gather(teams,RMSE,-BS_sample) %>%
    mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
    ggplot() +
    # individual performance:
    geom_violin(aes(teams,RMSE,color ="team prediction", fill = "team prediction"),draw_quantiles = 0.5) + 
    
    
    #geom_line(aes(teams,score_med,group=1), color = my_colors[[1]]) + 
    
    # combination of ordered predictions:
    geom_point(data =scores_ordered_teams, aes(n,Median,color = "combination top N") ) +
    geom_line(data =scores_ordered_teams, aes(n,Median,color = "combination top N")) +
    
    # combination of random teams
    geom_errorbar(data=sum_scores_random_teams, aes(Sample_size,ymin = score_min,ymax=score_max, color ="random combination" ) ) + 
    geom_point(data=sum_scores_random_teams, aes(Sample_size,score_med), color = my_colors[[3]]) + 
    geom_line(data=sum_scores_random_teams, aes(Sample_size,score_med,color = "random combination")) + 
    #geom_boxplot(aes(teams, RMSE, group=teams),size=0.2, outlier.size = 0.5) + 
    geom_hline(aes(yintercept = c(1.44),color="random")) +
    #geom_hline(yintercept = c(0.74,0.772, 0.91,1.44,2.21)) +
    #geom_text(data = tibble(x_coord=c(1,1,1,1,1),
    #							y_coord=c(0.74,0.772,1,1.44,2.21),
    #							data_type=c("best in condition prediction","RF all CL","reference","shuffle by condition","shuffle all")),
    #			  aes(x_coord,y_coord,label=data_type), hjust = 0,vjust = 1 )+
    geom_hline(data = NULL, aes(yintercept = mean(RF_combination$asv_RF_32), color = "RF combination"))+
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
    coord_cartesian(ylim = c(0.84,1.5)) + ylab("Score (RMSE)") +
    scale_color_manual(values = c("team prediction"= my_colors[[1]],
                                  "combination top N" = my_colors[[2]],
                                  "random combination" = my_colors[[3]],
                                  "RF combination" = my_colors[[5]],
                                  "random" = my_colors[[4]])) +
    scale_fill_manual(values = c("team prediction"= my_colors[[1]])) 
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

method_usage_plot <- method_usage %>% select(-5,-10,-11,-12,-13) %>%
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



data_usage_plot <- method_usage %>% select(1,5,11,12,13) %>% 
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
          axis = "lr", align = "v", rel_heights =  c(1.5, 1, 0.7), ncol = 1)
ggsave("./publication/figures/figure2/sc1_leaderboard_anonym.pdf", width = 6,height = 5)


########## Let's try to separate the ranking an combination predictions 



### leaderboard plot, teams by number, violin plot, with combinations
hide_names <- tibble(teams = ranked_teams, alt_name = as.character(1:length(ranked_teams)))
hide_names$alt_name = factor(hide_names$alt_name,levels =hide_names$alt_name)

rank_plot <- 
    # uncomment if for points: 
    # bootstrap_RMSE %>% 
    # gather(teams,RMSE,-BS_sample) %>%
    # group_by(teams) %>% summarise(score_min = min(RMSE),
    #                               score_max = max(RMSE),
    #                               score_med = median(RMSE)) %>% 
    # ungroup() %>%
    # mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
    # ggplot() +
    # # individual performance:
    # geom_point(aes(teams,score_med, color = "team prediction")) + 

# violin plot version: 
bootstrap_RMSE %>% 
    gather(teams,RMSE,-BS_sample) %>%
    mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
    left_join(hide_names,by="teams") %>% 
    select(-teams) %>%
    ggplot() +
    # individual performance:
    geom_violin(aes(alt_name,RMSE,color ="team prediction", fill = "team prediction"),draw_quantiles = 0.5) + 
    #geom_line(aes(teams,score_med,group=1), color = my_colors[[1]]) + 
    
    # combination of ordered predictions:
    geom_point(data =scores_ordered_teams, aes(n,Median,color = "combination top N") ) +
    geom_line(data =scores_ordered_teams, aes(n,Median,color = "combination top N")) +
    
    # combination of random teams
    geom_errorbar(data=sum_scores_random_teams, aes(Sample_size,ymin = score_min,ymax=score_max, color ="random combination" ) ) + 
    geom_point(data=sum_scores_random_teams, aes(Sample_size,score_med), color = my_colors[[3]]) + 
    geom_line(data=sum_scores_random_teams, aes(Sample_size,score_med,color = "random combination")) + 
    #geom_boxplot(aes(teams, RMSE, group=teams),size=0.2, outlier.size = 0.5) + 
    geom_hline(aes(yintercept = c(1.44),color="random")) +
    #geom_hline(yintercept = c(0.74,0.772, 0.91,1.44,2.21)) +
    #geom_text(data = tibble(x_coord=c(1,1,1,1,1),
    #							y_coord=c(0.74,0.772,1,1.44,2.21),
    #							data_type=c("best in condition prediction","RF all CL","reference","shuffle by condition","shuffle all")),
    #			  aes(x_coord,y_coord,label=data_type), hjust = 0,vjust = 1 )+
    geom_hline(data = NULL, aes(yintercept = mean(RF_combination$asv_RF_32), color = "RF combination"))+
    #ggtitle("Subchallenge I: performance over bootstrap samples") + 
    #geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        #legend.position = "left",
        legend.title = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        #legend.background = element_blank(),
        legend.box.margin =  margin(0, 0, 0, 0),
        legend.spacing.y = unit(0.05, 'cm'),
        legend.position = c(0.2,0.6)
        ) +
    labs(x="Teams") + 
    coord_cartesian(ylim = c(0.84,1.5)) + ylab("Score (RMSE)") +
    scale_color_manual(values = c("team prediction"= my_colors[[1]],
                                  "combination top N" = my_colors[[2]],
                                  "random combination" = my_colors[[3]],
                                  "RF combination" = my_colors[[5]],
                                  "random" = my_colors[[4]])) +
    scale_fill_manual(values = c("team prediction"= my_colors[[1]])) 
print(rank_plot)

ggsave("./publication/figures/figure2/sc1_leaderboard_anonym_rank_only.pdf", width = 6,height = 3)




### leaderboard plot, teams by number, violin plot, no combinations
hide_names <- tibble(teams = ranked_teams, alt_name = as.character(1:length(ranked_teams)))
hide_names$alt_name = factor(hide_names$alt_name,levels =hide_names$alt_name)

rank_plot <- 
    # violin plot version: 
    bootstrap_RMSE %>% 
    gather(teams,RMSE,-BS_sample) %>%
    mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
    left_join(hide_names,by="teams") %>% 
    select(-teams) %>%
    ggplot() +
    # individual performance:
    geom_violin(aes(alt_name,RMSE, fill = "team prediction"),draw_quantiles = 0.5) + 
    geom_hline(aes(yintercept = c(1.44),color="random")) +
    geom_hline(aes(yintercept = c(0.903103),color="reference")) +
    
    theme_bw() +
    theme(panel.grid = element_blank(),
          #legend.position = "left",
          legend.position = c(0.2,0.6),
          legend.title = element_blank(),
          #legend.margin = margin(0, 0, 0, 0),
          #legend.box.margin =  margin(0, 0, 0, 0),
          legend.spacing.y = unit(0.05, 'cm')
    ) +
    xlab("Teams") + 
    coord_cartesian(ylim = c(0.84,1.5)) + ylab("Score (RMSE)") +
    scale_color_manual(values = c("team prediction"= my_colors[[1]],
                                  "random" = my_colors[[4]],
                                  'reference' = my_colors[[3]])) +
    scale_fill_manual(values = c("team prediction"= my_colors[[1]])) 
print(rank_plot)

ggsave("./publication/figures/figure2/sc1_leaderboard_anonym_rank_only_nocombination.pdf", width = 6,height = 3)



method_usage_plot <- method_usage %>% select(-5,-8,-11,-12) %>%
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



data_usage_plot <- method_usage %>% select(1,5,11,12,13) %>% 
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
          axis = "lr", align = "v", rel_heights =  c(1.5, 1, 0.7), ncol = 1)




### Plot only the combinations:



ggplot() +
    # combination of ordered predictions:
    geom_point(data =scores_ordered_teams, aes(n,Median,color = "combination top N") ) +
    geom_line(data =scores_ordered_teams, aes(n,Median,color = "combination top N")) +
    
    # combination of random teams
    geom_errorbar(data=sum_scores_random_teams, aes(Sample_size,ymin = score_min,ymax=score_max, color ="random combination" ) ) + 
    geom_point(data=sum_scores_random_teams, aes(Sample_size,score_med), color = my_colors[[3]]) + 
    geom_line(data=sum_scores_random_teams, aes(Sample_size,score_med,color = "random combination")) + 
    #geom_boxplot(aes(teams, RMSE, group=teams),size=0.2, outlier.size = 0.5) + 
    #geom_hline(aes(yintercept = c(1.44),color="random")) +
    geom_hline(aes(yintercept = SC_leaderboard$score[[1]], color="best team")) +
    #geom_text(data = tibble(x_coord=c(1,1,1,1,1),
    #							y_coord=c(0.74,0.772,1,1.44,2.21),
    #							data_type=c("best in condition prediction","RF all CL","reference","shuffle by condition","shuffle all")),
    #			  aes(x_coord,y_coord,label=data_type), hjust = 0,vjust = 1 )+
    geom_hline(data = NULL, aes(yintercept = mean(RF_combination$asv_RF_32), color = "RF combination"))+
    #ggtitle("Subchallenge I: performance over bootstrap samples") + 
    #geom_point(data = SC_leaderboard %>% filter(submitterId %in% ranked_teams), aes(submitterId,score)) +
    theme_bw() +
    theme(legend.position = "left",
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0),
          #legend.background = element_blank(),
          legend.box.margin =  margin(0, 0, 0, 0),
          legend.spacing.y = unit(0.05, 'cm')) +
    xlab("Number of combined predictions") + 
    #coord_cartesian(ylim = c(0.84,1.5)) +
    ylab("Score (RMSE)") +
    scale_color_manual(values = c("best team"= my_colors[[1]],
                                  "combination top N" = my_colors[[2]],
                                  "random combination" = my_colors[[3]],
                                  "RF combination" = my_colors[[5]],
                                  "random" = my_colors[[4]])) +
    scale_fill_manual(values = c("team prediction"= my_colors[[1]])) 

ggsave("./publication/figures/figure2/sc1_teams_combinations_score.pdf", width = 6,height = 3)


best_team = tibble(methods = "best_team",score = SC_leaderboard$score[[1]])
RF        = tibble(methods = "RF",score = mean(RF_combination$asv_RF_32))

random_teams = tibble(methods = "random_combination_team",score = median(scores_random_teams$Median),
                      min = min(scores_random_teams$Median),
                      max= max(scores_random_teams$Median))

ordered_teams = tibble(methods = "ordered_combination_team",score = median(scores_ordered_teams$Median),
                       min = min(scores_ordered_teams$Median),
                       max= max(scores_ordered_teams$Median))

bind_rows(best_team,RF,random_teams,ordered_teams) %>%
    ggplot(aes(methods,score)) + geom_point() + geom_errorbar(aes(ymin=min,ymax=max))


###

scores_random_teams %>%
    group_by(Sample_size) %>%
    summarise(N_score = mean(Median)) %>%
    arrange(N_score) 
scores_ordered_teams$Median[[2]] 
random_teams
ordered_teams

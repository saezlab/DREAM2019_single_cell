# DREAM challenge post analysis
# SC2: compare the teams and combination of predictions, data-usage etc
# A. Gabor

library(tidyverse)
library(cowplot)



### SC2 -----------------------------------------------------------------------
subchallange <- 2
submission_folder <- "./submission_data/final_round/SC2"
bootstrap_RMSE <- read_rds("./submission_analysis/intermediate_data/sc2_bootstrap_stats.rds")
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc2_ranked_teams.rds")

conditional_error_stats <- read_rds("./submission_analysis/intermediate_data/sc2_stats_sumSquared_conditions.rds")


# load prediction combinations

#Combination by sampling from random teams:
#scores_ordered_teams <- read_rds("./submission_data/combination_scoring/SC2/SC2_ordered_combination.rds")
#scores_random_teams <- read_rds("./submission_data/combination_scoring/SC2/SC2_random_subs_scores.rds")


#Combination by stats from random teams:
scores_ordered_teams <- read_rds("./prediction_combinations/SC2/SC2_ordered_combination_stats.rds")
scores_random_teams <- read_rds("./prediction_combinations/SC2/SC2_random_subs_scores_stats.rds")


RF_combination <- read_rds("./submission_data/combination_scoring/SC2/L2O_CV_scores.rds") %>%
    summarise_at(3:9,~sqrt(sum(.^2)))


method_usage <- readxl::read_xlsx("./methods/sc2_methods_processed.xlsx",sheet = 2)
descriptor_order <- method_usage$Question

method_usage <-  method_usage %>%
    gather(submitterId,answer,-Question) %>% spread(Question,answer)


# read leaderboard
SC_leaderboard = read_csv(
    file.path(submission_folder,"leaderboard_final_sc2.csv")) %>%
    select(-writeUp, -createdOn) %>% 
    mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
    mutate(submitterId = make.names(submitterId))


# fix naming
ranked_teams <- factor(gsub("X.","",ranked_teams,fixed = T),levels = gsub("X.","",ranked_teams,fixed = T))
names(bootstrap_RMSE) <- gsub("X.","",names(bootstrap_RMSE),fixed = T)
names(conditional_error_stats)<- gsub("X.","",names(conditional_error_stats),fixed = T)

# Basic statistics:
conditional_error_stats %>% 
    gather(teams,cond_score,-cond_id) %>% 
    separate(cond_id,into = c("cell_line","treatment","time"),sep = "_")%>%
    group_by(teams) %>% 
    summarise(score = sqrt(sum(cond_score))) %>% 
    ungroup()  %>% summarise(min_score = min(score),
                             max_score = max(score),
                             mean_score = mean(score),
                             median_score = median(score),
                             std_score = sd(score))

bootstrap_RMSE %>% 
    gather(teams,RMSE,-BS_sample) %>%
    group_by(teams) %>% summarise(score_min = min(RMSE),
                                  score_max = max(RMSE),
                                  score_med = median(RMSE),
                                  score_mean = mean(RMSE),
                                  score_std = sd(RMSE)) 
    
    
# process combinations:
sum_scores_random_teams <- scores_random_teams %>%  group_by(Sample_size) %>%
    summarise(score_min = quantile(Median,0,25),
              score_max = quantile(Median,0.75),
              score_med = median(Median)) 


my_colors  = RColorBrewer::brewer.pal(8,"Dark2")


# Plot the performance of the team on the bootstrap samples
rank_plot <- conditional_error_stats %>% 
    gather(teams,cond_score,-cond_id) %>% 
    separate(cond_id,into = c("cell_line","treatment","time"),sep = "_")%>%
    group_by(teams) %>% 
    summarise(score = sqrt(sum(cond_score))) %>% 
    ungroup() %>%
    mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
# bootstrap_RMSE %>% 
    # gather(teams,RMSE,-BS_sample) %>%
    # group_by(teams) %>% summarise(score_min = min(RMSE),
    #                               score_max = max(RMSE),
    #                               score_med = median(RMSE)) %>% 
    # ungroup() %>%
    # mutate(teams = factor(teams,levels = levels(ranked_teams))) %>%
    ggplot() +
    # individual performance:
    geom_point(aes(teams,score,col="team prediction")) + 
    
    # combination of ordered predictions:
    geom_point(data =scores_ordered_teams, aes(n,median_stats,color="combination top N")) +
    geom_line(data =scores_ordered_teams, aes(n,median_stats,color="combination top N")) +
    
    # combination of random teams
    geom_errorbar(data=sum_scores_random_teams, aes(Sample_size,ymin = score_min,ymax=score_max) ,color = my_colors[[3]]) + 
    geom_point(data=sum_scores_random_teams, aes(Sample_size,score_med,color = "random combination")) + 
    geom_line(data=sum_scores_random_teams, aes(Sample_size,score_med), color = my_colors[[3]]) + 
    geom_hline(data=NULL,aes(yintercept = mean(RF_combination$stats_RF), color = "RF combination"))+
    geom_hline(data=NULL,aes(yintercept = 29.1559, color = "EGF prediction"))+
    
    theme_bw() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        #legend.position = c(0.2, 0.85),
        legend.position = "left",
        legend.title = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        #legend.background = element_blank(),
        legend.box.margin =  margin(0, 0, 0, 0),
        legend.spacing.y = unit(0.05, 'cm')) +
    labs(x=NULL) + 
    coord_cartesian(ylim = c(20,100)) +
    ylab("Score (Median, Cov)") +
    scale_color_manual(values = c("team prediction"=my_colors[[1]],
                                  "combination top N" = my_colors[[2]],
                                  "random combination" = my_colors[[3]],
                                  "EGF prediction" = my_colors[[4]],
                                  "RF combination" = my_colors[[5]]))
    
print(rank_plot)




method_usage_plot <- method_usage %>% 
    select(-6,-12,-14, -15) %>% # remove data usage entities
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
          axis.ticks = element_blank(),legend.position = "none") +
    labs(x=NULL)  



hide_names <- tibble(submitterId = ranked_teams, alt_name = as.character(1:length(ranked_teams)))
hide_names$alt_name = factor(hide_names$alt_name,levels =hide_names$alt_name)

data_usage_plot <- method_usage %>% 
    select(1,6,12,14,15) %>%  # select ID and  data usage entities
    mutate(submitterId = factor(submitterId,levels = levels(ranked_teams))) %>%
    left_join(hide_names,by="submitterId") %>% 
    select(-submitterId) %>%
    gather(method_info,method_value,-alt_name) %>%
    mutate(method_info = factor(method_info,levels = rev(descriptor_order))) %>%
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
          axis = "lr", align = "v", rel_heights =  c(1.5, 1, .7), ncol = 1)
ggsave("./publication/figures/figure3/sc2_leaderboard_anonym.pdf", width = 6,height = 5)


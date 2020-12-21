# DREAM challenge post analysis
# SC2: compare the predictions with biological replica
# A. Gabor

library(tidyverse)
# SC2


# load data
submission_folder = "./submission_data/final_round/SC2/"
SC_leaderboard = read_csv(file.path(submission_folder,"leaderboard_final_sc2.csv")) %>%
    select(-writeUp, -createdOn) %>% 
    mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
    mutate(submitterId = make.names(submitterId))

gs <- read_csv("./challenge_data/validation_data/sc2gold.csv")
source("./scoring_scripts/score_sc2.R")

# fileID table maps the fileIDs with the biological replica A and B
fileID_table <- read_csv("./challenge_data/FileID_table.csv") %>% 
    select(cell_line,treatment,time,time_course,fileID)

ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc2_ranked_teams.rds")
combined_statistics <- read_rds("./submission_analysis/intermediate_data/sc2_stats_conditions.rds")

# stats for gold standard data:: (A and B are merged)
gs_stats_original <- gs %>% select(-fileID,-cellID) %>% data_to_stats(.) %>%
    arrange(cell_line, treatment,  time, stat_variable) %>%
    rename(standard_test=stat_value)


# find conditions, where there is both A and B condition measured
gs %>%	left_join(fileID_table, by = c("cell_line", "treatment", "time", "fileID")) %>%
    group_by(cell_line, treatment, time) %>% summarise(n_replica = length(unique(time_course))) %>% 
    filter(n_replica>1)

# these happens only at time 0 in inhibitory experiments
# for each time_course (A/B) we compute the stats using the data_to_stats function, 
# this function internally considers the grouping by treatment, time and cell-line
gs_stats <- gs %>% filter(time==0) %>%
    left_join(fileID_table, by = c("cell_line", "treatment", "time", "fileID")) %>%
    select(-fileID,-cellID) %>%
    group_by(time_course) %>% nest() %>%
    mutate(stats = map(data,data_to_stats)) %>%
    select(-data) %>% 
    unnest(stats)	


# we have time values for both time_courses (A/B: 
gs_stats %>% spread(time_course,stat_value) 

# we just merge the time_course A and timecourse B statsitcs with the participants. 
# note that many conditions disappear from the combined_statistics (left_join!),
# where there is not both A B available
mixed_stats <- gs_stats %>% spread(time_course,stat_value) %>% 
    left_join(gs_stats_original,by = c("cell_line", "treatment", "time", "stat_variable")) %>%
    left_join(combined_statistics, by = c("cell_line", "treatment", "time", "stat_variable"))

# all(tmp$standard ==tmp$standard_test)
selected_predictions <- names(mixed_stats)[c(5:7,10:25)]

condition_stats_SumSquared <- mixed_stats %>% 
    mutate(cond_id= paste(cell_line,treatment,time,sep = "_")) %>%
    group_by(cond_id) %>% 
    summarise_at(selected_predictions,~ sum((standard - .)^2))

condition_stats_SumSquared <- condition_stats_SumSquared %>% select(-standard_test)

# basic stats: 
condition_stats_SumSquared %>% 
    summarise_at(2:19,~ sqrt(sum(.))) %>%
    gather(team,score)

# show with points
condition_stats_SumSquared %>% 
    summarise_at(2:19,~ sqrt(sum(.))) %>%
    gather(team,score) %>%
    mutate(team = gsub("X.","",team,fixed = T)) %>%
    mutate(source = ifelse(team %in% c("A","B"),"Biological replicas","Team predictions")) %>%
    ggplot() + 
    geom_point(aes(source,score)) + coord_cartesian(ylim = c(0,70))+
    ggrepel::geom_text_repel(aes(source,score,label=team)) +
    scale_y_continuous(expand = c(0, 0)) + 
    xlab("") + ylab("Score (Mean, Cov)")

ggsave(filename = "./publication/figures/figure3/score_team_vs_bio_replicas.pdf",width = 3,height = 4)

condition_stats_SumSquared %>% 
    summarise_at(2:19,~ sqrt(sum(.)/n())) %>%
    gather(team,score) %>%
    mutate(team = gsub("X.","",team,fixed = T)) %>%
    mutate(source = ifelse(team %in% c("A","B"),"Biological replicas","Team predictions")) %>%
    ggplot() + 
    geom_point(aes(source,score)) + 
    #coord_cartesian(ylim = c(0,70))+
    ggrepel::geom_text_repel(aes(source,score,label=team)) +
    scale_y_continuous(expand = c(0, 0)) + 
    xlab("") + ylab("Score (Mean, Cov)")


hide_names <- tibble(submitterId = ranked_teams, alt_name = paste0("#",as.character(1:length(ranked_teams))))



condition_stats_SumSquared %>%  
    summarise_at(2:19,~ sqrt(sum(.))) %>%
    gather(submitterId,score) %>%
    left_join(hide_names, by="submitterId") %>%
    mutate(alt_name = ifelse(submitterId %in% c("A","B"),submitterId,alt_name)) %>%
    mutate(source = ifelse(submitterId %in% c("A","B"),"Biological replicas","Team predictions")) %>%
    ggplot() + 
    geom_point(aes(source,score)) + coord_cartesian(ylim = c(0,70))+
    ggrepel::geom_text_repel(aes(source,score,label=alt_name)) +
    scale_y_continuous(expand = c(0, 0)) + 
    xlab("") + ylab("Score (Mean, Cov)")

ggsave(filename = "./publication/figures/figure3/score_team_vs_bio_replicas_anonym.pdf",width = 3,height = 4)


### Bootstrapping over these conditions:
N_bootstrap = 1000

bootstrap_stats <- tibble(BS_sample = seq(N_bootstrap)) %>%
    mutate( BS_score = map(BS_sample, .f = ~ condition_stats_SumSquared %>% 
                    sample_frac(size = 1, replace = TRUE) %>%
                    summarise_at(2:19,~ sqrt(sum(.))))
    ) %>% unnest(BS_score)

hide_names_bio = hide_names %>%
    add_row(submitterId = "A",alt_name = "A",.before = 1)%>%
    add_row(submitterId = "B",alt_name = "B",.after = 1)

bootstrap_stats %>% gather(source,score,-BS_sample)%>%

    rename(submitterId ="source")%>%
    left_join(hide_names_bio, by="submitterId") %>%
    mutate(alt_name = factor(alt_name,levels = hide_names_bio$alt_name)) %>%
    ggplot() + geom_boxplot(aes(alt_name,score))  +
    xlab("Teams") + 
    ylab("Bootstrap Score (Median, Cov)") + 
    theme_bw()+ 
    theme(axis.text.x = element_text(hjust=1,angle = 45)
          #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #axis.ticks = element_blank(),
          #panel.border = element_blank(),
          #legend.position = "none"
          ) 


if(FALSE) write_rds(bootstrap_stats,"./submission_analysis/intermediate_data/sc2_bootstrap_stats_bioreplica.rds")



### Let's compute the Root mean square error, so we could compare SC2 and SC3:

# here we just compute the RMSE across all statistical value
mixed_stats %>% 
    ungroup() %>%
    mutate_at(selected_predictions,~ (standard - .)^2) %>%
    select(-standard_test,-standard) %>% 
    summarise_at(c(5:6,8:21),~ sqrt(sum((.)^2,na.rm = TRUE)/sum(!is.na(.)))) %>%
    gather(submitterId,score) %>%
    left_join(hide_names, by="submitterId") %>%
    mutate(alt_name = ifelse(submitterId %in% c("A","B"),submitterId,alt_name)) %>%
    mutate(source = ifelse(submitterId %in% c("A","B"),"Biological replicas","Team predictions")) %>%
    ggplot() + 
    geom_point(aes(source,score)) + 
    coord_cartesian(ylim = c(0,0.2))+
    ggrepel::geom_text_repel(aes(source,score,label=alt_name)) +
    scale_y_continuous(expand = c(0, 0)) + 
    theme_bw() +
    xlab("") + ylab("RMSE error") +
    theme(panel.grid = element_blank())

ggsave(filename = "./publication/figures/figure3/global_RMSE_team_vs_bio_replicas_anonym.pdf",width = 3,height = 4)

# here we compute the RMSE across all statistical value for each cell-line/time then average across condition
mixed_stats %>% 
    mutate(cond_id= paste(cell_line,treatment,time,sep = "_")) %>%
    group_by(cond_id) %>% 
    summarise_at(selected_predictions,~ sqrt(sum((standard - .)^2,na.rm = TRUE)/sum(!is.na(.)))) %>%
    ungroup() %>% 
    summarise_at(selected_predictions,mean) %>%
    select(-standard_test) %>% 
    gather(submitterId,score) %>%
    left_join(hide_names, by="submitterId") %>%
    mutate(alt_name = ifelse(submitterId %in% c("A","B"),submitterId,alt_name)) %>%
    mutate(source = ifelse(submitterId %in% c("A","B"),"Biological replicas","Team predictions")) %>%
    ggplot() + 
    geom_point(aes(source,score)) + 
    coord_cartesian(ylim = c(0,0.4))+
    ggrepel::geom_text_repel(aes(source,score,label=alt_name)) +
    scale_y_continuous(expand = c(0, 0)) + 
    theme_bw() +
    xlab("") + ylab("mean RMSE error") +
    theme(panel.grid = element_blank())
ggsave(filename = "./publication/figures/figure3/mean_RMSE_team_vs_bio_replicas_anonym.pdf",width = 3,height = 4)


mixed_stats %>% 
    mutate(cond_id= paste(cell_line,treatment,time,sep = "_")) %>%
    group_by(cond_id) %>% 
    summarise_at(selected_predictions,~ sqrt(sum((standard - .)^2,na.rm = TRUE)/sum(!is.na(.)))) %>%
    ungroup() %>% 
    summarise_at(selected_predictions,mean) 

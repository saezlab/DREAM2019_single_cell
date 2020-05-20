
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

fileID_table <- read_csv("./challenge_data/FileID_table.csv") %>% 
    select(cell_line,treatment,time,time_course,fileID)

ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc2_ranked_teams.rds")
combined_statistics <- read_rds("./submission_analysis/intermediate_data/sc2_stats_conditions.rds")




score_by_stat_type <- combined_statistics %>% 
    mutate(stat_type = ifelse(grepl("cov_",stat_variable),"covariance","mean")) %>%
    group_by(cond_id, cell_line,treatment,time,stat_type) %>%
    summarise_at(as.character(ranked_teams), ~ sum((.-standard)^2)) %>%
    ungroup() %>%
    group_by(stat_type) %>%
    summarise_at(as.character(ranked_teams), ~ sqrt(sum(.))) %>%
    group_by(stat_type) %>%
    gather(team, score,as.character(ranked_teams)) %>% 
    spread(stat_type,score)


# testing: by the pythagorian law, the two scores should result in the final score    
tmp <- score_by_stat_type %>% 
    mutate(f_score = sqrt(covariance^2 + mean^2)) %>%
    arrange(f_score) %>%
    pull(f_score)
tmp - SC_leaderboard$score


equi_score_lines <- expand.grid(r = seq(10,200,20), alpha= seq(0,3.14159/2,length.out = 100)) %>%
    mutate(x = r*sin(alpha),
           y=r*cos(alpha)) %>% 
    select(-alpha)

    
hide_names <- tibble(team = ranked_teams, alt_name = paste0("#",as.character(1:length(ranked_teams))))

score_by_stat_type %>% 
    left_join(hide_names, by="team") %>%
    ggplot(aes(covariance,mean)) + geom_point() +
    geom_line(data = equi_score_lines,  aes(x,y,group=r),size=0.1) + theme_bw()  + 
    coord_equal(xlim=c(0,150) , ylim=c(0,150),expand = FALSE) + 
    xlab("Covariance score") + ylab("Mean score") + theme(panel.grid = element_blank())+
    geom_abline(slope = 1,intercept = 0,color="grey50",size=0.1)+
    geom_text_repel(aes(label=alt_name)) 

ggsave(filename = "./publication/figures/figure3/score_mean_vs_covar.pdf",width = 4,height = 4)


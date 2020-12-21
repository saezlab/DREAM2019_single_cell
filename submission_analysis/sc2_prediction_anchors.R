#1. we compare the stats computed on timecourse A and time course B to 
# the the participants statistics.  
# 2. compute the score of submitting the EGF instead of iInhib+EGF



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


gs_stats_original <- gs %>% select(-fileID,-cellID) %>% data_to_stats(.) %>%
	arrange(cell_line, treatment,  time, stat_variable) %>% rename(standard_test=stat_value)


# find conditions, where there is both A and B condition measured
gs %>%	left_join(fileID_table, by = c("cell_line", "treatment", "time", "fileID")) %>%
	group_by(cell_line, treatment, time) %>% summarise(n_replica = length(unique(time_course))) %>% 
	filter(n_replica>1)

# these happens only at time 0. 
# for each time_course we compute the stats using the data_to_stats function, 
# this function internally considers the grouping by treatment, time and cell-line
gs_stats <- gs %>% filter(time==0) %>%
	left_join(fileID_table, by = c("cell_line", "treatment", "time", "fileID")) %>%
	select(-fileID,-cellID) %>%
	group_by(time_course) %>% nest() %>%
	mutate(stats = map(data,data_to_stats)) %>%
	unnest(stats)	
	

# we have time values for both timecourses: 
gs_stats %>% spread(time_course,stat_value) 

# we just merge the time_course A and timecourse B statsitcs with the participants. 
# note that many conditions disappear from the combined_statistics (left_join!),
# where there is not both A B available
mixed_stats <- gs_stats %>% spread(time_course,stat_value) %>% 
	left_join(gs_stats_original,by = c("cell_line", "treatment", "time", "stat_variable")) %>%
	left_join(combined_statistics, by = c("cell_line", "treatment", "time", "stat_variable"))
	
# all(tmp$standard ==tmp$standard_test)
selected_predictions <- names(mixed_stats)[c(5:7,10:25)]

condition_stats_SumSquared <- mixed_stats %>% mutate(cond_id= paste(cell_line,treatment,time,sep = "_")) %>%
	group_by(cond_id) %>% 
	summarise_at(selected_predictions,~ sum((standard - .)^2))


# Compute the bootstraps on the subset of data: 
N_bootstrap = 1000

bootstrap_stats <- tibble(BS_sample = seq(N_bootstrap)) %>%
	mutate( map(BS_sample, .f = ~ condition_stats_SumSquared %>% 
					sample_frac(size = 1, replace = TRUE) %>%
					summarise_at(as.character(selected_predictions),~ sqrt(sum(.))))
	) %>% unnest()

bootstrap_stats %>% select(-standard_test) %>% 
	gather(teams,score,-BS_sample) %>%
	mutate(teams = factor(teams,levels = c("A","B",levels(ranked_teams)))) %>%
	ggplot() + geom_boxplot(aes(teams, score,color=teams)) + 
	ggtitle("Subchallenge II: performance over bootstrap samples") + 
	theme_bw() +
	theme(axis.text.x = element_text(hjust=1,angle = 45))



# Show the 2 types of scores separately : covariance and median

med_cov_stats <- mixed_stats %>% 
	mutate(cond_id= paste(cell_line,treatment,time,sep = "_")) %>%
	mutate(stat_type= ifelse(grepl("cov_",stat_variable),"cov","mean")) %>%
	group_by(stat_type) %>% 
	summarise_at(selected_predictions,~ sum((standard - .)^2))

med_cov_stats %>% select(-standard_test) %>%
	gather(teams,score,-stat_type) %>% spread(stat_type,score) %>%
	ggplot(aes(x=cov,y=mean)) + geom_point() +
	ggrepel::geom_text_repel(aes(label=teams)) +
	#coord_cartesian(xlim = c(9,3000),ylim = c(9,3000)) +
	scale_x_log10()+
	scale_y_log10() + 
	coord_equal(xlim = c(9,3000),ylim = c(9,3000)) +
	theme_bw() + xlab("score based on covariance") +ylab("score based on mean")


# show with barplot


condition_stats_SumSquared %>% select(-standard_test) %>% 
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


hide_names <- tibble(submitterId = ranked_teams, alt_name = paste0("#",as.character(1:length(ranked_teams))))



condition_stats_SumSquared %>% select(-standard_test) %>% 
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




condition_stats_SumSquared %>% summarise_at(selected_predictions,~ sqrt(sum((.)^2)))


### Idea 2. --------------------------------------------------------------------
# from Atta participant: what about using EGF as an estimate for the other kinase inhibition


library(tidyverse)
# SC2

combined_statistics <- read_rds("./submission_analysis/intermediate_data/sc2_stats_conditions.rds")
source("./scoring_scripts/score_sc2.R")
# load the single cell data
EGF_data <- list.files("./challenge_data/single_cell_phospho/subchallenge_2/", full.names = T) %>%
    lapply(., function(sc_file){
        sc_data <- read_csv(sc_file)
        sc_data %>% filter(treatment=="EGF")
    }) %>% bind_rows()

EGF_stats <- EGF_data %>%
    select(-cellID,-fileID) %>% 
    data_to_stats()	
    
EGF_stats <- EGF_stats %>% 
    rename("pred_treatment"="treatment") %>%
    rename("EGF_stat"="stat_value")


# we have to fix the timing:
# - EGF and inhibitors have different times
# - luckily only {53 MDAMB231  iPI3K      16   NA } is missing from the EGF, 
# because it is {107 MDAMB231  NA         17   EGF}
# run this to check:  
# full_join(combined_statistics %>% select(cell_line,treatment,time) %>% unique(),
# EGF_stats %>% select(cell_line,pred_treatment,time) %>% unique()) %>% print(n=121)
#
# to fix, we overwrite the EGF 17 to EGF 16. (that is the clsest timepoint)
EGF_stats <- EGF_stats %>% mutate(time = ifelse(cell_line == "MDAMB231" & time==17, 16, time))

selected_predictions <- names(combined_statistics)[c(7:22)]

all_stats <- left_join(combined_statistics, EGF_stats, by = c("cell_line", "time", "stat_variable")) %>%
    select(1:5, EGF_stat, everything()) %>%
    summarise_at(c("EGF_stat",selected_predictions),~ sqrt(sum((standard - .)^2)))
# %>%
#     filter(is.na(EGF_stat)) %>% select(cell_line, treatment, time) %>% unique()

all_stats$EGF_stat

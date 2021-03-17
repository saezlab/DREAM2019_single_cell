# In SC2 and SC3 we have bio replica on the subsets. 
# plot them on the same figure. 
library(tidyverse)
library(ggrepel)
library(ggplot2)


stats = read_rds("./submission_analysis/intermediate_data/sc2_stats_conditions.rds")
sc2_bootstrap_stats = read_rds("./submission_analysis/intermediate_data/sc2_bootstrap_stats_bioreplica.rds")
sc3_bootstrap_stats = read_rds("./submission_analysis/intermediate_data/sc3_bootstrap_stats_bioreplica.rds")

sc2_stats_long <- sc2_bootstrap_stats %>% gather(id,score,-BS_sample) %>% mutate(sc = 2)
sc3_stats_long <- sc3_bootstrap_stats %>% gather(id,score,-BS_sample) %>% mutate(sc = 3)


sc2_ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc2_ranked_teams.rds") %>% as.character()
sc3_ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc3_ranked_teams.rds") %>% as.character()

# combine the participants of the 2 subchallange
# notice how many they participated and what is their average rank
teams <- bind_rows(
    tibble(team = sc2_ranked_teams, rank = 1:length(sc2_ranked_teams),sc = 2),
    tibble(team = sc3_ranked_teams, rank = 1:length(sc3_ranked_teams),sc = 3 ))

team_ranks <- teams %>%
    group_by(team) %>%
    summarise(N_challenge = n(),
              mean_rank = mean(rank)) %>% arrange(mean_rank)




# plot the performance of the teams and biological replica with boxplots
# note that the number of conditions are different in SC2 and SC3

bind_rows(sc2_stats_long,sc3_stats_long) %>%
    mutate(id = factor(id,levels = (c("A","B",team_ranks$team)))) %>%
    ggplot() + geom_boxplot(aes(id,score,fill=factor(sc)))  +
    xlab("Teams") + 
    ylab("Bootstrap Score (Median, Cov)") + 
    theme_bw()+ 
    theme(axis.text.x = element_text(hjust=1,angle = 45)
          #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #axis.ticks = element_blank(),
          #panel.border = element_blank(),
          #legend.position = "none"
    ) 


# how does it look if we compensate for the number of conditions
# SC2: 12
# SC3: 38
bind_rows(sc2_stats_long %>% mutate(scaled_score = score/sqrt(12)),
          sc3_stats_long %>%mutate(scaled_score = score/sqrt(38))) %>%
    mutate(id = factor(id,levels = (c("A","B",team_ranks$team)))) %>%
    ggplot() + geom_boxplot(aes(id,scaled_score,fill=factor(sc)))  +
    xlab("Teams") + 
    ylab("Bootstrap Score (Median, Cov)") + 
    theme_bw()+ 
    theme(axis.text.x = element_text(hjust=1,angle = 45)
          #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #axis.ticks = element_blank(),
          #panel.border = element_blank(),
          #legend.position = "none"
    ) 



team_ranks <- team_ranks %>% mutate(alt_name = as.character(1:n())) %>% 
    add_row(team = "A",N_challenge=2,mean_rank = 0,alt_name="replica A") %>%
    add_row(team = "B",N_challenge=2,mean_rank = 0.5,alt_name="replica B") %>% 
    arrange(mean_rank) %>%
    mutate(alt_name = factor(alt_name,levels = alt_name))

# how does it look if we compensate for the different number of conditions (Cell_line and stats)
# SC2: 12 cell_line * 665 (medians and covariance pairs) = 7980
# SC3: 38 * 665 = 25270
bind_rows(sc2_stats_long %>% mutate(scaled_score = score/sqrt(7980)),
          sc3_stats_long %>%mutate(scaled_score = score/sqrt(25270))) %>%
    left_join(team_ranks,by=c("id"="team")) %>%
    rename(Subchallenge = "sc") %>%
    mutate(Subchallenge = factor(paste0("SC ",Subchallenge))) %>%
    ggplot() + 
    geom_boxplot(aes(alt_name,scaled_score,fill=Subchallenge),outlier.size = .5,size=0.1)  +
    xlab("Teams") + 
    ylab("Error per estimated statistics") + 
    theme_bw()+ 
    theme(axis.text.x = element_text(hjust=1,angle = 45),legend.position = c(0.2,0.8)
    ) 

ggsave("./publication/figures/sc2_sc3_error_per_stats.pdf",width = 6,height = 4)





bind_rows(sc2_stats_long %>% mutate(scaled_score = score/sqrt(7980)),
          sc3_stats_long %>%mutate(scaled_score = score/sqrt(25270))) %>% 
    filter(id %in% c("A","B","icx_bxai")) %>% group_by(id,sc) %>%
    summarise(mean_score = mean(scaled_score),
              sd_score = sd(scaled_score)) %>%
    arrange(sc)



sd((stats$standard))






bind_rows(sc2_stats_long %>% mutate(scaled_score = score/sqrt(12)),
          sc3_stats_long %>%mutate(scaled_score = score/sqrt(38))) %>%
    ggplot()  + geom_boxplot(aes(factor(sc),scaled_score))

bind_rows(sc2_stats_long %>% mutate(scaled_score = score),
          sc3_stats_long %>%mutate(scaled_score = score)) %>%
    ggplot()  + geom_boxplot(aes(factor(sc),scaled_score))

bind_rows(sc2_stats_long %>% mutate(scaled_score = score/sqrt(7980)),
          sc3_stats_long %>%mutate(scaled_score = score/sqrt(25270))) %>%
    ggplot()  + geom_boxplot(aes(factor(sc),scaled_score))

# scale the scores by the scaore of the respective biological replica

all_scores <- bind_rows(sc2_stats_long,sc3_stats_long) 

mean_bio_score <- all_scores %>% 
    filter(id %in% c("A","B")) %>% 
    group_by(sc) %>% 
    summarise(mean_bio = mean(score))

all_scores %>% left_join(mean_bio_score) %>%
    
    mutate(scaled_score = score/mean_bio-1) %>%
    mutate(id = factor(id,levels = (c("A","B",team_ranks$team)))) %>%
    ggplot() + geom_boxplot(aes(id,scaled_score,fill=factor(sc)))  +
    xlab("Teams") + 
    ylab("Bootstrap Score (Median, Cov)") + 
    theme_bw()+ 
    theme(axis.text.x = element_text(hjust=1,angle = 45)
          #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #axis.ticks = element_blank(),
          #panel.border = element_blank(),
          #legend.position = "none"
    ) 


# Plot with density SC2
dens_plt <- all_scores %>% 
    mutate(id = factor(id,levels = (c("A","B",team_ranks$team)))) %>%
    filter(sc==2)

teams_score <- filter(dens_plt,! id %in% c("A","B")) %>% group_by(id) %>% summarise(ms = mean(score))

ggplot() + 
    geom_density(data = filter(dens_plt,id %in% c("A","B")), aes(score,fill=1),alpha=.3)  +
    geom_point(data = teams_score , aes(ms,0)) +
    geom_text_repel(data = teams_score , aes(ms,0,label=id),ylim=c(0,1)) +
    xlab("Bootstrap Score (Median, Cov)") + 
    theme_bw()

# Plot with density SC3
dens_plt <- all_scores %>% 
    mutate(id = factor(id,levels = (c("A","B",team_ranks$team)))) %>%
    filter(sc==3)

teams_score <- filter(dens_plt,! id %in% c("A","B")) %>% group_by(id) %>% summarise(ms = mean(score))

ggplot() + 
    geom_density(data = filter(dens_plt,id %in% c("A","B")), aes(score),fill="pink",alpha=.3)  +
    geom_point(data = teams_score , aes(ms,0)) +
    geom_text_repel(data = teams_score , aes(ms,0,label=id),ylim=c(0,1)) +
    xlab("Bootstrap Score (Median, Cov)") + 
    theme_bw()


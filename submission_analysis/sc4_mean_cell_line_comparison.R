## 


library(tidyverse)
library(ggpubr)

median_predictions <- read_rds("./submission_analysis/intermediate_data/sc4_combined_data.rds")

median_training <- read_csv("challenge_data/median_phospho/median_phospho_data.csv")

ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc4_ranked_teams.rds")

marker_columns <- unique(median_predictions$marker)



average_cell_line <- median_training %>%
    gather(marker,value,-cell_line, -treatment, -time) %>%
    group_by(treatment,time,marker) %>%
    summarise(ave_cell_line = median(value,na.rm=TRUE),
              sd_cell_line = sd(value,na.rm=TRUE)) 

## Calculate statistics:

RMSE <- median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    group_by(cell_line, treatment, marker) %>%
    summarise_at(c("ave_cell_line",as.character(ranked_teams)),~ sqrt(sum((standard - .)^2) / n()))


RMSE_all <- median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    select(-sd_cell_line) %>%
    gather(team, prediction, -cell_line, -treatment, -marker, -time,-standard) %>% 
    group_by(cell_line,treatment,marker,team) %>%
    summarise(RMSE = sqrt(sum((standard - prediction)^2) / n()),
              trajectory_sd = sd(standard),
              m_measurement = mean(standard))


ranked_teams = RMSE_all %>% group_by(team) %>% summarise(mean_RMSE = mean(RMSE)) %>% arrange(mean_RMSE) %>% pull(team)

#Range of data
range(median_predictions$standard)
median(median_predictions$standard)
hist(median_predictions$standard)


# RMSE across treatments: uniform
RMSE %>% select(cell_line, treatment, marker, ave_cell_line) %>%
    group_by(treatment) %>% summarise(med_rmse = mean(ave_cell_line))

RMSE %>% select(cell_line, treatment, marker, ave_cell_line) %>%
    ggplot() + geom_violin(aes(treatment,ave_cell_line),draw_quantiles = .5)



# RMSE across cell-lines and nodes: 

RMSE %>% select(cell_line, treatment, marker, ave_cell_line) %>%
    ggplot() + geom_boxplot(aes(marker, ave_cell_line)) + facet_grid(cell_line~.) +
    theme(axis.text.x = element_text(angle = 90,hjust = 1))


# Badly approximated celllin/node: RMSE is larger than 0.5
RMSE %>% ungroup() %>% select(cell_line, marker, ave_cell_line) %>%
    group_by(marker) %>%
    summarise(mean_rmse = mean(ave_cell_line)) %>%
    arrange(desc(mean_rmse))
    
# well approximated celllin/node: RMSE is larger than 0.5
RMSE %>% ungroup() %>% select(cell_line, marker, ave_cell_line) %>%
    group_by(marker) %>%
    summarise(mean_rmse = mean(ave_cell_line)) %>%
    arrange(mean_rmse)

## Heatmap

RMSE %>% ungroup() %>% select(cell_line, marker, ave_cell_line) %>%
    group_by(marker,cell_line) %>%
    summarise(mean_rmse = mean(ave_cell_line)) %>%
    spread(marker,mean_rmse)    %>% column_to_rownames("cell_line") %>% 
    pheatmap::pheatmap(color = RColorBrewer::brewer.pal(7,"Reds"))
    
    
## IS RMSE correlating with the variance of the data? no
RMSE_all %>% filter(team == "ave_cell_line") %>%
    mutate(lab = ifelse( (RMSE<0.5 & trajectory_sd>0.9) | (RMSE>2 & trajectory_sd<0.3) ,paste(cell_line, treatment, marker),"")) %>%
    ggplot(aes(RMSE,trajectory_sd)) + geom_point() 

RMSE_all %>% filter(team %in%  c("ave_cell_line",ranked_teams[1:3])) %>%
    mutate(lab = ifelse( (RMSE<0.5 & trajectory_sd>0.9) | (RMSE>2 & trajectory_sd<0.3) ,paste(cell_line, treatment, marker),"")) %>%
    ggplot(aes(RMSE,trajectory_sd,color = team)) + geom_point() + facet_wrap(~team)



RMSE_all %>% filter(team == "ave_cell_line") %>%
    ggplot(aes(RMSE,m_measurement)) + geom_point()


## The RMSE strongly correlates with the variance of nodes between cell-line. 
RMSE_all %>% filter(team == "ave_cell_line") %>%
    group_by(treatment,marker) %>%
    summarise(mean_RMSE = mean(RMSE),
              max_RMSE = max(RMSE),
              sd_over_cell_line = sd(m_measurement)) %>%
    ggplot() + geom_point(aes(mean_RMSE,sd_over_cell_line))
    


RMSE_all %>% filter(team %in% c("ave_cell_line",ranked_teams[1:3])) %>%
    group_by(treatment,marker, team) %>%
    summarise(mean_RMSE = mean(RMSE),
              max_RMSE = max(RMSE),
              sd_over_cell_line = sd(m_measurement)) %>%
    ggplot() + geom_point(aes(mean_RMSE,sd_over_cell_line,col=marker)) + facet_wrap(~team)





### When the predictions improved the most
selected_teams = c("ave_cell_line","icx_bxai","AMbeRland","orangeballs")

marker_order = RMSE_all %>% filter(team=="ave_cell_line") %>% 
    group_by(marker) %>%
    summarise(mean_rmse = mean(RMSE)) %>%
             arrange(desc(mean_rmse)) %>% pull(marker)

    
RMSE_all %>% ungroup() %>% select(cell_line, marker, team,RMSE) %>%
    #filter(team %in% selected_teams) %>%
    group_by(marker,team) %>%
    summarise(mean_rmse = mean(RMSE),.groups="drop") %>%
    arrange(desc(mean_rmse)) %>% 
    mutate(marker = factor(marker,levels = marker_order)) %>%
    mutate(team = factor(team,levels = ranked_teams)) %>%
    mutate(data = ifelse(team %in% ranked_teams[c(1:3)],"top 3 teams","other teams"),
           data = ifelse(team == "ave_cell_line","average cell line",data)) %>%
    ggplot(aes(marker,mean_rmse,col=data)) + 
    geom_point(aes(alpha = data)) + 
    geom_line(aes(group = team,alpha = data)) +
    theme_bw()+ 
    theme(axis.text = element_text(angle = 90,hjust=1)) + 
    scale_alpha_manual(values = c(`top 3 teams`=0.5,`average cell line`=1,`other teams`=0.15))


# exmaple where the predictions are much better than the average cell-line: p53
median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    select(-sd_cell_line) %>%
    gather(team, prediction, -cell_line, -treatment, -marker, -time,-standard) %>%
    filter(marker == "p.p53") %>%
    filter(team %in% ranked_teams[c(1:3,14)]) %>%
    ggplot() + 
    geom_line(aes(time,prediction,group=team,col=team)) +
    geom_line(aes(time,standard),col="red",size=1) +
    facet_grid(cell_line~treatment)


### Do the predictions following the dynamics? 

# lets take the top 3 most varying marker in EGF
traj_sd <- RMSE_all %>% filter(treatment=="EGF") %>%
    arrange(desc(trajectory_sd)) %>% select(cell_line,treatment,marker,trajectory_sd) %>%
    pull(marker) %>% unique() %>% head(5)


median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    select(-sd_cell_line) %>%
    gather(team, prediction, -cell_line, -treatment, -marker, -time,-standard) %>%
    filter(team %in% ranked_teams[c(1:3,14)]) %>%
   filter(treatment %in% "EGF") %>% 
    filter(marker %in% traj_sd) %>%
    ggplot() + 
    geom_line(aes(time,prediction,group=team,col=team)) +
    geom_line(aes(time,standard),col="red",size=1) +
    facet_grid(marker~cell_line,scales = "free")




###  MeK indepenendent ERK signaling? Not really present
median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    select(-sd_cell_line) %>%
    gather(team, prediction, -cell_line, -treatment, -marker, -time,-standard) %>%
    filter(marker == "p.ERK") %>%
    filter(treatment %in% c("EGF","iMEK")) %>%
    filter(team %in% ranked_teams[c(1:3,14)]) %>%
    ggplot() + 
    geom_line(aes(time,prediction,group=team,col=team)) +
    geom_line(aes(time,standard),col="red",size=1) +
    facet_grid(cell_line~treatment,scales = "free") 






# variance of the marker across cell-line:
marker_variance_cell_line <- RMSE_all %>% group_by(marker,treatment) %>% summarise(sd_cell_line = sd(m_measurement)) %>%
    group_by(marker) %>% summarise(msd_cell_line = max(sd_cell_line))

# variance of marker across time: 
marker_variance_time <- RMSE_all %>% group_by(marker) %>% 
    summarise(timesd = max(trajectory_sd)) 



### RMSE, time and cell-line variance 

marker_order = RMSE_all %>% filter(team=="ave_cell_line") %>% 
    group_by(marker) %>%
    summarise(mean_rmse = mean(RMSE)) %>%
    arrange(desc(mean_rmse)) %>% pull(marker)


rmse_plot <- RMSE_all %>% ungroup() %>% select(cell_line, marker, team,RMSE) %>%
    #filter(team %in% selected_teams) %>%
    group_by(marker,team) %>%
    summarise(mean_rmse = mean(RMSE),.groups="drop") %>%
    mutate(marker = factor(marker,levels = marker_order)) %>%
    mutate(team = factor(team,levels = ranked_teams)) %>%
    mutate(data = ifelse(team %in% ranked_teams[c(1:3)],"top 3 teams","other teams"),
           data = ifelse(team == "ave_cell_line","average cell line",data)) %>%
    ggplot(aes(marker,mean_rmse,col=data)) + 
    geom_point(aes(alpha = data)) + 
    geom_line(aes(group = team,alpha = data)) +
    theme_bw()+ 
    theme(axis.text.x = element_blank(),legend.position = c(0.8,0.7),legend.title = element_blank()) + 
    scale_alpha_manual(values = c(`top 3 teams`=0.5,`average cell line`=1,`other teams`=0.15)) +
    xlab("")+ylab("average RMSE")

sd_cell_line_plot <- marker_variance_cell_line %>%
    mutate(marker = factor(marker,levels = marker_order)) %>%
    ggplot() + geom_col(aes(marker,msd_cell_line)) + theme_bw()+ 
    theme(axis.text.x = element_blank())  +xlab("") +ylab("variance between \n cell lines")

sd_time_plot <- marker_variance_time %>%
    mutate(marker = factor(marker,levels = marker_order)) %>%
    ggplot() + geom_col(aes(marker,timesd)) + theme_bw()+ 
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))  + ylab("variance in time")

cowplot::plot_grid(rmse_plot,sd_cell_line_plot,sd_time_plot, axis = "lr", align = "v", rel_heights =  c(0.5,0.3,0.5), ncol = 1)

ggsave("./publication/figures/figure5/RMSE_comparison_ave_cellline.pdf",width = 6,height = 6)




RMSE_cd_correlation <- RMSE_all %>% group_by(marker,team) %>%summarise(RMSE = mean(RMSE) ) %>%
    left_join(marker_variance_cell_line) %>% left_join(marker_variance_time) %>%
    filter(team %in% c("ave_cell_line",ranked_teams[1:3])) %>%
    group_by(team) %>%
    summarise(cor_cell_line_sd = cor(RMSE,msd_cell_line),
              cor_time_sd = cor(RMSE,timesd))

print(RMSE_cd_correlation)




marker_variance_cell_line %>% left_join(marker_variance_time) %>%
    ggplot() + geom_col(aes(marker,msd_cell_line))+
    geom_col(aes(marker,-timesd)) + geom_hline(yintercept = 0,size=2,col="white")







# improvement comapred to average: 

RMSE_improvement <- RMSE_all %>%  filter(team %in% ranked_teams[c(1:3,14)])  %>%
    select(cell_line,treatment,marker, team,RMSE) %>%
    spread(team,RMSE) %>%
    mutate(mean_prediction_error = pmap_dbl(list(AMbeRland,icx_bxai,orangeballs),function(a,b,c){
        mean(c(a,b,c))
    })) %>% group_by(marker) %>%
    summarise(RMSEave_cell_line = mean(ave_cell_line),
              RMSEpredictions = mean(mean_prediction_error)) %>%
    mutate(RMSE_improvement = RMSEave_cell_line - RMSEpredictions )


RMSE_improvement %>% left_join(marker_variance_time) %>%
    summarise(corr_improv_sdtime = cor(timesd,RMSE_improvement))



RMSE_improvement %>% left_join(marker_variance_cell_line) %>%
    left_join(marker_variance_time) %>%
    ggplot() + geom_point(aes(msd_cell_line,timesd,col=RMSE_improvement)) +
    scale_color_distiller(palette = "Spectral", limits =c(-0.5,0.5) ) + theme_bw()



RMSE_improvement %>% left_join(marker_variance_cell_line) %>%
    left_join(marker_variance_time) %>%
    ggplot() + geom_point(aes(msd_cell_line,timesd,col=RMSE_improvement)) +
    geom_label(aes(msd_cell_line,timesd,label = marker, col=RMSE_improvement)) +
    scale_color_distiller(palette = "YlOrRd",direction = 1 ) + theme_bw()

RMSE_improvement %>% left_join(marker_variance_cell_line) %>%
    left_join(marker_variance_time) %>%
    ggplot() + geom_point(aes(msd_cell_line,timesd,col=RMSEpredictions)) +
    geom_label(aes(msd_cell_line,timesd,label = marker, col=RMSEpredictions)) +
    scale_color_distiller(palette = "Spectral", oob = scales::squish ) + theme_bw()

RMSE_improvement %>% left_join(marker_variance_cell_line) %>%
    left_join(marker_variance_time) %>%
    ggplot() + geom_point(aes(msd_cell_line,timesd,col=RMSEave_cell_line)) +
    geom_label(aes(msd_cell_line,timesd,label = marker, col=RMSEave_cell_line)) +
    scale_color_distiller(palette = "Spectral", oob = scales::squish ) + theme_bw()

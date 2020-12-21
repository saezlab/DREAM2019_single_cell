# can we use variance explained instead of RMSE? 

library(tidyverse)
library(ranger)



# Participating team names
submissions <-  readRDS("./submission_analysis/intermediate_data/sc1_ranked_teams.rds") %>% as.character()


sub_data_all <- read_rds("./submission_analysis/intermediate_data/sc1_all_NP_predictions.rds")  %>%
    mutate_if(is.character, as.factor) 


error <- sub_data_all %>% select(1:7,"standard",submissions ) %>% 
    mutate_at(submissions,~ .-standard)


# showing that the error is correlated
#error %>% select(submissions) %>% cor ()
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc1_ranked_teams.rds")



# compute the variance exaplained for the teams per condition




conditional_corr <- sub_data_all %>% select(1:7,"standard",submissions ) %>%
    group_by(cell_line,treatment,time,marker) %>%
    summarise_at(submissions, ~cor(x = .,y=standard))


conditional_corr%>% ggplot() + geom_boxplot(aes(1,icx_bxai)) 

# format to table: 
# orangeballs submitted constant prediction --> correlation is NA
conditional_corr %>%  gather(submitterId,correlation,-cell_line, -treatment,  -time, -marker) %>%
    filter(!is.na(correlation)) %>% 
    group_by(submitterId) %>% summarise(min = min(correlation),
                                        `1st Qu` = quantile(correlation,0.25),
                                        mean = mean(correlation),
                                        median  = median(correlation),
                                        `3rd Qu` = quantile(correlation,0.75),
                                        max = max(correlation)) %>%
    print(.,n=22)



conditional_R2 <- error %>% group_by(cell_line,treatment,time,marker) %>%
    # center the measurements: 
    mutate(standard = standard - mean(standard)) %>%
    # compute the residual and total sum of squares:
    summarise_at(c("standard",submissions),~ sum(.^2)) %>%
    # compute the variance explained: 1- SS_res/SS_tot
    group_by(cell_line,treatment,time,marker) %>%
    mutate_at(submissions,~ 1 - ./standard)



conditional_R2 %>% ggplot() + 
    geom_boxplot(aes(treatment,icx_bxai)) + 
    facet_grid(cell_line~marker) + coord_cartesian(ylim = c(-1,1))+
    theme(axis.text.x = element_text(angle = 90)) + 
    ggtitle("R2 of best performing model across conditions") 

# format to table: 

conditional_R2 %>% select(-standard) %>% gather(submitterId,R2,-cell_line, -treatment,  -time, -marker) %>%
    group_by(submitterId) %>% summarise(min = min(R2),
                                        `1st Qu` = quantile(R2,0.25),
                                        mean = mean(R2),
                                        median  = median(R2),
                                        `3rd Qu` = quantile(R2,0.75),
                                        max = max(R2)) %>%
    print(.,n=22)




hide_names <- tibble(submitterId = c("data",as.character(ranked_teams)),
                     alt_name = c("data", as.character(1:length(ranked_teams))))

hide_names$alt_name = factor(hide_names$alt_name,levels =hide_names$alt_name)



full_join(conditional_R2 %>% select(-standard) %>% gather(submitterId,R2,-cell_line, -treatment,  -time, -marker),
          conditional_corr  %>% gather(submitterId,corr,-cell_line, -treatment,  -time, -marker),
          by = c("cell_line", "treatment", "time", "marker", "submitterId")) %>%
    left_join(hide_names,by = "submitterId") %>% 
    gather(stat,value,R2,corr) %>% 
    ggplot() + geom_boxplot(aes(alt_name,value,fill = alt_name, group=paste(alt_name,stat))) + 
    coord_cartesian(ylim=c(-1,1)) + facet_wrap(~stat) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90,hjust=1)) + guides(fill="none") +
    ylab("") + xlab("teams")


    
ggsave("./publication/figures/figure2/sc1_corr_R2.pdf", height = 3, width = 5)


conditional_corr  %>%
    gather(submitterId,corr,-cell_line, -treatment,  -time, -marker)%>%
    left_join(hide_names,by = "submitterId") %>% 
    gather(stat,value,corr) %>% 
    ggplot() +
    geom_boxplot(aes(alt_name,value, group=paste(alt_name,stat)),outlier.size = 0.2,) + 
    coord_cartesian(ylim=c(-1,1)) + 
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    guides(fill="none") +
    ylab("Correlation") + 
    xlab("Ranked teams")

ggsave("./publication/figures/figure2/sc1_correlation.pdf", height = 2, width = 4)

conditional_R2 %>% 
    select(-standard) %>%
    gather(submitterId, R2, -cell_line, -treatment,  -time, -marker) %>%
    left_join(hide_names,by = "submitterId") %>% 
    gather(stat,value,R2) %>% 
    ggplot() +
    geom_boxplot(aes(alt_name,value, group=paste(alt_name,stat)),outlier.size = 0.2,) + 
    coord_cartesian(ylim=c(-2,1)) + 
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    guides(fill="none") +
    ylab("Explained variance") + 
    xlab("Ranked teams")

ggsave("./publication/figures/figure2/sc1_R2.pdf", height = 2, width = 4)

conditional_R2 %>% ungroup() %>%
    select(submissions) %>% 
    pheatmap::pheatmap(mat=.,
                       color = RColorBrewer::brewer.pal(8,"RdYlBu"),
                       breaks = seq(-1,1,length.out = 7),
                       
                       show_rownames = FALSE)

# compute the variance explained for the teams 
error %>% ungroup() %>%
    # center the measurements: 
    mutate(standard = standard - mean(standard)) %>%
    # compute the residual and total sum of squares:
    summarise_at(c("standard",submissions),~ sum(.^2)) %>%
    # compute the variance explained: 1- SS_res/SS_tot
    #group_by(marker) %>%
    mutate_at(submissions,~ 1 - ./standard) %>% 
    ungroup() %>% gather(team, R2, submissions) %>%
    group_by(team) %>% 
    summarise(mean_r2 = mean(R2)) %>% 
    arrange(desc(mean_r2)) %>%
    mutate(team = factor(team,levels = team)) %>%
    ggplot() + geom_point(aes(team,mean_r2)) + ylim(-1,1)



error %>% ungroup() %>%
    # center the measurements: 
    mutate(standard = standard - mean(standard)) %>%
    group_by(cell_line,treatment,time,marker) %>%
    # compute the residual and total sum of squares:
    summarise_at(c("standard",submissions),~ sum(.^2)) %>%
    # compute the variance explained: 1- SS_res/SS_tot
    group_by(cell_line,treatment,time,marker) %>%
    mutate_at(submissions,~ 1 - ./standard) %>% 
    ungroup() %>% gather(team, R2, submissions) %>%
    # group_by(team) %>% 
    # summarise(mean_r2 = mean(R2)) %>% 
    #arrange(desc(mean_r2)) %>%
    mutate(team = factor(team,levels = submissions)) %>%
    ggplot() + geom_violin(aes(team,R2),draw_quantiles = c(0.5)) +
    coord_cartesian(ylim = c(-1, 1))




error %>% ungroup() %>%
    # center the measurements: 
    mutate(standard = standard - mean(standard)) %>%
    group_by(cell_line,treatment,time,marker) %>%
    # compute the residual and total sum of squares:
    summarise_at(c("standard",submissions),~ sqrt(sum(.^2)/n())) 



conditional_R2 %>% 
    select(-standard) %>%
    gather(submitterId, R2, -cell_line, -treatment,  -time, -marker) %>%
    left_join(hide_names,by = "submitterId") %>%
    
    ggplot() + geom_violin(aes(cell_line,R2),draw_quantiles = .5) + ylim(c(0,1))



conditional_R2 %>% 
    select(-standard) %>%
    gather(submitterId, R2, -cell_line, -treatment,  -time, -marker) %>%
    filter(R2>0.5) %>% pull(cell_line) %>% table()


conditional_R2 %>% 
    select(-standard) %>%
    gather(submitterId, R2, -cell_line, -treatment,  -time, -marker) %>%
    filter(R2>0.5) %>% pull(marker) %>% table()


conditional_R2 %>% 
    select(-standard) %>%
    gather(submitterId, R2, -cell_line, -treatment,  -time, -marker) %>%
    filter(R2>0.5) %>% pull(treatment) %>% table()




#### Plot RMSE example where the RMSE is around the scored values
error %>% ungroup() %>%
    # center the measurements: 
    mutate(standard = standard - mean(standard)) %>%
    group_by(cell_line,treatment,time,marker) %>%
    # compute the residual and total sum of squares:
    summarise_at(submissions,~ sqrt(sum(.^2)/n())) 

tmp <- .Last.value

tmp %>% ungroup() %>% select(submissions) %>% colMeans(.)

# select 3 representative conditions: 
# cond A: close to the mean score
tmp %>% ungroup() %>%
    filter(icx_bxai < mean(tmp$icx_bxai)*1.005, icx_bxai > mean(tmp$icx_bxai)*0.998)
# cond B: around .75 percentile (worse RMSE)
tmp %>% ungroup() %>%
    filter(icx_bxai == quantile(icx_bxai,0.75,type = 1)) %>% arrange(icx_bxai)
# cond C: around .25 percentile (good RMSE)
tmp %>% ungroup() %>% 
    filter(icx_bxai < quantile(icx_bxai,0.25)) %>% arrange(desc(icx_bxai))






bind_rows(sub_data_all %>% filter(cell_line == 'AU565',treatment=="EGF",time ==7,marker == "p.S6") %>%
              mutate(cond = "mean: 0.853"),
          sub_data_all %>% filter(cell_line == 'LY2',treatment=="iPI3K",time ==40,marker == "p.PLCg2")%>%
              mutate(cond = "0.25 percentile: 0.619"),
          sub_data_all %>% filter(cell_line == 'MACLS2',treatment=="iPI3K",time ==40,marker == "p.HER2")%>%
              mutate(cond = "0.75 percentile: 0.965")) %>% 
    select(cellID, standard,icx_bxai,cond) %>% 
    group_by(cond) %>% 
    sample_n(200) %>% 
    gather(source,value,-cellID,-cond) %>% 
    ungroup() %>%
    mutate(cond = factor(cond)) %>% 
    mutate(cond = factor(cond,levels = levels(cond)[c(1,3,2)])) %>%
    ggplot(aes(source,value)) +
    geom_violin() +
    geom_line(aes(group=cellID), color="grey80")  +
    geom_point() + facet_wrap(~cond)



bind_rows(sub_data_all %>% filter(cell_line == 'AU565',treatment=="EGF",time ==7,marker == "p.S6") %>%
              mutate(cond = "mean: 0.853"),
          sub_data_all %>% filter(cell_line == 'LY2',treatment=="iPI3K",time ==40,marker == "p.PLCg2")%>%
              mutate(cond = "0.25 percentile: 0.619"),
          sub_data_all %>% filter(cell_line == 'MACLS2',treatment=="iPI3K",time ==40,marker == "p.HER2")%>%
              mutate(cond = "0.75 percentile: 0.965")) %>% 
    select(cellID, standard,icx_bxai,cond) %>% 
    group_by(cond) %>% 
    sample_n(200) %>% 
    ungroup() %>%
    mutate(cond = factor(cond)) %>% 
    mutate(cond = factor(cond,levels = levels(cond)[c(1,3,2)])) %>%
    ggplot(aes(standard,icx_bxai)) +
    geom_point() + facet_wrap(~cond) + geom_abline(slope = 1) + coord_equal()







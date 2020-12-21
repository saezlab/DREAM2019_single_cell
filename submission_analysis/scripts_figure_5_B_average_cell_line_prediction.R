# SC4: compute the median of traiing cell-line and submit as a candidate.
# compare the predictions with participants

library(tidyverse)
library(ggpubr)

median_predictions <- read_rds("./submission_analysis/intermediate_data/sc4_combined_data.rds")

median_training <- read_csv("challenge_data/median_phospho/median_phospho_data.csv")

ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc4_ranked_teams.rds")

marker_columns <- unique(median_predictions$marker)


# compute the average cell-line by summarise the training



average_cell_line <- median_training %>%
    gather(marker,value,-cell_line, -treatment, -time) %>%
    group_by(treatment,time,marker) %>%
    summarise(ave_cell_line = median(value,na.rm=TRUE),
              sd_cell_line = sd(value,na.rm=TRUE)) 



average_cell_line %>% group_by(marker) %>% summarise(msd = mean(sd_cell_line)) %>% arrange(desc(msd))


# plot the response of the average cell-line
average_cell_line %>%
    ggplot(aes(time,ave_cell_line,col=treatment)) + 
    geom_point() + 
    geom_line()  + 
    facet_wrap(~marker)
    #facet_wrap(~marker, scales = "free_y")


# compute the RMSE score for all participant and the ave_cell_line

median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    group_by(cell_line, treatment, marker) %>%
    summarise_at(c("ave_cell_line",as.character(ranked_teams)),~ sqrt(sum((standard - .)^2) / n())) %>%
    ungroup() %>%
    summarise_at(c("ave_cell_line",as.character(ranked_teams)),mean) 
# ave_cell_line: 0.4283408


# Correlation:
median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    group_by(cell_line, treatment, marker) %>%
    gather(team, prediction,c("ave_cell_line",as.character(ranked_teams))) %>%
    #summarise_at(c("ave_cell_line",as.character(ranked_teams)),~ 1 - sum((standard - .)^2)/sum((standard - mean(standard))^2) ) %>%
    group_by(cell_line, treatment, marker,team) %>%
    summarise(cor =  cor(standard, prediction)) %>%
    ungroup() %>%
    ggplot() + geom_boxplot(aes(team,cor))+ coord_cartesian(ylim = c(-1,1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))





median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    select(cell_line,treatment,time, marker, standard,icx_bxai,ave_cell_line) %>%
    gather(source,value,standard,icx_bxai,ave_cell_line) %>%
    #filter(treatment=="EGF",cell_line=="CAL120") %>%
    filter(treatment=="iMEK") %>%
    filter(marker %in% c("p.ERK","p.S6","p.Akt.Ser473.")) %>%
    ggplot(aes(time, value,col=source)) +
    geom_point() +
    geom_line() +
    facet_grid(cell_line~marker,scales = "free_y") +
    theme_bw() +
    scale_color_manual(values = c(ave_cell_line = "grey",
                                  standard = RColorBrewer::brewer.pal(7,"Dark2")[[1]],
                                  icx_bxai = RColorBrewer::brewer.pal(7,"Dark2")[[2]]))




gg_Data <-  median_predictions %>%
    left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    #select(cell_line,treatment,time, marker, standard,icx_bxai,ave_cell_line,) %>%
    gather(source,value,standard,icx_bxai,ave_cell_line,as.character(ranked_teams[1:13]))  %>%
    #filter(treatment=="EGF",cell_line=="CAL120") %>%
    #filter(cell_line=="CAL120") %>%
    filter(treatment=="EGF") %>%
    filter(marker %in% c("p.ERK","p.S6","p.Akt.Ser473.","p.p90RSK")) 


gg1 <- ggplot() +
    geom_point(data = filter(gg_Data,!source %in% c("standard","average_cell_line")),
               aes(time, value,col="predictions")) +
    geom_line(data = filter(gg_Data,!source %in% c("standard","average_cell_line")),
              aes(time, value,group=source,col="predictions"),col="grey80") +
    geom_point(data = filter(gg_Data,source=="standard"),aes(time, value, col="data")) +
    geom_line(data = filter(gg_Data,source=="standard"),aes(time, value, col="data")) +
    geom_point(data = filter(gg_Data,source=="ave_cell_line"),aes(time, value, col="average cell line")) +
    geom_line(data = filter(gg_Data,source=="ave_cell_line"),aes(time, value,col="average cell line")) +
    facet_grid(cell_line+treatment~marker) +
    theme_bw() +
     scale_color_manual(values = c(predictions = "grey80",
                                   `average cell line` = RColorBrewer::brewer.pal(7,"Dark2")[[1]],
                                   data = RColorBrewer::brewer.pal(7,"Dark2")[[2]]))  +
    theme(panel.grid = element_blank(),
          legend.title = element_blank()) +
    xlab("time [min]") + ylab("signal internsity")
# scale_color_brewer(palette = "Dark2")
print(gg1)

ggsave("publication/figures/figure5/sc4_prediction_examples_wlegend_cond_EGF.pdf",plot = gg1, height = 8,width = 5)

gg1_legend <- get_legend(gg1)

ggsave("publication/figures/figure5/sc4_prediction_examples_legend_only.pdf",plot = gg1_legend, height = 1,width = 2)

# without legend
ggplot() +
    geom_point(data = filter(gg_Data,!source %in% c("standard","average_cell_line")),
               aes(time, value,col="predictions")) +
    geom_line(data = filter(gg_Data,!source %in% c("standard","average_cell_line")),
              aes(time, value,group=source,col="predictions"),col="grey80") +
    geom_point(data = filter(gg_Data,source=="standard"),aes(time, value, col="data")) +
    geom_line(data = filter(gg_Data,source=="standard"),aes(time, value, col="data")) +
    geom_point(data = filter(gg_Data,source=="ave_cell_line"),aes(time, value, col="average cell line")) +
    geom_line(data = filter(gg_Data,source=="ave_cell_line"),aes(time, value,col="average cell line")) +
    facet_grid(cell_line~marker) +
    theme_bw() +
    guides(color="none") +
    scale_color_manual(values = c(predictions = "grey80",
                                  `average cell line` = RColorBrewer::brewer.pal(7,"Dark2")[[1]],
                                  data = RColorBrewer::brewer.pal(7,"Dark2")[[2]]))  +
    theme(panel.grid = element_blank(),
          legend.title = element_blank()) +
    xlab("time [min]") + ylab("signal internsity")
# scale_color_brewer(palette = "Dark2")


ggsave("publication/figures/figure5/sc4_prediction_examples_cond_EGF.pdf", height = 8,width = 5)




## Emphesizing the top 3 teams: 
team_rank = tibble(teams = as.character(ranked_teams), position = 1:length(ranked_teams))


gg_Data <-  median_predictions %>%
    left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    gather(source,value,standard,icx_bxai,ave_cell_line,as.character(ranked_teams))  %>%
    #filter(treatment=="EGF",cell_line=="CAL120") %>%
    filter(cell_line=="CAMA1") %>%
    #filter(treatment=="EGF") %>%
    filter(marker %in% c("p.ERK","p.S6","p.Akt.Ser473.","p.p90RSK")) 

gg_Data <- gg_Data %>% left_join(team_rank, by = c("source"="teams"))

gg1 <- ggplot() +
    #geom_point(data = filter(gg_Data,! source %in% c("standard","average_cell_line")),
    #           aes(time, value,col="predictions",alpha = position<3 )) +
    
    geom_line(data = filter(gg_Data,!source %in% c("standard","average_cell_line",position>3)),
              aes(time, value,group=source),col="grey90",alpha=0.9) +
    geom_line(data = filter(gg_Data,!source %in% c("standard","average_cell_line"),position<=3),
              aes(time, value,group=source),col="black") +
    
    geom_point(data = filter(gg_Data,source=="standard"),aes(time, value, col="data")) +
    geom_line(data = filter(gg_Data,source=="standard"),aes(time, value, col="data")) +
    
    #geom_point(data = filter(gg_Data,source=="ave_cell_line"),aes(time, value, col="average cell line")) +
    geom_line(data = filter(gg_Data,source=="ave_cell_line"),aes(time, value, col="average cell line"),size=1) +
    
    facet_grid(cell_line+treatment~marker) +
    theme_bw() +
    scale_color_manual(values = c(predictions = "black",
                                  `average cell line` = RColorBrewer::brewer.pal(7,"Dark2")[[1]],
                                  data = RColorBrewer::brewer.pal(7,"Dark2")[[2]]))  +
    theme(panel.grid = element_blank(),
          legend.title = element_blank()) +
    xlab("time [min]") + ylab("signal internsity")
# scale_color_brewer(palette = "Dark2")
print(gg1)

ggsave("publication/figures/figure5/sc4_CAMA1_examples_wlegend_cond_ALL.pdf",plot = gg1, height = 8,width = 8)

####

# compare the RMSE for the top 3 teams and average cell-line

RMSE_all <- median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    group_by(cell_line, treatment, marker) %>%
    summarise_at(c("ave_cell_line",as.character(ranked_teams)),~ sqrt(sum((standard - .)^2) / n())) %>%
    ungroup() %>%
    select(cell_line,treatment,marker,ave_cell_line,icx_bxai,AMbeRland,orangeballs) %>% 
    gather(source, RMSE,ave_cell_line,icx_bxai,AMbeRland,orangeballs) %>%
    mutate(source = ifelse(source=="ave_cell_line","average","predictions"))
    

library(ggpubr)

# (1) based on cell lines 

RMSE_all %>% ggplot() + 
    geom_boxplot(aes(cell_line,RMSE,group=paste(cell_line,source),fill=source),draw_quantiles = 0.5) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust = 1)) + labs(fill = "")

ggsave("publication/figures/figure5/sc4_prediction_vs_averageCL_celllines.pdf", height = 3,width = 5)

### 
# (2) based on markers
# markers_prdered by rmse

makers_by_rmse <- RMSE_all %>% group_by(marker) %>%
    summarise(med_rmse = median(RMSE)) %>% arrange(desc(med_rmse)) %>% pull(marker)

RMSE_all %>% mutate(marker = factor(marker,levels = makers_by_rmse)) %>%
    ggplot() + 
    geom_boxplot(aes(marker,RMSE,group=paste(marker,source),fill=source),draw_quantiles = 0.5) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90,hjust = 1))

ggsave("publication/figures/figure5/sc4_prediction_vs_averageCL_markers.pdf", height = 3,width = 10)


# which sites are predicted better than average
rmse_ttest_markers <- RMSE_all %>% group_by(marker) %>% nest() %>%
    mutate(ttest = map(data,function(df){
    
        x = df %>% filter(source == "average") %>% pull(RMSE)
        y = df %>% filter(source != "average") %>% pull(RMSE)
        t = t.test(x = x,y=y)
        broom::tidy(t) %>% rename(rmse_diff  = estimate,average_rmse = estimate1, top3_rmse = estimate2)
        
    })) %>% select(-data) %>% unnest(ttest)


makers_by_ttest <- rmse_ttest_markers %>% arrange(desc(rmse_diff)) %>% pull(marker)
    

RMSE_all <- RMSE_all %>% mutate(marker = factor(marker,levels = makers_by_ttest)) 
RMSE_all %>%    ggplot() + 
    geom_boxplot(aes(marker,RMSE,group=paste(marker,source),fill=source),draw_quantiles = 0.5) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90,hjust = 1))

makers_by_ttest <- rmse_ttest_markers %>% arrange(p.value) %>% pull(marker)
p <- RMSE_all %>% filter(marker %in% makers_by_ttest[1:10] ) %>%
    ggboxplot(data = ., x = "source", y = "RMSE",
              fill = "source", 
              add = "jitter",
              add.params = list(size=0.4),
              facet.by = "marker", short.panel.labs = FALSE) + ylim(0,1.6)
# Use only p.format as label. Remove method name.
p + stat_compare_means(comparisons = list(c("average","predictions")), label = "p.format",label.y = 1.25)





## (3) based on treatment
RMSE_all %>% ggplot() + 
    geom_boxplot(aes(treatment,RMSE,group=paste(treatment,source),fill=source),draw_quantiles = 0.5) +
    theme_bw()




gg_Data <-  median_predictions %>%
    left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    #select(cell_line,treatment,time, marker, standard,icx_bxai,ave_cell_line,) %>%
    gather(source,value,standard,icx_bxai,ave_cell_line,as.character(ranked_teams))  %>%
    #filter(treatment=="EGF",cell_line=="CAL120") %>%
    #filter(cell_line=="CAL120") %>%
    filter(treatment=="EGF") %>%
    filter(marker %in% c("p.p53","p.MKK3.MKK6","p.AMPK.","p.SRC","p.4EBP1","p.BTK")) 

gg_Data <- gg_Data %>% left_join(team_rank, by = c("source"="teams"))

gg1 <- ggplot() +
    #geom_point(data = filter(gg_Data,! source %in% c("standard","average_cell_line")),
    #           aes(time, value,col="predictions",alpha = position<3 )) +
    geom_line(data = filter(gg_Data,!source %in% c("standard","average_cell_line")),
              aes(time, value,group=source,col="predictions",alpha =  position<3 ),col="black") +
    geom_point(data = filter(gg_Data,source=="standard"),aes(time, value, col="data")) +
    geom_line(data = filter(gg_Data,source=="standard"),aes(time, value, col="data")) +
    geom_point(data = filter(gg_Data,source=="ave_cell_line"),aes(time, value, col="average cell line")) +
    geom_line(data = filter(gg_Data,source=="ave_cell_line"),aes(time, value,col="average cell line")) +
    facet_grid(cell_line+treatment~marker,scales = "free_y") +
    theme_bw() +
    scale_color_manual(values = c(predictions = "black",
                                  `average cell line` = RColorBrewer::brewer.pal(7,"Dark2")[[1]],
                                  data = RColorBrewer::brewer.pal(7,"Dark2")[[2]]))  +
    theme(panel.grid = element_blank(),
          legend.title = element_blank()) +
    xlab("time [min]") + ylab("signal internsity")
# scale_color_brewer(palette = "Dark2")
print(gg1)






median_predictions %>%
    left_join(average_cell_line, by = c("treatment", "time", "marker"))  %>%
    select(cell_line, treatment,time,marker,ave_cell_line,icx_bxai,AMbeRland,orangeballs) %>%
    filter(cell_line=="CAMA1",marker=="p.p90RSK")








#### Other interesting cases based on the heatmap

gg_Data <-  median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    #select(cell_line,treatment,time, marker, standard,icx_bxai,ave_cell_line,) %>%
    gather(source,value,standard,icx_bxai,ave_cell_line,as.character(ranked_teams)) %>%
    #filter(treatment=="EGF",cell_line=="CAL120") %>%
    filter(treatment=="iMEK") %>%
    filter(marker %in% c("Ki.67","p.CREB","p.GSK3b","p.RB","p.STAT5")) 



ggplot() +
    geom_point(data = filter(gg_Data,!source %in% c("standard","average_cell_line")),
               aes(time, value,col="predictions")) +
    geom_line(data = filter(gg_Data,!source %in% c("standard","average_cell_line")),
              aes(time, value,group=source,col="predictions"),col="grey80") +
    geom_point(data = filter(gg_Data,source=="standard"),aes(time, value, col="data")) +
    geom_line(data = filter(gg_Data,source=="standard"),aes(time, value, col="data")) +
    geom_point(data = filter(gg_Data,source=="ave_cell_line"),aes(time, value, col="average cell line")) +
    geom_line(data = filter(gg_Data,source=="ave_cell_line"),aes(time, value,col="average cell line")) +
    facet_grid(cell_line~marker) +
    theme_bw() +
    guides(color="none") +
    scale_color_manual(values = c(predictions = "grey80",
                                  `average cell line` = RColorBrewer::brewer.pal(7,"Dark2")[[1]],
                                  data = RColorBrewer::brewer.pal(7,"Dark2")[[2]]))  +
    theme(panel.grid = element_blank(),
          legend.title = element_blank()) +
    xlab("time [min]") + ylab("signal internsity")

ggsave("publication/figures/figure5/sc4_prediction_examples_high_variance.pdf", height = 5,width = 7)




### based on Marco's feedback




gg_Data <-  median_predictions %>% 
    
    gather(source,value,standard,icx_bxai,as.character(ranked_teams)) %>%
    #filter(treatment=="EGF",cell_line=="CAL120") %>%
    filter(treatment=="EGF") %>%
    filter(marker %in% c("p.ERK","p.S6","p.Akt.Ser473.")) 


gg1 <- ggplot() +
    #geom_point(data = filter(gg_Data,!source %in% c("standard","average_cell_line")),
    #           aes(time, value,col="predictions")) +
    geom_line(data = filter(gg_Data,!source %in% c("standard","average_cell_line")),
              aes(time, value,group=source,col="predictions"),col="black",alpha=0.2) +
    geom_point(data = filter(gg_Data,source=="standard"),aes(time, value, col="data")) +
    geom_line(data = filter(gg_Data,source=="standard"),aes(time, value, col="data")) +
    facet_grid(cell_line~marker) +
    theme_bw() +
    scale_color_manual(values = c(predictions = "grey80",
                                  `average cell line` = RColorBrewer::brewer.pal(7,"Dark2")[[1]],
                                  data = RColorBrewer::brewer.pal(7,"Dark2")[[2]]))  +
    theme(panel.grid = element_blank(),
          legend.title = element_blank()) +
    xlab("time [min]") + ylab("signal internsity")
# scale_color_brewer(palette = "Dark2")


##### How different are the cell-lines from the average?


# plot the difference between the mean and the cell-lines:
median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    select(cell_line, treatment,time, marker,standard,ave_cell_line) %>%
    ggplot() + geom_line(aes(time,y = standard-ave_cell_line,col=treatment)) + facet_grid(cell_line~marker)


ave_standards_diff = median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>% 
    group_by(cell_line, treatment,marker) %>%
    summarise(rmse_diff = sqrt(sum((standard-ave_cell_line)^2/n())),
              sd_data = sd(standard))

ave_standards_diff %>% arrange(sd_data)


ave_standards_diff %>% ggplot() + geom_violin(aes(cell_line,rmse_diff))
ave_standards_diff %>% ggplot() + geom_violin(aes(treatment,rmse_diff))

ave_standards_diff %>% group_by(cell_line) %>% summarise(mean_diff= mean(rmse_diff))

ave_standards_diff %>% group_by(marker) %>%
    summarise(mean_diff= mean(rmse_diff),
              mean_sd= mean(sd_data)) %>%
    mutate(rank_sd = rank(mean_sd)) %>%  arrange(desc(mean_diff))


# this figure shows that highly variable markers are also have larger difference from the median cell-line: 
ave_standards_diff %>% group_by(marker) %>%
    summarise(mean_diff= mean(rmse_diff),
              mean_sd= mean(sd_data)) %>%
    ggplot(aes(mean_sd,mean_diff)) + geom_point() + geom_text(aes(label=marker))


ave_standards_diff %>% group_by(marker) %>% summarise(mean_diff= mean(rmse_diff)) %>% arrange(desc(mean_diff))

rmse_aov <- aov(diff~marker+treatment+cell_line,data = ave_standards_diff)
summary(rmse_aov)



#### Top team
top_standards_diff = median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>% 
    gather(team, predicted_value, c("ave_cell_line",as.character(ranked_teams)) )%>%
    group_by(cell_line, treatment,marker,team) %>%
    summarise(rmse_diff = sqrt(sum((standard-predicted_value)^2/n())),
              sd_data = sd(standard))

top_standards_diff %>% arrange(sd_data)




# this figure shows that highly variable markers are also have larger difference from the median cell-line: 
top_standards_diff %>% group_by(marker,team) %>%
    summarise(mean_diff= mean(rmse_diff),
              mean_sd= mean(sd_data)) %>% 
    filter(team %in% c("ave_cell_line","icx_bxai")) %>%
    ggplot(aes(mean_sd,mean_diff,col=team)) + geom_point() + geom_text(aes(label=marker))


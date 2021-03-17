# DREAM challenge post analysis
# What is the most difficult to predict? 
# Violin plots based on conditions - ANOVA analysis
# A. Gabor


library(tidyverse)
# SC2

combined_statistics <- read_rds("./submission_analysis/intermediate_data/sc2_stats_conditions.rds")


ranked_teams = names(combined_statistics)[7:22]


# convert the statictics to error (signed)
error_terms <- combined_statistics %>% 
    mutate_at(teams,~.-standard) %>% 
    gather(team,error,7:22)

# we go with the top 4 which are better than reference
#error_terms <- error_terms %>% filter(team %in% ranked_teams[1:4] )


fill_color <- RColorBrewer::brewer.pal(8,"Dark2")[[1]]
# mean square error by CELL-LINE:
plot_rmse_cell_line <- error_terms %>% group_by(cell_line,treatment, time,team) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    mutate(top = team %in% ranked_teams[1:4] ) %>% 
    ggplot() + geom_violin(aes(cell_line,rmse),draw_quantiles = .5,fill= fill_color)  +
    theme_bw() + 
    guides(fill = "none") + 
    labs( ylab= NULL) + xlab("") +
    theme(axis.text.x = element_text(angle = 90)) + ggtitle("Cell lines")


plot_rmse_treatment <- error_terms %>% group_by(cell_line,treatment, time,team) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    mutate(top = team %in% ranked_teams[1:4] ) %>% 
    ggplot() + geom_violin(aes(treatment,rmse),draw_quantiles = .5,fill= fill_color)  +
    theme_bw() + 
    guides(fill = "none") + 
    labs( ylab= NULL) + xlab("") + ggtitle("Treatment")

plot_rmse_time <- error_terms %>% group_by(cell_line,treatment, time,team) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    mutate(time = ifelse(time==16,17,time)) %>%
    mutate(top = team %in% ranked_teams[1:4] ) %>% 
    ggplot() + geom_violin(aes(factor(time),rmse),draw_quantiles = .5,fill= fill_color)  +
    theme_bw() + 
    guides(fill = "none") + 
    labs( ylab= NULL) + xlab("") + ggtitle("Time")


cowplot::plot_grid(plot_rmse_cell_line,plot_rmse_time,NULL, plot_rmse_treatment,align = "vh")

ggsave("./publication/figures/supp_figures/sc2_prediction_error_per_condition.pdf",width = 8,height = 6)


# with statistics:
mean_error <- error_terms %>% group_by(cell_line,treatment, time,team) %>%
    summarise(score = sqrt(sum(error^2))) %>%
    mutate(top = team %in% ranked_teams[1:4] ) 

plot_rmse_time <- mean_error %>% ggviolin(x="time",y="score",draw_quantiles = 0.5,fill= fill_color)+
    stat_compare_means(method = "anova",label.y = 27,label.x = 2) + 
    geom_hline(yintercept = median(mean_error$score), linetype = 2,color="red",size=1)+
    stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.", label.y = 25)+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60,hjust = 1))+
    ggtitle("Time points") +
    xlab("") + ylab("") + 
    coord_cartesian(ylim = c(0,27))

plot_rmse_treatment <- mean_error %>% ggviolin(x="treatment",y="score",draw_quantiles = 0.5,fill= fill_color)+
    stat_compare_means(method = "anova",label.y = 27,label.x = 1) + 
    geom_hline(yintercept = median(mean_error$score), linetype = 2,color="red",size=1)+
    stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.", label.y = 25)+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60,hjust = 1))+
    ggtitle("Treatments") +
    xlab("")+ ylab("") + 
    coord_cartesian(ylim = c(0,27))

plot_rmse_cell_line <- mean_error %>% ggviolin(x="cell_line",y="score",draw_quantiles = 0.5,fill= fill_color)+
    stat_compare_means(method = "anova",label.y = 27,label.x = 2) + 
    geom_hline(yintercept = median(mean_error$score), linetype = 2,color="red",size=1)+
    stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.", label.y = 25)+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60,hjust = 1))+
    ggtitle("Cell lines") +
    xlab("")+ ylab("condition score") + 
    coord_cartesian(ylim = c(0,27))



cowplot::plot_grid(plot_rmse_cell_line,plot_rmse_time, plot_rmse_treatment,nrow=1,align = "vh")
ggsave("./publication/figures/supp_figures/sc2_prediction_error_per_condition.pdf",width = 12,height = 4)


# T-test: cell line group vs global
mean_error %>% ungroup() %>%
    select(cell_line,score) %>%
    group_nest(cell_line) %>%
    mutate(t_test = map(data,function(x){
        t.test(x = x$score,y = mean_error$score)
    })) %>%
    mutate(res = map(t_test,broom::tidy)) %>%
    select(cell_line,res) %>% unnest(res)%>%
    rename(groupMeanScore = "estimate1",
           globalMeanScore = "estimate2",
           deltaScore = "estimate") %>%
    mutate(relScore = deltaScore/globalMeanScore) %>%
    select(-method,-alternative)


    
## TOP 4 teams only
mean_error <- error_terms %>% group_by(cell_line,treatment, time,team) %>%
    summarise(score = sqrt(sum(error^2))) %>%
    filter(team %in% ranked_teams[1:4] ) 

plot_rmse_time <- mean_error %>% ggviolin(x="time",y="score",draw_quantiles = 0.5,fill= fill_color)+
    stat_compare_means(method = "anova",label.y = 27,label.x = 2) + 
    geom_hline(yintercept = median(mean_error$score), linetype = 2,color="red",size=1)+
    stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.", label.y = 25)+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60,hjust = 1))+
    ggtitle("Time points") +
    xlab("") + ylab("") + 
    coord_cartesian(ylim = c(0,27))

plot_rmse_treatment <- mean_error %>% ggviolin(x="treatment",y="score",draw_quantiles = 0.5,fill= fill_color)+
    stat_compare_means(method = "anova",label.y = 27,label.x = 1) + 
    geom_hline(yintercept = median(mean_error$score), linetype = 2,color="red",size=1)+
    stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.", label.y = 25)+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60,hjust = 1))+
    ggtitle("Treatments") +
    xlab("")+ ylab("") + 
    coord_cartesian(ylim = c(0,27))

plot_rmse_cell_line <- mean_error %>% ggviolin(x="cell_line",y="score",draw_quantiles = 0.5,fill= fill_color)+
    stat_compare_means(method = "anova",label.y = 27,label.x = 2) + 
    geom_hline(yintercept = median(mean_error$score), linetype = 2,color="red",size=1)+
    stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.", label.y = 25)+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60,hjust = 1))+
    ggtitle("Cell lines") +
    xlab("")+ ylab("condition score") + 
    coord_cartesian(ylim = c(0,27))



cowplot::plot_grid(plot_rmse_cell_line,plot_rmse_time, plot_rmse_treatment,nrow=1,align = "vh")
ggsave("./publication/figures/supp_figures/sc2_prediction_error_per_condition_top4.pdf",width = 12,height = 4)


# T-test: cell line group vs global
mean_error %>% ungroup() %>%
    select(cell_line,score) %>%
    group_nest(cell_line) %>%
    mutate(t_test = map(data,function(x){
        t.test(x = x$score,y = mean_error$score)
    })) %>%
    mutate(res = map(t_test,broom::tidy)) %>%
    select(cell_line,res) %>% unnest(res)%>%
    rename(groupMeanScore = "estimate1",
           globalMeanScore = "estimate2",
           deltaScore = "estimate") %>%
    mutate(relScore = deltaScore/globalMeanScore) %>%
    select(-method,-alternative) %>% 












# sort by max of error:
sorted_stats <- error_terms %>% group_by(stat_variable) %>% summarise(max = max(abs(error))) %>%
    arrange(desc(max)) 
error_terms %>% 
    mutate(stat_variable = factor(stat_variable,levels = sorted_stats$stat_variable)) %>%
    filter(stat_variable %in% sorted_stats$stat_variable[1:20]) %>%
    ggplot() + geom_boxplot(aes(stat_variable,error))  +
    theme_bw() + 
    guides(fill = "none") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90)) +     
    labs( ylab= NULL) + xlab("variables")
    
ggsave("./publication/figures/supp_figures/sc2_prediction_error_per_statistics.pdf",width = 5,height = 6)

error_terms %>% ggplot() + geom_histogram(aes(error),bins = 100) + theme_bw()

ggsave("./publication/figures/supp_figures/sc2_prediction_error_distribution.pdf",width = 3,height = 3)


error_terms %>% group_by(cell_line,treatment, time) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    aov(formula = rmse ~ cell_line+treatment+time, data = .) %>% summary()

error_terms %>% group_by(cell_line,treatment, time) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    aov(formula = rmse ~ cell_line, data = .) %>% summary()

error_terms %>% group_by(cell_line,treatment, time) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    aov(formula = rmse ~ treatment, data = .) %>% summary()

error_terms %>% group_by(cell_line,treatment, time) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    aov(formula = rmse ~ time, data = .) %>% summary()


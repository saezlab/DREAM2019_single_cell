# DREAM challenge post analysis
# What is the most difficult to predict? 
# Violin plots based on conditions - ANOVA analysis
# A. Gabor


library(tidyverse)
# SC2

combined_statistics <- read_rds("./submission_analysis/intermediate_data/sc3_stats_conditions.rds")


teams = names(combined_statistics)[7:20]


# convert the statictics to error (signed)
error_terms <- combined_statistics %>% 
    mutate_at(teams,~.-standard) %>% 
    gather(team,error,7:20)


# fix timing:
error_terms <- error_terms %>%
    mutate(time = ifelse(time ==14, 13, time)) %>% 
    mutate(time = ifelse(time ==16, 17, time)) 

fill_color <- RColorBrewer::brewer.pal(8,"Dark2")[[1]]
# mean square error by CELL-LINE:
plot_rmse_cell_line <- error_terms %>% group_by(cell_line,treatment, time) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    ggplot() + geom_violin(aes(cell_line,rmse),draw_quantiles = .5,fill= fill_color)  +
    theme_bw() + 
    guides(fill = "none") + 
    labs( ylab= NULL) + xlab("Cell Lines") +
    theme(axis.text.x = element_text(angle = 90))


plot_rmse_time <- error_terms %>% group_by(cell_line,treatment, time) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    ggplot() + geom_violin(aes(factor(as.numeric(time)),rmse),draw_quantiles = .5,fill= fill_color)  +
    theme_bw() + 
    guides(fill = "none") + 
    labs( ylab= NULL) + xlab("Time")


cowplot::plot_grid(plot_rmse_cell_line,plot_rmse_time,rel_widths = c(3,1))

ggsave("./publication/figures/supp_figures/sc3_prediction_error_per_condition.pdf",width = 8,height = 3)


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


error_terms %>% group_by(cell_line, time) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    aov(formula = rmse ~ cell_line+time, data = .) %>% summary()

error_terms %>% group_by(cell_line,treatment, time) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    aov(formula = rmse ~ cell_line, data = .) %>% summary()

error_terms %>% group_by(cell_line,treatment, time) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    aov(formula = rmse ~ treatment, data = .) %>% summary()

error_terms %>% group_by(cell_line,treatment, time) %>%
    summarise(rmse = sqrt(sum(error^2)/n())) %>%
    aov(formula = rmse ~ time, data = .) %>% summary()


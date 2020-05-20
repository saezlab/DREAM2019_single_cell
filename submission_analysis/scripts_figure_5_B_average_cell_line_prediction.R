# SC4: compute the median of traiing cell-line and submit as a candidate.
# compare the predictions with participants

library(tidyverse)
library(ggpubr)

median_predictions <- read_rds("./submission_analysis/intermediate_data/sc4_combined_data.rds")

median_training <- read_csv("challenge_data/median_phospho/median_phospho_data.csv")

ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc4_ranked_teams.rds")

marker_columns <- unique(median_predictions$marker)


# compute teh average cell-line by summarise the training

average_cell_line <- median_training %>% 
    group_by(treatment,time) %>%
    summarise_at(marker_columns,median,na.rm=TRUE) %>%
    gather(marker, ave_cell_line, marker_columns) 

# plot the response of the average cell-line
average_cell_line %>%
    ggplot(aes(time,ave_cell_line,col=treatment)) + 
    geom_point() + 
    geom_line()  + 
    facet_wrap(~marker)


# compute the RMSE score for all participant and the ave_cell_line

median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    group_by(cell_line, treatment, marker) %>%
    summarise_at(c("ave_cell_line",as.character(ranked_teams)),~ sqrt(sum((standard - .)^2) / n())) %>%
    ungroup() %>%
    summarise_at(c("ave_cell_line",as.character(ranked_teams)),mean) 
# ave_cell_line: 0.4283408


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




gg_Data <-  median_predictions %>% left_join(average_cell_line, by = c("treatment", "time", "marker")) %>%
    #select(cell_line,treatment,time, marker, standard,icx_bxai,ave_cell_line,) %>%
    gather(source,value,standard,icx_bxai,ave_cell_line,as.character(ranked_teams)) %>%
    #filter(treatment=="EGF",cell_line=="CAL120") %>%
    filter(treatment=="EGF") %>%
    filter(marker %in% c("p.ERK","p.S6","p.Akt.Ser473.")) 


gg1 <- ggplot() +
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
     scale_color_manual(values = c(predictions = "grey80",
                                   `average cell line` = RColorBrewer::brewer.pal(7,"Dark2")[[1]],
                                   data = RColorBrewer::brewer.pal(7,"Dark2")[[2]]))  +
    theme(panel.grid = element_blank(),
          legend.title = element_blank()) +
    xlab("time [min]") + ylab("signal internsity")
# scale_color_brewer(palette = "Dark2")


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




#### Other interesting casess based on the heatmap

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



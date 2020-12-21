# find what is the most difficult to predict, like in SC1
library(tidyverse)
library(ggplot2)

rmse_data <- read_rds("./submission_analysis/intermediate_data/sc4_rmse_conditions.rds")
median_data <- read_rds("./submission_analysis/intermediate_data/sc4_combined_data.rds")


# accuracy is computed from the RMSE
rmse_data <- rmse_data %>% gather(team, rmse,-1:-4)


# accuracy and diversity of the prediction errors:
# computed for each condition (cl, tr) and marker. Median taken across teams. 

rmse_summary_data <- rmse_data %>% 
    group_by(cell_line, treatment, marker) %>%
    summarise(med_rmse = median(rmse),
              QCV_rmse = (quantile(rmse,0.75)-quantile(rmse,0.25))/quantile(rmse,0.5)) 

# Diversity is computed from the medians
# computed for each condition (cl, tr, time) and marker. Computed across teams. 
# finally averaged over time
div_summary_data <- median_data %>%  gather(team, pred_signal,-1:-5) %>%
    group_by(cell_line, treatment, marker, time) %>%
    summarise(QCV_pred = (quantile(pred_signal,0.75)-quantile(pred_signal,0.25))/quantile(pred_signal,0.5),
              IQR_pred = quantile(pred_signal,0.75)-quantile(pred_signal,0.25)) %>%
    group_by(cell_line, treatment, marker) %>%
    summarise(QCV_pred = mean(QCV_pred),
              IQR_pred = mean(IQR_pred))


custom_colors = colorRampPalette(RColorBrewer::brewer.pal(9,"Greens"))


# updated: accuracy and QCV of predictions!
gg_main <- rmse_summary_data %>% full_join(div_summary_data, by = c("cell_line", "treatment", "marker"))  %>% 
    ggplot(aes(marker,treatment)) +
    geom_point(aes(size=med_rmse, fill = QCV_pred),na.rm = T,shape=21) +
    theme_bw() +
    theme(axis.text.x = element_text(hjust=0,angle = 90)) +
    xlab("") + ylab("") + 
    scale_fill_gradientn(colors = custom_colors(100)) +
    scale_size_area(max_size = 5)+
    facet_grid(cell_line~.) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    coord_equal()
ggsave("./publication/figures/figure5/conditional_prediction_accuracy_QCV_prediction.pdf",plot = gg_main,width = 7,height = 7)


# summarise RMSE by marker and treatment

fill_color <- RColorBrewer::brewer.pal(8,"Dark2")[[1]]

ordered_markers = rmse_markerwise_plot <- rmse_summary_data %>% group_by(marker) %>%
    summarise(rmse = mean(med_rmse)) %>% arrange(desc(rmse)) %>%
    mutate(marker = factor(marker,levels = marker)) %>% pull(marker)

rmse_markerwise_plot <- rmse_summary_data %>% group_by(marker) %>%
    summarise(rmse = mean(med_rmse)) %>% arrange(desc(rmse)) %>%
    mutate(marker = factor(marker,levels = marker)) %>%
    ggplot() + geom_col(aes(marker,rmse)) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))

rmse_markerwise_plot_violin <- rmse_summary_data %>% 
    mutate(marker = factor(marker,levels = ordered_markers)) %>%
    ggplot() +
    geom_violin(aes(marker,med_rmse),draw_quantiles = 0.5,fill=fill_color) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,hjust = 1)) +ylab("median RMSE")


ggsave("./publication/figures/figure5/rmse_by_markers.pdf",height = 3,width = 5,plot = rmse_markerwise_plot_violin)

rmse_cellline_plot <- rmse_summary_data %>% group_by(cell_line) %>%
    summarise(rmse = mean(med_rmse)) %>% 
    ggplot() + geom_col(aes(cell_line,rmse)) + theme_bw()

rmse_summary_data %>% 
    ggplot() + geom_violin(aes(cell_line,med_rmse)) + theme_bw()


rmse_treatment_plot <- rmse_summary_data %>% group_by(cell_line,treatment) %>%
    summarise(rmse = mean(med_rmse)) %>% 
    ggplot() + geom_col(aes(treatment,rmse)) + theme_bw() +facet_grid(~cell_line)

rmse_summary_data %>% 
    ggplot() + geom_violin(aes(treatment,med_rmse)) + theme_bw() +facet_grid(~cell_line)


rmse_aov <- aov(med_rmse~marker+treatment+cell_line,data = rmse_summary_data)
summary(rmse_aov)




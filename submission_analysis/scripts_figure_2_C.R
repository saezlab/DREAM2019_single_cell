# scripts_figure_2_C
# plots the diversity of the predictions and prediction accuracy by conditons.

# accuracy of the predictions is shown by the median RMSE of the teams.
# diversity of the predictions: is the QCV of the medians
# diversity of the errors: 




submission_folder <- "./submission_data/final_round/SC1"
rmse_data <- read_rds("./submission_analysis/intermediate_data/sc1_rmse_conditions.rds")
median_data <- read_rds("./submission_analysis/intermediate_data/sc1_median_conditions.rds")


# accuracy is computed from the RMSE
rmse_data <- rmse_data %>% gather(team, rmse,-1:-4) %>%
    # fix time for plotting
    mutate(time = ifelse(time==18,17,time)) %>%
    mutate(time = ifelse(time==14,13,time))


# accuracy and diversity of the prediction errors:
# computed for each condition (cl, tr, time) and marker. Median taken across teams. 
# then we average over time. 
rmse_summary_data <- rmse_data %>% 
    group_by(cell_line, treatment, marker, time) %>%
    summarise(med_rmse = median(rmse),
              QCV_rmse = (quantile(rmse,0.75)-quantile(rmse,0.25))/quantile(rmse,0.5)) %>%
    group_by(cell_line, treatment, marker) %>%  # take the average over time. 
    summarise(med_rmse = median(med_rmse),
              QCV_rmse = mean(QCV_rmse))
    #gather(grp,value,med_rmse,mad_rmse) %>% 

# Diversity is computed from the medians
# computed for each condition (cl, tr, time) and marker. Computed across teams. 
# then we average over time. 
div_summary_data <- median_data %>%  gather(team, pred_signal,-1:-4) %>%
    # fix time for plotting
    mutate(time = ifelse(time==18,17,time)) %>%
    mutate(time = ifelse(time==14,13,time)) %>%
    group_by(cell_line, treatment, marker, time) %>%
    summarise(QCV_pred = (quantile(pred_signal,0.75)-quantile(pred_signal,0.25))/quantile(pred_signal,0.5),
              IQR_pred = quantile(pred_signal,0.75)-quantile(pred_signal,0.25)) %>%
    group_by(cell_line, treatment, marker) %>%
    summarise(QCV_pred = mean(QCV_pred),
              IQR_pred = mean(IQR_pred))


custom_colors = colorRampPalette(RColorBrewer::brewer.pal(9,"Greens"))


# updated: accuracy and QCV of predictions!
rmse_summary_data %>% full_join(div_summary_data, by = c("cell_line", "treatment", "marker"))  %>% 
    ggplot(aes(marker,treatment)) +
    geom_point(aes(size=med_rmse, fill = QCV_pred),na.rm = T,shape=21) +
    theme_bw() +
    theme(axis.text.x = element_text(hjust=0,angle = 90)) +
    xlab("") + ylab("") + 
    scale_fill_gradientn(colors = custom_colors(100)) +
    scale_size_area(max_size = 5)+
    facet_grid(~cell_line) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    coord_equal()
ggsave("./publication/figures/figure2/conditional_prediction_accuracy_QCV_prediction.pdf",width = 8.3,height = 7.4)



rmse_summary_data %>% full_join(div_summary_data, by = c("cell_line", "treatment", "marker"))  %>% 
    ggplot(aes(marker,treatment)) +
    geom_point(aes(size=med_rmse, fill = IQR_pred),na.rm = T,shape=21) +
    theme_bw() +
    theme(axis.text.x = element_text(hjust=0,angle = 90)) +
    xlab("") + ylab("") + 
    scale_fill_gradientn(colors = custom_colors(100)) +
    scale_size_area(max_size = 5)+
    facet_grid(~cell_line) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    coord_equal()





rmse_summary_data %>% full_join(div_summary_data, by = c("cell_line", "treatment", "marker")) %>%
    filter(cell_line == "HCC2218") %>% arrange(med_rmse) %>% print(n=30)


rmse_summary_data %>% ggplot() + geom_violin(aes(marker,med_rmse))

# Summary tables of the rmse scores:
rmse_summary_data %>% group_by(marker) %>% summarise(median = median(med_rmse))


rmse_summary_data %>% group_by(cell_line) %>% 
    summarise(median = median(med_rmse)) %>%
    arrange(median)





# Chekc if diversity corelates with accuracy
rmse_summary_data %>% full_join(div_summary_data, by = c("cell_line", "treatment", "marker")) %>%
        ggplot() + geom_point(aes(med_rmse/mean(med_rmse),QCV_pred/mean(QCV_pred))) +facet_grid(marker~.) +
    coord_equal()

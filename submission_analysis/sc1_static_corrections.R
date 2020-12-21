
library(ggplot2)
library(tidyverse)

library(ggpubr)


# Participating team names
submissions <-  readRDS("./submission_analysis/intermediate_data/sc1_ranked_teams.rds") %>% as.character()
sub_data_all <- readRDS("./submission_analysis/intermediate_data/sc1_all_NP_predictions.rds")  %>%
    mutate_if(is.character, as.factor) 


sub_data_all <- sub_data_all %>% select(1:8,all_of(submissions))


sub_data_all_long = sub_data_all %>% gather(team, prediction,all_of(submissions))

# original RMSE by conditions:
RMSE <- sub_data_all_long %>% 
    group_by(team,cell_line,treatment,time,marker)  %>% 
    summarise(RMSE = sqrt(sum((prediction-standard)^2)/n()) )

# check, we get the same result: 0.849 for icx_bxai
RMSE %>% filter(team=="icx_bxai") %>% pull(RMSE) %>% mean()

# original explained variance:
R2 <- sub_data_all_long %>% 
    group_by(team,cell_line,treatment,time,marker)  %>% 
    summarise(R2 = 1 - (sum((prediction-standard)^2))/(sum((standard-mean(standard))^2)) )


R2 %>% mutate(team = factor(team,level = submissions)) %>%
    ggplot(aes(team,R2,fill=team)) + geom_boxplot() + ylim(-1,1)

# correction I.: correaction based on "full"
# let's say if we correct based on the distance of the mean predicted marker values in 
# EGF condition

# the correction is for each participant, cell-line and 
# marker, the difference in mean between their prediction and measurement
# in the full, unperturbed condition. 

correction <- sub_data_all_long %>% filter(treatment=="full") %>%
    group_by(team, cell_line, marker) %>%
    summarise(diff = mean(standard) - mean(prediction)) %>% ungroup()

# use this correction and recompute the RMSE
RMSE_corrected <- sub_data_all_long %>%
    left_join(correction,by = c("team","marker","cell_line")) %>%
    group_by(team,cell_line,treatment,time,marker)  %>% 
    summarise(RMSE = sqrt(sum((prediction+diff-standard)^2)/n()) ) %>%
    ungroup()

RMSE_corrected %>% mutate(team = factor(team,levels = submissions)) %>%
    ggplot() + geom_violin(aes(team,RMSE))

RMSE_corrected %>% filter(team=="icx_bxai") %>% pull(RMSE) %>% mean()

# compare new and old scores: 

correction_effect <- RMSE_corrected %>% rename(corrected = RMSE) %>% 
    left_join(RMSE %>% rename(original = RMSE) ,by = c("team", "cell_line", "treatment", "time", "marker")) %>% 
    group_by(team) %>% summarise(corrected_mean = mean(corrected),
                                 original_mean = mean(original)) %>%
    mutate(diff = corrected_mean - original_mean) %>%
    ungroup() %>% 
    summarise(corrected_mean = mean(corrected_mean),
              original_mean = mean(original_mean),
              mean_reduction = mean(diff),
              sd_reduction = sd(diff))



## Compare RMSE and R2 corrected vs uncorrected: 
RMSE_corrected %>% rename(corrected = RMSE) %>% 
    left_join(RMSE %>% rename(original = RMSE) ,by = c("team", "cell_line", "treatment", "time", "marker")) %>% 
    gather(predictions,RMSE,corrected,original) %>% 
    group_by(team,predictions) %>% 
    summarise(RMSE = mean(RMSE)) %>% 
    ungroup() %>%
    ggviolin(data = .,
             x = "predictions", y = "RMSE",
          fill = "predictions",
          palette = "jco",
          draw_quantiles = .5,xlab = FALSE) +
    geom_point(size = 0.2)  +
    geom_line(aes(predictions,RMSE,group=team),size = 0.5, alpha = 0.2)  +
    stat_compare_means(comparisons = list(c("original","corrected")), 
                       label.y = c(2),
                       paired = TRUE,
                       label = "p.format",
                       method = "wilcox.test")

ggsave ("./publication/figures/figure2/RMSE_original_vs_correction.pdf",width = 4,height = 4)



R2_corrected %>% rename(corrected = R2) %>% 
    left_join(R2%>% rename(original = R2) ,by = c("team", "cell_line", "treatment", "time", "marker")) %>% 
    filter(team =="icx_bxai") %>%
    gather(predictions,R2,corrected,original) %>%
    ggviolin(data = .,
             x = "predictions", y = "R2",
             fill = "predictions",
             palette = "jco",
             draw_quantiles = .5,xlab = FALSE)+ 
    stat_compare_means(comparisons = list(c("original","corrected")), label.y = c(1.2),paired = TRUE) +coord_cartesian(ylim = c(-6,2))

ggsave ("./publication/figures/figure2/RMSE_original_vs_correction.pdf",width = 4,height = 4)

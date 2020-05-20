#####################   Edge inference point ###############################

library(tidyverse)
# SC2

combined_statistics <- read_rds("./submission_analysis/intermediate_data/sc2_stats_conditions.rds")


covariance <- combined_statistics %>% gather(source,value,6:22) %>% filter(grepl("cov_",stat_variable)) %>%
    separate(col = stat_variable,into = c(NA,"var1","var2"),sep = "_",remove = FALSE) 

# create variance terms
variance_1 <- covariance %>% filter(var1==var2) %>% select(-var2)
variance_2 <- covariance %>% filter(var1==var2) %>% select(-var1)

# merge the covariance matrix and variance terms then calculate correlations
correlation <- covariance %>% left_join(variance_1, 
                                        by = c("cond_id", "cell_line", "treatment", "time", "source", "var1"),
                                        suffix = c("","_var1")) %>%
    left_join(variance_2, 
              by = c("cond_id", "cell_line", "treatment", "time", "source", "var2"),
              suffix = c("","_var2")) %>%
    select(-stat_variable_var1,-stat_variable_var2) %>%
    mutate(corr = value/sqrt(value_var1)/sqrt(value_var2))


hist(correlation$corr)
correlation %>% filter(var1!=var2) %>% 
    filter(source == "standard")%>%
    arrange(desc(corr)) %>%
    ggplot() + geom_density(aes(corr),fill=fill_color) +
    theme_bw() +
    xlab("correlations")

ggsave("./publication/figures/supp_figures/sc2_marker_correlation_data.pdf",width = 3,height = 3)


# create precision recall



ground_T_correlations <- correlation %>% filter(var1!=var2) %>% filter(source=="standard")
teams_correlation <- correlation %>% filter(var1!=var2) %>% filter(source!="standard")

# rearrange so that the correlation in the data is in a separate column (standard)
correlation <- correlation %>% filter(var1!=var2) %>% 
    select(cell_line,treatment,time,stat_variable,source,corr) %>%
    spread(source,corr) %>%
    gather(source,corr, -cell_line,-treatment,-time,-stat_variable,-standard)




# we have to binarise the covariance matrix of the real data: 
# for this we use a fix threshold: 
eps_true = 0.4

# for the binarisation of the predicted data, we use a running threshold 
eps = seq(0,1,by = 0.1)
n_total_edge = 595
ROC = list()
for(i in 1:length(eps)){
    
    ROC[[i]] = correlation %>% group_by(cell_line,treatment,time,source) %>%
        summarise(P = sum(standard > eps_true),
                  N = sum(standard < eps_true),
                  TP = sum(standard > eps_true & corr > eps[[i]]),
                  TN = sum(standard < eps_true & corr < eps[[i]]),
                  FP = sum(standard < eps_true & corr > eps[[i]])) %>%
        group_by(cell_line,treatment,time,source) %>%
        mutate(TPR = TP/P,
               FPR = FP/(FP + TN)) %>% 
        group_by(source) %>%
        summarise(mean_TPR = mean(TPR,na.rm = TRUE),
                  max_TPR = max(TPR,na.rm = TRUE),
                  min_TPR = min(TPR,na.rm = TRUE),
                  mean_FPR = mean(FPR,na.rm = TRUE),
                  max_FPR = max(FPR,na.rm = TRUE),
                  min_FPR = min(FPR,na.rm = TRUE),
                  mean_P = mean(P)) %>%
        mutate(threshold = eps[[i]])
    
    
}


bind_rows(ROC) %>% ggplot(aes(mean_FPR,mean_TPR,,col=source)) +
    geom_point() + geom_line() + theme_bw() + 
    xlab("False positive rate (averaged across conditions)") +
    ylab("True positive rate (averaged across conditions)")


bind_rows(ROC)  %>% 
    ggplot(aes(mean_FPR,mean_TPR,,col=source)) +
    geom_pointrange(aes(ymin = min_TPR, ymax = max_TPR),alpha=0.3) +
    geom_errorbarh(aes(xmax = min_FPR, xmin = max_FPR, height = 0),alpha=.3) + 
    geom_line() +
    theme_bw() + 
    xlab("False positive rate (averaged across conditions)") +
    ylab("True positive rate (averaged across conditions)")



ggsave("./publication/figures/supp_figures/sc2_ROC_curve.pdf",width = 7,height = 5)

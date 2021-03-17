# for SC3 we check how the score depends on the used methods by t-test. 


library(tidyverse)
library(cowplot)
library(broom)
library(ggpubr)
### A) Leaderboard overview ------------------------------------------------------

submission_folder <- "./submission_data/final_round/SC3"
bootstrap_RMSE <- read_rds("./submission_analysis/intermediate_data/sc3_bootstrap_stats.rds")
conditional_error_stats <- read_rds("./submission_analysis/intermediate_data/sc3_stats_sumSquared_conditions.rds")
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc3_ranked_teams.rds")

method_usage <- readxl::read_xlsx("./challenge_methods/sc3_methods_processed.xlsx",sheet = 2)
descriptor_order <- method_usage$Question

method_usage <-  method_usage %>%
    gather(submitterId,answer,-Question) %>% spread(Question,answer)

# fix naming
ranked_teams <- factor(gsub("X.","",ranked_teams,fixed = T),levels = gsub("X.","",ranked_teams,fixed = T))
names(bootstrap_RMSE) <- gsub("X.","",names(bootstrap_RMSE),fixed = T)
names(conditional_error_stats) <- gsub("X.","",names(conditional_error_stats),fixed = T)

my_colors  = RColorBrewer::brewer.pal(8,"Dark2")

robust_score <- bootstrap_RMSE %>% 
    gather(teams,RMSE,-BS_sample) %>% 
    group_by(teams) %>% summarise(score_min = min(RMSE),
                                  score_max = max(RMSE),
                                  score_med = median(RMSE))

method_res <- method_usage %>% 
    mutate(teams = factor(submitterId,levels = levels(ranked_teams))) %>%
    select(-submitterId) %>%
    gather(method_info,method_value,-teams) %>%
    mutate(method_info = factor(method_info,levels = rev(descriptor_order))) 

# T-test:
# H0: using a specific method is as good as other method
# Halt: not using a specific method gives worse (higher) RMSE
full_join(robust_score,method_res, by="teams") %>%
    group_by(method_info) %>%
    filter(!method_info %in% c("Prior knowledge network","Ensemble of models") ) %>%
    nest() %>%
    mutate(t_res = map(data,function(df){
        
        tres <- t.test(score_med~method_value,data = df,alternative = "greater")
        tidy(tres)
    } )) %>% unnest(t_res) %>% arrange(p.value)



full_join(robust_score,method_res, by="teams") %>%
    group_by(method_info,method_value) %>%
    mutate(method_value = ifelse(method_value=="y","used","not used")) %>%
    filter(!method_info %in% c("Prior knowledge network","Ensemble of models") ) %>%
    summarise(mean_score = mean(score_med),
              min_score = min(score_med),
              max_score = max(score_med)) %>%
    ggplot(aes(method_value,mean_score,fill=method_value)) + 
    geom_bar(stat = "identity",col="black") +
    
    geom_errorbar(aes(ymin= min_score, ymax=max_score)) + 
    facet_wrap(~method_info) + theme_bw() + xlab('') + ylab('Score (median,cov)') + guides(fill=FALSE)
ggsave("./publication/figures/figure3/sc2_methods_score_improvement.pdf")


full_join(robust_score,method_res, by="teams") %>%
    group_by(method_info,method_value) %>%
    mutate(method_value = ifelse(method_value=="y","infers population statistic","other method")) %>%
    filter(method_info %in% c("Infer statistics") ) %>%
    ggplot(aes(method_value,score_med)) + 
    geom_boxplot(aes(col=method_value),width=0.3) +
    geom_jitter(aes(fill=method_value),stroke = 0, shape = 21,size=5,alpha=0.5,width = 0.05,height = 0,) +
    stat_compare_means()+
    theme_bw() + xlab('') + ylab('Score (median,cov)') + guides(fill=FALSE,color=FALSE)

ggsave("./publication/figures/figure3/sc2_infer_stat_vs_other.pdf",width = 4,height = 4)

full_join(robust_score,method_res, by="teams") %>%
    group_by(method_info,method_value) %>%
    filter(method_info %in% c("Resampling","Infer statistics")) %>%
    spread(method_info,method_value) %>%
    mutate(method = ifelse(`Infer statistics` == "y","statistics", ifelse(Resampling == "y","resamples","other"))) %>%
    select(1:4,method) %>%
    ggplot(aes(method,score_med)) + 
    geom_boxplot(aes(col=method),width=0.3) +
    geom_jitter(aes(fill=method),stroke = 0, shape = 21,size=5,alpha=0.5,width = 0.05,height = 0,) +
    stat_compare_means()+
    theme_bw() + xlab('') + ylab('Score (median,cov)') + guides(fill=FALSE,color=FALSE)


full_join(robust_score,method_res, by="teams") %>%
    group_by(method_info,method_value) %>%
    mutate(method_value = ifelse(method_value=="y","infers population statistic","other method")) %>%
    filter(method_info %in% c("Infer statistics") ) %>%
    group_by(method_value) %>% summarise(m_score = mean(score_med))



## Barplot showing the improvement of the RMSE: 
bar_data <- full_join(robust_score,method_res, by="teams") %>%
    ungroup() %>%
    mutate(method_value = ifelse(method_value=="y","used","not used")) %>%
    filter(!method_info %in% c("Single cell data") ) %>% # "Prior knowledge network",
    mutate(improvement = score_med - mean(score_med))

ggplot(bar_data, aes(paste(method_info,method_value),score_med,fill=method_info)) + 
    geom_bar( position = "dodge", stat = "summary", fun = "mean") +
    geom_point() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 60,hjust = 1))+
    xlab('') +
    ylab('Change in score (delta RMSE)') +
    guides(fill=FALSE) + 
    ggpubr::stat_compare_means(data =bar_data,
                               comparisons = list(c("Cell line metadata used","Cell line metadata not used"),
                                                  c("Cell line similarity used","Cell line similarity not used"),
                                                  c("Custom models used","Custom models not used"),
                                                  #c("Ensemble of models used","Ensemble of models not used"),
                                                  c("Infer statistics used", "Infer statistics not used"),
                                                  c("Linear model used","Linear model not used"),
                                                  #c("Prior knowledge network used","Prior knowledge network not used"),
                                                  c("Resampling used","Resampling not used" ),
                                                  c("Neural network used","Neural network not used"),
                                                  c("Preprocessing used","Preprocessing not used"),
                                                  c("Tree based model used","Tree based model not used")),
                               method = "wilcox.test",label.y = 100,na.rm = TRUE) 









ggplot(bar_data, aes(paste(method_info,method_value),score_med,fill=method_info)) + 
    geom_bar( position = "dodge", stat = "summary", fun = "mean") +
    geom_point() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 60,hjust = 1))+
    xlab('') +
    ylab('Change in score (delta RMSE)') +
    guides(fill=FALSE) + 
    ggpubr::stat_compare_means(data =bar_data,
                               comparisons = list(c("Cell line metadata used","Cell line metadata not used"),
                                                  c("Cell line similarity used","Cell line similarity not used"),
                                                  c("Custom models used","Custom models not used"),
                                                  #c("Ensemble of models used","Ensemble of models not used"),
                                                  c("Infer statistics used", "Infer statistics not used"),
                                                  c("Linear model used","Linear model not used"),
                                                  #c("Prior knowledge network used","Prior knowledge network not used"),
                                                  c("Resampling used","Resampling not used" ),
                                                  c("Neural network used","Neural network not used"),
                                                  c("Preprocessing used","Preprocessing not used"),
                                                  c("Tree based model used","Tree based model not used")),
                               method = "wilcox.test",label.y = 100) 



full_join(robust_score,method_res, by="teams") %>% 
    arrange(score_med)  %>%
    filter(method_info =="Prior knowledge network" ) %>%
    mutate(teams = factor(teams,levels = unique(teams) )) %>% 
    ggplot() + 
    geom_point(aes(teams,score_med,col = method_value)) +
    theme(axis.text.x = element_text(angle = 90,hjust = 1))

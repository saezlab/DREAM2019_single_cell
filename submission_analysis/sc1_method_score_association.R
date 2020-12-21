# for SC1 we check how the score depends on the used methods by t-test. 


library(tidyverse)
library(cowplot)
library(broom)
library(ggpubr)
### A) Leaderboard overview ------------------------------------------------------

submission_folder <- "./submission_data/final_round/SC1"
bootstrap_RMSE <- read_rds("./submission_analysis/intermediate_data/sc1_bootstrap_rmse.rds")
RMSE <- read_rds("./submission_analysis/intermediate_data/sc1_rmse_conditions.rds")
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc1_ranked_teams.rds")

method_usage <- readxl::read_xlsx("./methods/sc1_methods_processed.xlsx",sheet = 2)
descriptor_order <- method_usage$Question

method_usage <-  method_usage %>%
    gather(submitterId,answer,-Question) %>% spread(Question,answer)

# we remove the Anand_1812 from this analysis: 
method_usage <- method_usage[-which(method_usage$submitterId =="anand_1812/CSBL"),]



# fix naming
ranked_teams <- factor(gsub("X.","",ranked_teams,fixed = T),levels = gsub("X.","",ranked_teams,fixed = T))
names(bootstrap_RMSE) <- gsub("X.","",names(bootstrap_RMSE),fixed = T)
names(RMSE) <- gsub("X.","",names(RMSE),fixed = T)



my_colors  = RColorBrewer::brewer.pal(8,"Dark2")



robust_score <- bootstrap_RMSE %>% 
    gather(teams,RMSE,-BS_sample) %>% 
#%>% filter(BS_sample %in% sample(1:1000,10,replace = F))

    group_by(teams) %>% summarise(score_min = min(RMSE),
                                  score_max = max(RMSE),
                                  score_med = median(RMSE))


method_res <- method_usage %>% select(-5,-11,-12, -13) %>%
    mutate(teams = factor(submitterId,levels = levels(ranked_teams))) %>%
    select(-submitterId) %>%
    gather(method_info,method_value,-teams) %>%
    mutate(method_info = factor(method_info,levels = rev(descriptor_order))) 


# T-test:
# H0: using a specific method is as good as other method
# Halt: not using a specific method gives worse (higher) RMSE
full_join(robust_score,method_res, by="teams") %>%
    group_by(method_info) %>%
    #filter(teams!="Huiyuan") %>% # remove the team worse than random
    filter(!method_info %in% c("Prior knowledge network","Single cell data") ) %>%
    nest() %>%
    mutate(t_res = map(data,function(df){
        
        #tres <- t.test(score_med~method_value,data = df,alternative = "less")
        tres <- t.test(score_med~method_value,data = df,alternative = "greater")
        #tres <- t.test(score_med~method_value,data = df,alternative = "two.sided")
        tidy(tres)
    } )) %>% unnest(t_res) %>% arrange(p.value)



full_join(robust_score,method_res, by="teams") %>%
    group_by(method_info,method_value) %>%
    mutate(method_value = ifelse(method_value=="y","used","not used")) %>%
    filter(!method_info %in% c("Prior knowledge network","Single cell data") ) %>%
    summarise(mean_score = mean(score_med),
              min_score = min(score_med),
              max_score = max(score_med)) %>%
    ggplot(aes(method_value,mean_score,fill=method_value)) + 
    geom_bar(stat = "identity",col="black") +
    
    geom_errorbar(aes(ymin= min_score, ymax=max_score)) + 
    facet_wrap(~method_info) + theme_bw() + xlab('') + ylab('Score (mean RMSE)') + guides(fill=FALSE)
ggsave("./publication/figures/figure2/sc1_methods_score_improvement.pdf")

full_join(robust_score,method_res, by="teams") %>%
    group_by(method_info,method_value) %>%
    mutate(method_value = ifelse(method_value=="y","used","not used")) %>%
    filter(!method_info %in% c("Prior knowledge network","Single cell data") ) %>%
    ggbarplot(x = "method_value", y = "score_med", add = "mean_ci",facet.by="method_info")+
    stat_compare_means(label.y = 1.5,vjust = 1,paired = TRUE) 



full_join(robust_score,method_res, by="teams") %>%
    group_by(method_info,method_value) %>%
    mutate(method_value = ifelse(method_value=="y","used","not used")) %>%
    filter(!method_info %in% c("Prior knowledge network","Single cell data") ) %>%
    summarise(mean_score = mean(score_med),
              min_score = min(score_med),
              max_score = max(score_med)) %>%
    ggplot(aes(paste(method_info,method_value),mean_score-mean(mean_score),fill=method_info)) + 
    geom_bar(stat = "identity",col="black",position = position_dodge()) +
#    geom_errorbar(aes(ymin= min_score-mean(mean_score), ymax=max_score-mean(mean_score))) + 
     theme_bw() + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1))+
    xlab('') +
    ylab('Score (mean RMSE)') + guides(fill=FALSE)

ggsave("./publication/figures/figure2/sc1_methods_score_improvement_baseline_based.pdf", width = 6,height = 3)

bar_data <- full_join(robust_score,method_res, by="teams") %>%
    ungroup() %>%
    mutate(method_value = ifelse(method_value=="y","used","not used")) %>%
    filter(!method_info %in% c("Prior knowledge network","Single cell data") ) %>%
    mutate(method_info = ifelse(method_info == "Enseble of models","Ensemble of models",as.character(method_info)))


ranked_methods <- bar_data %>% group_by(method_info) %>%
    filter(!method_info %in% c("Prior knowledge network","Single cell data") ) %>%
    nest() %>%
    mutate(t_res = map(data,function(df){
        
        #tres <- t.test(score_med~method_value,data = df,alternative = "less")
        tres <- t.test(score_med~method_value,data = df,alternative = "greater")
        #tres <- t.test(score_med~method_value,data = df,alternative = "two.sided")
        tidy(tres)
    } )) %>% unnest(t_res) %>% arrange(p.value) %>% pull(method_info)

methods_yn_ordered <- paste(rep(ranked_methods,each=2),rep(c("used","not used"),length(ranked_methods)))


bar_data <- bar_data %>% mutate(method = factor(paste(method_info,method_value),levels = methods_yn_ordered))
#ggplot(bar_data, aes(paste(method_info,method_value),score_med,fill=method_info)) + 
ggplot(bar_data, aes(method,score_med,fill=method_info)) + 
    geom_bar( position = "dodge", stat = "summary", fun = "mean") +
    geom_point() +
    #geom_errorbar(aes(ymin= min_score, ymax=max_score)) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 60,hjust = 1))+
    xlab('') +
    ylab('Score (mean RMSE)') + guides(fill=FALSE) + 
    ggpubr::stat_compare_means(data =bar_data,
                               comparisons = list(c("Cell line metadata used","Cell line metadata not used"),
                                                  c("Cell line similarity used","Cell line similarity not used"),
                                                  c("Ensemble of models used","Ensemble of models not used"),
                                                  c("Linear model used","Linear model not used"),
                                                  c("Models by conditions used","Models by conditions not used"),
                                                  c("Neural network used","Neural network not used"),
                                                  c("Preprocessing used","Preprocessing not used"),
                                                  c("Tree based model used","Tree based model not used")),
                               method = "t.test",method.args = list(alternative="less"),label.y = 2 ) +
    ylim(0,3.2)

ggsave("./publication/figures/figure2/sc1_methods_score_improvement_pvalue.pdf", width = 6,height = 4)




## Barplot showing the improvement of the RMSE: 
bar_data <- full_join(robust_score,method_res, by="teams") %>%
    #filter(teams!="Huiyuan") %>% # remove the team worse than random
    ungroup() %>%
    mutate(method_value = ifelse(method_value=="y","used","not used")) %>%
    filter(!method_info %in% c("Prior knowledge network","Single cell data") ) %>%
    mutate(method_info = ifelse(method_info == "Enseble of models","Ensemble of models",as.character(method_info))) %>%
    mutate(improvement = score_med - mean(score_med))

ggplot(bar_data, aes(paste(method_info,method_value),improvement,fill=method_info)) + 
    geom_bar( position = position_dodge(width=0.5), stat = "summary", fun = "mean") +
    #geom_boxplot( position = "dodge") +
    geom_point() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 60,hjust = 1))+
    xlab('') +
    ylab('Change in score (delta RMSE)') +
    guides(fill=FALSE) + 
    ggpubr::stat_compare_means(data =bar_data,
                               comparisons = list(c("Cell line metadata used","Cell line metadata not used"),
                                                  c("Cell line similarity used","Cell line similarity not used"),
                                                  c("Ensemble of models used","Ensemble of models not used"),
                                                  c("Linear model used","Linear model not used"),
                                                  c("Models by conditions used","Models by conditions not used"),
                                                  c("Neural network used","Neural network not used"),
                                                  c("Preprocessing used","Preprocessing not used"),
                                                  c("Tree based model used","Tree based model not used")),
                               method = "t.test",label.y = .25,method.args = list(alternative = "less")) 
    ylim(-0.1,0.1)

ggsave("./publication/figures/figure2/sc1_methods_score_improvement_pvalue_deltaRMSE.pdf", width = 6,height = 4)

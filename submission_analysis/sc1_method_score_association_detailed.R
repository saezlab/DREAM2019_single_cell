library(tidyverse)
library(cowplot)
library(broom)
library(ggpubr)
### A) Leaderboard overview ------------------------------------------------------

submission_folder <- "./submission_data/final_round/SC1"
bootstrap_RMSE <- read_rds("./submission_analysis/intermediate_data/sc1_bootstrap_rmse.rds")
RMSE <- read_rds("./submission_analysis/intermediate_data/sc1_rmse_conditions.rds")
ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc1_ranked_teams.rds")

method_usage <- readxl::read_xlsx("./methods/sc1_methods_processed.xlsx",sheet = 1)
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

#robust_score <- robust_score %>% filter(score_med<1.3)

method_res <- method_usage %>% 
    mutate(teams = factor(submitterId,levels = levels(ranked_teams))) %>%
    select(-submitterId) %>%
    gather(method_info,method_value,-teams) %>%
    mutate(method_info = factor(method_info,levels = rev(descriptor_order))) 

annotate = robust_score %>% select(teams,score_med) %>% column_to_rownames("teams") %>%
    mutate(score_med = cut(score_med,breaks = c(seq(0.84,0.99,by=0.02),2),labels = FALSE))
rownames(annotate) = robust_score$teams
method_res %>% filter(teams %in%robust_score$teams ) %>%
    mutate(method_value = ifelse(is.na(method_value),"n",method_value)) %>%
    mutate(method_value = ifelse(method_value=="y",1,method_value)) %>%
    mutate(method_value = ifelse(method_value=="n",0,method_value)) %>%
    mutate(method_value = ifelse(method_value%in% c(1,0),method_value,2)) %>%
    mutate(method_value = as.numeric(method_value)) %>% 
    spread(method_info,method_value) %>% column_to_rownames("teams") %>%
    pheatmap::pheatmap(annotation_row =annotate )

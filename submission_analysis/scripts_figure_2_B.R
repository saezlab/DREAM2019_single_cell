
#### B) representative predictions and measured values.-------------------------
library(tidyverse)
library(ranger)



# Participating team names
submissions <-  readRDS("./submission_analysis/intermediate_data/sc1_ranked_teams.rds") %>% as.character()
sub_data_all <- readRDS("./submission_analysis/intermediate_data/sc1_all_NP_predictions.rds")  %>%
    mutate_if(is.character, as.factor) 



# some single cell statistics
# the average standard deviation of a marker in a condition
sub_data_all %>% group_by(cell_line, treatment,time, marker) %>%
    summarise(sd = sd(standard),
              mean = mean(standard)) %>% 
    ungroup() %>%
    summarise(mean_sd=mean(sd),
              `1st qu` = quantile(sd,.25),
              median = median(sd),
              `3st qu` = quantile(sd,.75))

# plot: 
bind_rows(sub_data_all %>% filter(cell_line == 'AU565',treatment=="EGF",time ==7,marker == "p.S6") %>%
              mutate(cond = "mean: 0.853")
          #sub_data_all %>% filter(cell_line == 'LY2',treatment=="iPI3K",time ==40,marker == "p.PLCg2")%>%
          #    mutate(cond = "0.25 percentile: 0.619"),
          #sub_data_all %>% filter(cell_line == 'MACLS2',treatment=="iPI3K",time ==40,marker == "p.HER2")%>%
          #    mutate(cond = "0.75 percentile: 0.965")
) %>% 
    select(cellID, standard,icx_bxai,cond) %>% 
    filter(cellID < 200) %>% 
    group_by(cond) %>%
    mutate(cellID = rank(standard,ties.method ="first")) %>% 
    group_by(cond) %>% 
    gather(source,value,-cellID,-cond) %>% 
    ungroup() %>%
    mutate(source = ifelse(source=="standard","measured","predicted")) %>%
    mutate(cond = factor(cond)) %>% 
    mutate(cond = factor(cond,levels = levels(cond)[c(1,3,2)])) %>%
    ggplot(aes(cellID,value)) +
    geom_line(aes(group=cellID),color="grey90") +
    geom_point(aes(col=source),size=1) +
    ylab("marker intensity") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# + facet_wrap(~cond)
ggsave("./publication/figures/figure2/representative_distribution.pdf",width = 4,height = 3)


bind_rows(sub_data_all %>% filter(cell_line == 'MACLS2',treatment=="EGF",time ==40,marker == "p.PLCg2") %>%
              mutate(cond = "mean: 0.910")
) %>% 
    select(cellID, standard,icx_bxai,cond) %>% 
    filter(cellID < 200) %>% 
    group_by(cond) %>%
    mutate(cellID = rank(standard,ties.method ="first")) %>% 
    group_by(cond) %>% 
    gather(source,value,-cellID,-cond) %>% 
    ungroup() %>%
    mutate(source = ifelse(source=="standard","measured","predicted")) %>%
    mutate(cond = factor(cond)) %>% 
    mutate(cond = factor(cond,levels = levels(cond)[c(1,3,2)])) %>%
    ggplot(aes(cellID,value)) +
    geom_line(aes(group=cellID),color="grey90") +
    geom_point(aes(col=source),size=1) +
    ylab("marker intensity") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# + facet_wrap(~cond)
ggsave("./publication/figures/figure2/representative_distribution_2.pdf",width = 4,height = 3)







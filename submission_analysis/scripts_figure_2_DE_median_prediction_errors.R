
library(tidyverse)
library(cowplot)

median_data <- read_rds("./submission_analysis/intermediate_data/sc1_median_conditions.rds")

### MEK independent ERK signaling 
# we compare 2 cell-lines ("AU565","HCC2218") 
plt_data <- median_data %>%
    #filter(cell_line %in% c("AU565","HCC2218")) %>%
    filter(marker=="p.ERK") %>%
    filter(treatment %in% c("iMEK","EGF")) %>%
    gather(team,prediction,6:27)


# Marco suggestion
plt_data %>% ggplot() +
    geom_line(aes(time,prediction,group=paste(team,treatment), color= treatment, size="predicted"), alpha = 0.2) +
    geom_point(aes(time,standard,col=treatment)) +
    geom_line(aes(time,standard,col=treatment,size="measured")) +
    facet_grid(~cell_line) + theme_bw() + 
    ylab("p-ERK intensity") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank()) +
    scale_size_manual(values = c(predicted=0.5, measured = 1))
ggsave("./publication/figures/figure2/sc1_mek_indep_ERK.pdf",width = 8,height = 2.5)


### EGF treatment + inhibitor
median_data %>% filter(treatment %in% c("EGF","iMEK"), marker=="p.ERK") %>%
    gather(team,prediction,6:27) %>%
    ggplot() +
    geom_line(aes(time,prediction,group=paste(team,treatment), color= marker, size="predicted"), alpha = 0.2) +
    geom_point(aes(time,standard,col=marker)) +
    geom_line(aes(time,standard,col=marker,group=paste(team,treatment), size="measured")) +
    facet_grid(cell_line~marker) + theme_bw() + 
    ylab("Marker intensity") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank()) +
    scale_size_manual(values = c(predicted=0.5, measured = 1))


### iEGFR independent AKT signaling  ----------
# we compare 2 cell-lines ("AU565","LY2") 
plt_data3 <- median_data %>%
    filter(cell_line %in% c("AU565","LY2")) %>%
    filter(marker=="p.Akt.Ser473.") %>%
    filter(treatment %in% c("iEGFR","EGF")) %>%
    gather(team,prediction,6:27)


# Marco suggestion
plt_data3 %>% ggplot() +
    geom_line(aes(time,prediction,group=paste(team,treatment), color= treatment, size="predicted"), alpha = 0.2) +
    geom_point(aes(time,standard,col=treatment)) +
    geom_line(aes(time,standard,col=treatment,size="measured")) +
    facet_grid(~cell_line) + theme_bw() + 
    ylab(bquote('p-'*~Akt^Ser473*' intensity' ) ) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank()) +
    scale_size_manual(values = c(predicted=0.5, measured = 1))

ggsave("./publication/figures/figure2/sc1_egfr_indep_Akt.pdf",width = 5,height = 2.5)



### all prediction in one plot 

gg_data <- median_data %>%
    filter(treatment != "full") %>%
    gather(team,prediction,6:27) %>%
    mutate(marker = ifelse(marker == "p.Akt.Ser473.","pAkt.S473",marker))
    
gg_data %>%
    filter(marker %in% c("p.ERK","p.HER2")) %>%
    ggplot() +
    geom_line(aes(time,prediction,group=paste(team,treatment), color=marker, size="predicted"), alpha = 0.2) +
    geom_point(aes(time,standard,col=marker)) +
    geom_line(aes(time,standard,col=marker,size="measured")) +
    facet_grid(marker+treatment~cell_line,scales = "free_y") + theme_bw() + 
    ylab("Signaling intensity") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank()) +
    scale_size_manual(values = c(predicted=0.5, measured = 1))

ggsave("./publication/figures/figure2/median_prediction_errors_p1.pdf",width = 8,height = 6) 

gg_data %>%
    filter(!marker %in% c("p.ERK","p.HER2")) %>%
    ggplot() +
    geom_line(aes(time,prediction,group=paste(team,treatment), color=marker, size="predicted"), alpha = 0.2) +
    geom_point(aes(time,standard,col=marker)) +
    geom_line(aes(time,standard,col=marker,size="measured")) +
    facet_grid(marker+treatment~cell_line,scales = "free_y") + theme_bw() + 
    ylab("Signaling intensity") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank()) +
    scale_size_manual(values = c(predicted=0.5, measured = 1))

ggsave("./publication/figures/figure2/median_prediction_errors_p2.pdf",width = 8,height = 9) 



#### p.HER 2 in HCC2218 looks very bad:

plt_data_pHER2_HCC2218 <- median_data %>%
    #filter(cell_line %in% c("HCC2218")) %>%
    filter(marker=="p.HER2") %>%
    #filter(treatment %in% c("iEGFR","EGF")) %>%
    filter(treatment %in% c("EGF")) %>%
    gather(team,prediction,6:27)


plt_data_pHER2_HCC2218 %>% ggplot() +
    geom_line(aes(time,prediction,group=paste(team,treatment), color= treatment, size="predicted"), alpha = 0.2) +
    geom_point(aes(time,standard,col=treatment)) +
    geom_line(aes(time,standard,col=treatment,size="measured")) +
    facet_grid(~cell_line) + theme_bw() + 
    ylab('p-HER2' ) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank()) +
    scale_size_manual(values = c(predicted=0.5, measured = 1))

ggsave("./publication/figures/figure2/sc1_HER2_acrosss_cell_lines_EGF.pdf",width = 8,height = 2.5)


#### p.S6 in all look bad

plt_data_pS6 <- median_data %>%
    filter(cell_line %in% c("HCC2218")) %>%
    filter(marker=="p.S6") %>%
    #filter(treatment %in% c("iEGFR","EGF")) %>%
    filter(treatment %in% c("EGF")) %>%
    gather(team,prediction,6:27)


plt_data_pS6 %>% ggplot() +
    geom_line(aes(time,prediction,group=paste(team,treatment), color= treatment, size="predicted"), alpha = 0.2) +
    geom_point(aes(time,standard,col=treatment)) +
    geom_line(aes(time,standard,col=treatment,size="measured")) +
    facet_grid(~cell_line) + theme_bw() + 
    ylab('p-HER2' ) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank()) +
    scale_size_manual(values = c(predicted=0.5, measured = 1))


# EFM19 p.AKT iEGFR: seems high disagreement between participants
plt_data_AKT1 <- median_data %>%
    filter(marker %in% c("p.S6","p.Akt.Ser473."), treatment=="iEGFR",cell_line=="EFM19") %>%
    gather(team,prediction,6:27)


plt_data_AKT1 %>% ggplot() +
    geom_line(aes(time,prediction,group=paste(team,treatment), color= treatment, size="predicted"), alpha = 0.2) +
    geom_point(aes(time,standard,col=treatment)) +
    geom_line(aes(time,standard,col=treatment,size="measured")) +
    facet_grid(~marker) + theme_bw() + 
    ylab('p.Akt.Ser473.' ) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank()) +
    scale_size_manual(values = c(predicted=0.5, measured = 1))



#### PLCg2 in PI3K : large variance across teams

plt_data4 <- median_data %>% 
    filter(cell_line %in% c("HCC2218")) %>%
    filter(marker %in% c("p.PLCg2")) %>%
    filter(treatment %in% c("iPKC","iMEK")) %>%
    gather(team,prediction,6:27) 


plt_data4 %>% ggplot() +
    geom_line(aes(time,prediction,group=paste(team,treatment), color= paste0(treatment,", predicted")), alpha = 0.2, size = 0.5) +
    geom_point(aes(time,standard,col=paste0(treatment,", measured"))) +
    geom_line(aes(time,standard,col=paste0(treatment,", measured")),size=1) +
    facet_grid(treatment~marker) + theme_bw() + 
    ylab('signal intensity') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank())




plt_data4 <- median_data %>% 
    filter(cell_line %in% c("MDAMB436")) %>%
    filter(marker %in% c("p.Akt.Ser473.")) %>%
    filter(treatment %in% c("iPKC","iPI3K")) %>%
    gather(team,prediction,6:27) 


plt_data4 %>% ggplot() +
    geom_line(aes(time,prediction,group=paste(team,treatment), color= paste0(treatment,", predicted")), alpha = 0.2, size = 0.5) +
    geom_point(aes(time,standard,col=paste0(treatment,", measured"))) +
    geom_line(aes(time,standard,col=paste0(treatment,", measured")),size=1) +
    facet_grid(~treatment) + theme_bw() + 
    ylab('p-S6 intensity') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank())





plt_data4 <- median_data %>% 
    filter(cell_line %in% c("MDAMB436")) %>%
    filter(marker %in% c("p.S6","p.PLCg2")) %>%
    filter(treatment %in% c("iEGFR")) %>%
    gather(team,prediction,6:27) 


plt_data4 %>% ggplot() +
    geom_line(aes(time,prediction,group=paste(team,treatment), color= paste0(treatment,", predicted")), alpha = 0.2, size = 0.5) +
    geom_point(aes(time,standard,col=paste0(treatment,", measured"))) +
    geom_line(aes(time,standard,col=paste0(treatment,", measured")),size=1) +
    facet_grid(marker~treatment) + theme_bw() + 
    ylab('p-S6 intensity') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank())


####  Cluster participatns by solution
names(median_data)[6:27] =paste(names(median_data)[6:27],1:22) 
median_data %>% select(-1:-5) %>% pheatmap::pheatmap(cluster_rows = FALSE,show_rownames = FALSE)

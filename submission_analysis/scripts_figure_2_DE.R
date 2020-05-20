
library(tidyverse)
library(cowplot)

median_data <- read_rds("./submission_analysis/intermediate_data/sc1_median_conditions.rds")

### MEK independent ERK signaling 
# we compare 2 cell-lines ("AU565","HCC2218") 
plt_data <- median_data %>%
    filter(cell_line %in% c("AU565","HCC2218")) %>%
    filter(marker=="p.ERK") %>%
    filter(treatment %in% c("iMEK","EGF")) %>%
    gather(team,prediction,6:27)
# #mutate(time= as_factor(time)) %>% 
# plt_data %>% ggplot() +
#     geom_boxplot(data = filter(plt_data,treatment %in% c("iMEK")), aes(time,prediction,group=time,fill="iMEK prediction"),
#                  color=RColorBrewer::brewer.pal(8,"Dark2")[[2]],outlier.size = 0.5) +
#     geom_line(data=filter(plt_data,treatment %in% c("iMEK"),team=="icx_bxai"), aes(time,prediction,col="iMEK prediction")) +
#     geom_point(data=filter(plt_data,treatment %in% c("iMEK")), aes(time,standard,col=treatment)) +
#     geom_line(aes(time,standard,col=treatment)) +
#     facet_grid(~cell_line) + theme_bw() + 
#     scale_color_manual(name = "", values = c("EGF" = "grey60", 
#                                              "iMEK"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],
#                                              "iMEK prediction"=RColorBrewer::brewer.pal(8,"Dark2")[[2]]), ) +
#     scale_fill_manual(name = "",values = c("iMEK prediction"=RColorBrewer::brewer.pal(8,"Dark2")[[2]])) +
#     ylab("p-ERK intensity") +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# ggsave("./publication/figures/figure2/sc1_mek_indep_ERK.pdf",width = 5,height = 2.5)
# 
# # alternative
# plt_data %>% ggplot() +
#     geom_boxplot(data = plt_data, aes(time,prediction,group=paste(time,treatment),fill=treatment, col=treatment),
#                  outlier.size = 0.5) +
#     #geom_line(data=filter(plt_data,treatment %in% c("iMEK"),team=="icx_bxai"), aes(time,prediction,col="iMEK prediction")) +
#     geom_point(data=filter(plt_data,treatment %in% c("iMEK")), aes(time,standard,col=treatment)) +
#     geom_line(aes(time,standard,col=treatment)) +
#     facet_grid(~cell_line) + theme_bw() + 
#     scale_color_manual(name = "", values = c("EGF" = "grey60", 
#                                              "iMEK"=RColorBrewer::brewer.pal(8,"Dark2")[[1]])) +
#     scale_fill_manual(name = "",values = c("iMEK"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],
#                                            "EGF" = "grey")) +
#     ylab("p-ERK intensity") +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 


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
ggsave("./publication/figures/figure2/sc1_mek_indep_ERK.pdf",width = 5,height = 2.5)



### iEGFR independent AKT signaling  ----------
# we compare 2 cell-lines ("AU565","LY2") 
plt_data3 <- median_data %>%
    filter(cell_line %in% c("AU565","LY2")) %>%
    filter(marker=="p.Akt.Ser473.") %>%
    filter(treatment %in% c("iEGFR","EGF")) %>%
    gather(team,prediction,6:27)
# #mutate(time= as_factor(time)) %>% 
# plt_data3 %>% ggplot() +
#     geom_boxplot(data = filter(plt_data3,treatment %in% c("iEGFR")),
#                  aes(time,prediction,group=time,fill="iEGFR prediction"),
#                  color=RColorBrewer::brewer.pal(8,"Dark2")[[2]],
#                  outlier.size = 0.5) +
#     geom_line(data=filter(plt_data3,treatment %in% c("iEGFR"),team=="icx_bxai"),
#               aes(time,prediction,col="iEGFR prediction")) +
#     geom_point(data=filter(plt_data3,treatment %in% c("iEGFR")),
#                aes(time,standard,col=treatment)) +
#     geom_line(aes(time,standard,col=treatment)) +
#     facet_grid(~cell_line) + 
#     theme_bw() + 
#     scale_color_manual(values = c("EGF" = "grey60",
#                                   "iEGFR"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],
#                                   "iEGFR prediction"=RColorBrewer::brewer.pal(8,"Dark2")[[2]])) +
#     scale_fill_manual(name = "",values = c("iEGFR prediction"=RColorBrewer::brewer.pal(8,"Dark2")[[2]])) +
#     ylab(bquote('p-'*~Akt^Ser473*' intensity' ) ) +
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank())
# 
# 
# ggsave("./publication/figures/figure2/sc1_egfr_indep_Akt.pdf",width = 5,height = 2.5)
# 
# plt_data3 %>% ggplot() +
#     geom_boxplot(data =plt_data3,
#                  aes(time,prediction,group=paste(time,treatment),fill=treatment,col=treatment),
#                  outlier.size = 0.5) +
#     #geom_line(data=filter(plt_data3,treatment %in% c("iEGFR"),team=="icx_bxai"),
#     #          aes(time,prediction,col="iEGFR prediction")) +
#     geom_point(data=filter(plt_data3,treatment %in% c("iEGFR")),
#                aes(time,standard,col=treatment)) +
#     geom_line(aes(time,standard,col=treatment)) +
#     facet_grid(~cell_line) + 
#     theme_bw() + 
#     scale_color_manual(values = c("EGF" = "grey60",
#                                   "iEGFR"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],
#                                   "iEGFR prediction"=RColorBrewer::brewer.pal(8,"Dark2")[[2]])) +
#     scale_fill_manual(name = "",
#                       values = c("EGF"= "grey",
#                                      "iEGFR" = RColorBrewer::brewer.pal(8,"Dark2")[[2]])) +
#     ylab(bquote('p-'*~Akt^Ser473*' intensity' ) ) +
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank())
# 


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






#### p.HER 2 in HCC2218 looks very bad:

plt_data_pHER2_HCC2218 <- median_data %>%
    #filter(cell_line %in% c("HCC2218")) %>%
    filter(marker=="p.HER2") %>%
    filter(treatment %in% c("iEGFR","EGF")) %>%
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

ggsave("./publication/figures/figure2/sc1_HER2_acrosss_cell_lines.pdf",width = 8,height = 2.5)


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


# DREAM challenge post analysis
# biological insight
# SC2: compare the the medians in different conditions
# A. Gabor


library(tidyverse)
# SC2

combined_statistics <- read_rds("./submission_analysis/intermediate_data/sc2_stats_conditions.rds")
sserr <- read_rds("./submission_analysis/intermediate_data/sc2_stats_sumSquared_conditions.rds")
source("./scoring_scripts/score_sc2.R")


required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 
#' read_to_stats
#' 
#' reads data and compute stats --> we only work with the stats afterwards
read_to_stats <- function(file_name){
    read_csv(file_name) %>% 
        select(required_columns) %>%
        data_to_stats(.)
}

# read in the data to get access to the median of the EGF condition (not predicteD)
gs_stats <- tibble(files = list.files("./challenge_data/single_cell_phospho/subchallenge_2/",full.names = TRUE)) %>%
    mutate(data = map(files,read_to_stats))

gs_stats <- gs_stats %>% unnest(data) %>% select(-files)

gs_stats <- gs_stats %>% rename(standard=stat_value)



combined_statistics <- combined_statistics %>% gather(teams,prediction,7:22) %>%
    bind_rows(gs_stats)




### MEK independent ERK signaling ----------------------------------------------


combined_statistics %>% 
    filter(stat_variable=="mean_p.ERK") %>%
    filter(treatment %in% c("EGF","iMEK")) %>% 
    filter(cell_line %in% c("184B5","HCC202","ZR751")) %>%
    ggplot() +
    geom_line(aes(time,prediction,group=paste(teams,treatment), color= treatment, size="predicted"), alpha = 0.4) +
    geom_point(aes(time,standard,col=treatment)) +
    geom_line(aes(time,standard,col=treatment,size="measured")) +
    facet_grid(~cell_line) + theme_bw() + 
    ylab("p.ERK" ) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank()) +
    scale_size_manual(values = c(predicted=0.3, measured = 1))

ggsave("./publication/figures/figure3/sc2_mek_indep_ERK.pdf",width = 5,height = 2.5)




### EGFR indep AKT signaling ----------------------------------------------


combined_statistics %>% 
    filter(stat_variable=="mean_p.Akt.Ser473.") %>%
    filter(treatment %in% c("EGF","iEGFR")) %>% 
    filter(cell_line %in% c("BT483","MCF12A","MDAMB468")) %>%
    ggplot() +
    geom_line(aes(time,prediction,group=paste(teams,treatment), color= treatment, size="predicted"), alpha = 0.4) +
    geom_point(aes(time,standard,col=treatment)) +
    geom_line(aes(time,standard,col=treatment,size="measured")) +
    facet_grid(~cell_line) + theme_bw() + 
    ylab(bquote('p-'*~Akt^Ser473*' intensity' ) ) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),legend.title = element_blank()) +
    scale_size_manual(values = c(predicted=0.3, measured = 1))

ggsave("./publication/figures/figure3/sc2_egfr_indep_akt.pdf",width = 5,height = 2.5)




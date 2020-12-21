## DREAM challenge post analysis
# 
# SC3: comapre the team's error vs EGF
# do this for both SC2 and SC3
# A. Gabor


# - read in the single cell data used in SC2 and SC3 and compute the stats
# - read in the stats computed for the test: both for the data and the predictions
# - compare them 

library(tidyverse)

source("./scoring_scripts/score_sc2.R")

#### 1. compute the data stats for the EGF conditions --------------------------
# part of complete cell-lines are used in SC3:
SC3_cell_lines <- c( 'BT20', 'BT474', 'BT549', 'CAL148', 'CAL51', 'CAL851',
                     'DU4475', 'EFM192A', 'EVSAT', 'HBL100', 'HCC1187',
                     'HCC1395', 'HCC1419', 'HCC1500', 'HCC1569', 'HCC1599',
                     'HCC1937', 'HCC1954', 'HCC2185', 'HCC3153', 'HCC38',
                     'HCC70', 'HDQP1', 'JIMT1', 'MCF7', 'MDAMB134VI', 'MDAMB157',
                     'MDAMB175VII', 'MDAMB361', 'MDAMB415', 'MDAMB453',
                     'MFM223', 'MPE600', 'MX1', 'OCUBM', 'T47D', 'UACC812',
                     'UACC893', 'ZR7530')

complete_CL <- tibble(path = list.files("challenge_data/single_cell_phospho/complete_cell_lines/",
                                        full.names = T)) %>%
    mutate(fbasename = basename(path),
           cell_line = gsub(".csv","",fbasename))

# check that all appears among the complete cell-lines: 
stopifnot(length(SC3_cell_lines[!SC3_cell_lines %in% complete_CL$cell_line])==0)

# import the cell-lines 1-by-1, filter for EGF treatment only and compute the stats

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
    
    read_csv(file_name) %>% filter(treatment =="EGF") %>%
        select(required_columns) %>% data_to_stats(.)
}


sc3_cl_stats <- complete_CL %>%
    filter(cell_line %in% SC3_cell_lines) %>% # filter out the cell-lines not used
    #slice(1:1) %>% 
    mutate(data_stats = map(path,read_to_stats)) 
    
sc3_cl_stats <- sc3_cl_stats %>% select(data_stats) %>%
    unnest(data_stats) 



###############################################################################
### 2. Compare the (1) EGF stats (2) predictions to the test stats ------------

sc3_combined_statistics <- read_rds("./submission_analysis/intermediate_data/sc3_stats_conditions.rds") 

# there are some timepoint missmatches that we have to fix before merging
# adjust times to :  0   7   9  13  17  40  60
sc3_combined_statistics <- sc3_combined_statistics %>% 
    mutate(time = ifelse(time ==14,13,time ),
           time = ifelse(time ==16,17,time ))

sc3_cl_stats <- sc3_cl_stats%>% 
    mutate(time = ifelse(time %in% c(12,14),13,time),
           time = ifelse(time ==35,40,time ))


# now we can merge:
sc3_all <- sc3_combined_statistics %>%
    left_join(sc3_cl_stats %>% rename(EGF_stat = stat_value),by = c("cell_line", "time", "stat_variable"))


# compute the Euclidean distance of the predictions and EGF data stats from the 
# test data (column: standard)

sc3_data_cols <- c( 'icx_bxai',  'orangeballs',  'AMbeRland', 
                    'X.msinkala',  'X.GaoGao199694',  'NAD',  'PaL',  'hipathia', 
                    'X.pqiu',  'CSBL',  'SanGuo',  'X.Huiyuan',  'KAUST_RSS', 
                    'Raghava_India_SCS',  'EGF_stat')

# compute the sum of squared distance from real stats for each condition
# if we sum up all we get the values from the leaderboard --> sanity check OK
sc3_all %>% 
    group_by(cell_line,time) %>% 	
    summarise_at(sc3_data_cols,~ sum((standard - .)^2)) %>%
    ungroup() %>%
    summarise_at(sc3_data_cols,~ sqrt(sum(.))) %>% unlist()


# Which stats change the most across all cell lines? 
# compute the mean absolute diff of the variables across cell-lines:
stats_per_var <- sc3_all %>% mutate(EGF_diff = abs(standard - EGF_stat)) %>% 
    select( cond_id, cell_line, treatment.x,  time, stat_variable, EGF_stat, EGF_diff, everything()) %>%
    group_by(stat_variable) %>%
    summarise_at(c("EGF_stat","standard","EGF_diff",sc3_data_cols), ~mean(.)) %>%
    arrange(desc(EGF_diff))


# select the top 20 most changed statistics between EGF and mTOR
# show the percantage of change: real and predicted by team 

stats_per_var %>% 
    filter(EGF_diff>0.1) %>%
    # top_n(40,EGF_diff) %>%
    gather(teams, prediction, 5:18) %>% 
    select(stat_variable,EGF_stat,standard,teams,prediction ) %>%
    mutate(EGF_stat = (EGF_stat-standard)/standard,
           prediction = (prediction-standard)/standard,
           standard = (standard-standard)/standard) %>%
    # select(-EGF_stat) %>%
    arrange(desc(EGF_stat)) %>% 
    filter( (EGF_stat > .1) | (EGF_stat < -0.1)) %>%
    select(-standard) %>%
    mutate(stat_variable=factor(stat_variable,levels = unique(stat_variable))) %>%
    ggplot() + 
    geom_line(aes(stat_variable,prediction,col=teams,group=teams)) +
    geom_hline(yintercept = 0) +
    theme_bw()  + 
    theme(axis.text.x = element_text(angle = 90,hjust=1),legend.title = element_blank()) +
    xlab("Variable") + 
    ylab("percentage change (EGF -> imTOR)") +
    coord_cartesian(ylim = c(-1,1))

    ###
    ggplot() + 
    geom_col(aes(stat_variable,value*100,fill=source),position="dodge") +
    geom_hline(yintercept = 100) +
    theme_bw()  + 
    theme(axis.text.x = element_text(angle = 90,hjust=1),legend.title = element_blank()) +
    xlab("Variable") + 
    ylab("percentage change (EGF -> imTOR)") +
    coord_cartesian(ylim = c(0,200))
    #ggplot() +geom_point(aes(EGF_stat,standard,col=diff))

ggsave(filename = "publication/figures/figure4/stat_change_egf_vs_team.pdf",height = 5,width = 5)

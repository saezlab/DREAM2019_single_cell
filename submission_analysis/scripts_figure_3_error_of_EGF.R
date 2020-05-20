## DREAM challenge post analysis
# 
# SC2: comapre the team's error, A and B replica vs EGF
# do this for both SC2 and SC3
# A. Gabor


# - read in the single cell data used in SC2 and SC3 and compute the stats
# - read in the stats computed for the test: both for the data and the predictions
# - compare them 

library(tidyverse)

source("./scoring_scripts/score_sc2.R")

#### 1. compute the data stats for the EGF conditions --------------------------
# cell-lines of SC2

SC2_CL <- tibble(path = list.files("challenge_data/single_cell_phospho/subchallenge_2/",
                                        full.names = T)) %>%
    mutate(fbasename = basename(path),
           cell_line = gsub(".csv","",fbasename))


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


sc2_cl_stats <- SC2_CL %>%
    mutate(data_stats = map(path,read_to_stats)) %>%
    select(cell_line,data_stats) %>%
    unnest(data_stats) %>%
    select(-cell_line1)


###############################################################################
### 2. Compare the (1) EGF stats (2) predictions to the test stats ------------

sc2_combined_statistics <- read_rds("./submission_analysis/intermediate_data/sc2_stats_conditions.rds") %>%
    add_column(SC=2)

# there are some timepoint missmatches that we have to fix before merging
# adjust times to :  0   7   9  13  17  40  60

sc2_combined_statistics <- sc2_combined_statistics %>%
    mutate(time = ifelse(time ==16,17,time )) 



sc2_all <- sc2_combined_statistics %>%
    left_join(sc2_cl_stats %>% rename(EGF_stat = stat_value),by = c("cell_line", "time", "stat_variable"))


# compute the Euclidean distance of the predictions and EGF data stats from the 
# test data (column: standard)

sc2_data_cols <- c('EGF_stat',  'icx_bxai',  'X.pqiu',  'orangeballs',  'PaL',  'NAD',  
                     'SCG',  'AMbeRland',  'Sleeping',  'hipathia',  'hulab.SCS',
                     'CSBL',  'X.msinkala',  'X.Huiyuan',  'SanGuo',  
                     'Raghava_India_SCS',  'KAUST_RSS')
# compute the sum of squared distance from real stats for each condition
# if we sum up all we get the values from the leaderboard --> sanity check OK


sc2_all %>% 
    group_by(cell_line,time,treatment.x) %>% 	
    summarise_at(sc2_data_cols,~ sum((standard - .)^2)) %>%
    ungroup() %>%
    summarise_at(sc2_data_cols,~ sqrt(sum(.))) %>%
    select(EGF_stat,everything())

# EGF_stat icx_bxai X.pqiu orangeballs   PaL   NAD   SCG AMbeRland Sleeping hipathia hulab.SCS  CSBL X.msinkala X.Huiyuan SanGuo
# <dbl>    <dbl>  <dbl>       <dbl> <dbl> <dbl> <dbl>     <dbl>    <dbl>    <dbl>     <dbl> <dbl>      <dbl>     <dbl>  <dbl>
# 1     29.2     23.2   25.7        26.6  27.1  30.3  33.5      35.9     58.8     60.3      61.5  66.8       75.9      86.2   90.5
# … with 2 more variables: Raghava_India_SCS <dbl>, KAUST_RSS <dbl>


# check the score for the subsets where Relica A and B overlaps: 
sc2_all %>% 
    filter(time==0) %>%
    group_by(cell_line,time,treatment.x) %>% 	
    summarise_at(sc2_data_cols,~ sum((standard - .)^2)) %>%
    ungroup() %>%
    summarise_at(sc2_data_cols,~ sqrt(sum(.))) %>%
    select(EGF_stat,everything())




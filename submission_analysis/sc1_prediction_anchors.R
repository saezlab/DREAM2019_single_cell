# naive and good models

library(tidyverse)
# SC1 

# load data
all_reporters <- c(	  'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
					  'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
					  'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
					  'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
					  'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
					  'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
					  'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
					  'p.STAT3', 'p.STAT5') 

reporters <- c("p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")
gs <- read_csv("./challenge_data/validation_data/sc1gold.csv") %>%
	arrange(glob_cellID) %>%
	gather(key = "marker", value = "standard", reporters)


rmse <- function(x,y){
	sqrt(sum((x-y)^2)/length(x))
}


### 1: Shuffling by condition---------------------------------------------------
# shuffles the golden standard data by groups and computes the mean RMSE across 
# conditions
iter <- 0
randomised_conditional_score <- function(i,gs){
	iter <<- iter + 1
	print(iter)
	gs %>% group_by(cell_line,treatment, time, marker) %>%
		mutate(random_prediction = sample(standard)) %>% 
		summarise(RMSE_random = rmse(standard,random_prediction)) %>%
		pull(RMSE_random) %>% mean()
		
	
}

random_pred_score <- tibble(samp = seq(1:5)) %>% rowwise() %>% 
	mutate( random_score = map_dbl(samp,randomised_conditional_score,gs))

random_pred_score %>% summarise(mean = mean(random_score),
 								sd = sd(random_score))
# A tibble: 1 x 2
# mean       sd
# <dbl>    <dbl>
# 	1  1.44 0.000262


### 2: Shuffling across all  ---------------------------------------------------
# shuffles the golden standard data overall and computes the mean RMSE across
# conditions
iter <- 0
randomised_score <- function(i,gs){
	iter <<- iter + 1
	print(iter)
	gs %>% mutate(random_prediction = sample(standard)) %>% 
		summarise(RMSE_random = rmse(standard,random_prediction)) %>%
		pull(RMSE_random) %>% mean()
	
	
}

random_pred_score <- tibble(samp = seq(1:5)) %>% rowwise() %>% 
	mutate( random_score = map_dbl(samp,randomised_score,gs))

random_pred_score %>% ungroup() %>%
	summarise(mean = mean(random_score),
			  sd = sd(random_score))
# A tibble: 1 x 2
# mean       sd
# <dbl>    <dbl>
# 	1  2.21 0.000370


### 3: neirest neighbours ------------------------------------------------------
sc1_data = map(list.files("./challenge_data/single_cell_phospho/subchallenge_1/", full.names = T),read_csv) %>%
	bind_rows()
gs_orig <- read_csv("./challenge_data/validation_data/sc1gold.csv")	
	
sc1_total <- sc1_data %>% select(-reporters) %>% arrange(cell_line, treatment, time, cellID, fileID) %>%
	bind_cols(gs_orig %>% arrange(cell_line, treatment, time, cellID, fileID) %>% select(reporters))

temp_data = sc1_total %>% filter(cell_line=="AU565",treatment=="EGF",time==0)


library(caret)
knn_model <- FNN::knn.reg(train = temp_data %>% select(all_reporters,-p.ERK),
						  y =  temp_data %>% pull(p.ERK),
						  test = NULL,
						  k = 5)

tibble(pred = knn_model$pred, data= temp_data$p.ERK) %>% ggplot() + geom_point(aes(pred,data))

# this would take for long

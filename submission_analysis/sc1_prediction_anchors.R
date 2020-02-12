# 
# we pick some metrics on the golden standard, to show how good are the models. 

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
	mutate( random_score = map_dbl(samp,randomised_conditional_score,gs)) %>%
	ungroup()

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
	mutate( random_score = map_dbl(samp,randomised_score,gs)) %>%
	ungroup()


random_pred_score %>% ungroup() %>%
	summarise(mean = mean(random_score),
			  sd = sd(random_score))
# A tibble: 1 x 2
# mean       sd
# <dbl>    <dbl>
# 	1  2.21 0.000370


### 3: predict in-cell-line ----------------------------------------------------
# the idea here is the following: 
# we measure the variance 

ch_data <- list.files("./challenge_data/single_cell_phospho/subchallenge_1/", full.names = T) %>%
	map_dfr(read_csv)

gs_wide <- gs %>% spread(marker,standard)


gs_all <- full_join(ch_data %>% select(-reporters),
					gs_wide, by = c("treatment", "cell_line", "time", "cellID", "fileID"))
library(caret)

compute_RMSE <- function(scdata){
	RMSE <- c()
	for(i in 1:length(reporters)){
		model_data <- scdata %>% select( reporters[[i]],setdiff(all_reporters,reporters))
		
		train_control<- trainControl(method="cv", number=10, savePredictions = TRUE)
		model<- train(reformulate(termlabels =setdiff(all_reporters,reporters),
								  response = reporters[[i]]), data=model_data, trControl=train_control, method="lm")
		RMSE[[i]] <- model$results$RMSE
		
	}
	names(RMSE) <- reporters
	return(enframe(RMSE,name = "marker"))
}


estim_RMSE_cl_tr_tm <- gs_all %>%
	group_by(cell_line, treatment, time) %>%
	nest() %>% 
	mutate(RMSE = map(data,compute_RMSE)) %>% select(-data) %>% unnest(RMSE) 
estim_RMSE_cl_tr_tm %>% pull(value) %>% mean()

estim_RMSE_global <- gs_all %>%
	nest() %>% 
	mutate(RMSE = map(data,compute_RMSE)) %>% select(-data) %>% unnest(RMSE) 

estim_RMSE_tr_tm <- gs_all %>%
	group_by(treatment, time) %>%
	nest() %>% 
	mutate(RMSE = map(data,compute_RMSE)) %>% select(-data) %>% unnest(RMSE) 

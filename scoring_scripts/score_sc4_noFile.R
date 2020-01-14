
library(dplyr)
library(tidyr)
library(readr)


# scores the subchallenge aim 2
#' @param prediction_data prediction data  (tibble)
#' @param validation_data validation data (tibble)
#' @description checks input for missing columns
#' computes root-mean square error by conditions, then averages these


score_sc4 <- function(prediction_data,validation_data){
	
	# HER2, PLCg2, cellID and fileID not included on purpose ! 
	required_columns <- c('cell_line','treatment', 'time',
						  'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
						  'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
						  'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
						  'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
						  'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
						  'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
						  'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
						  'p.STAT3', 'p.STAT5') 
	
	prediction_data = prediction_data %>% select(required_columns)
	# as we agreed, we remove plcg and her2 from the validation data:
	validation_data = validation_data %>% select(required_columns)

		### Formating -------------------------
	# join the test and validation data
	
	combined_data = full_join(prediction_data %>% gather(marker,prediction,-cell_line, -treatment, -time ),
							  validation_data %>% gather(marker,test,-cell_line, -treatment, -time ),
							  by=c("cell_line", "treatment", "time","marker"))
	
	### Calculate score --------------------
	# calculate the  distance over all stats
	RMSE_cond = combined_data %>% group_by(cell_line,treatment,marker) %>% 
		summarise(RMSE = sqrt(sum((test - prediction)^2)/n())) 
	
	final_score = mean(RMSE_cond$RMSE)
}






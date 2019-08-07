
# How about MI instead of covariance ?
# how about weighting stats stronger if it is further from zero?

library(dplyr)
library(tidyr)





# scores the subchallenge aim 1.2.1 and 1.2.2
#' @param prediction_data_file path to prediction data file (.csv)
#' @param validation_data_file path to validation data file (.csv)
#' @description checks input for missing columns
#' check input for missing conditions (missing predicted cells)
#' computes root-mean square error by conditions, then averages these


validate_aim_1_2 <- function(prediction_data_file,validation_data_file){
	# to be returned:
	error_status = list(state=0,message="")
	
	# load validation data
	validation_data <- read_csv (validation_data_file) %>% select(-fileID,-cellID)
	prediction_data <- read_csv(prediction_data_file)
	
	### Checking inputs -------------------
	# checking columns of input data table
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
	
	if(!all(required_columns %in% names(prediction_data))) {
		stop(paste0("missing columns detected. Required columns: ", paste(required_columns,collapse = ", ")))
	}
	prediction_data = prediction_data %>% select(required_columns)
	# as we agreed, we remove plcg and her2 from the validation data:
	validation_data = validation_data %>% select(required_columns)
	
	# checking for any missing conditions
	required_conditions <- validation_data %>% select(cell_line,treatment,time) %>% unique()
	predicted_conditions <- prediction_data %>% select(cell_line,treatment,time)
	
	missing_conditions = anti_join(required_conditions,predicted_conditions,by = c("cell_line", "treatment", "time"))
	
	if(nrow(missing_conditions)>0){
		print(missing_conditions %>% select(c("cell_line", "treatment", "time")))
		stop("missing predictions detected for above conditions")	
	} 
	return(error_status)
}







library(dplyr)
library(tidyr)
library(readr)

#' validate_sc4
#' @param prediction_data_file path to prediction data file (.csv)
#' @param validation_data_file path to validation data file (.csv)
#' @description  validates the predictions for subchallenge 4
#' checking:
#'  - missing columns
#'  - missing conditions
#'  - NAs in the prediction
#'  @return error_status list(), with fields state and message. 
#'  state = 0 if validation was successful
validate_sc4 <- function(prediction_data_file,validation_data_file){
	# to be returned:
	error_status = list(state=0,message="")
	
	# load validation data
	validation_data <- read_csv (validation_data_file) 
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
		error_status$state = -1
		error_status$message = paste0("missing columns detected. Required columns: ", paste(required_columns,collapse = ", "))
		return(error_status)
	}
	
	prediction_data = prediction_data %>% select(required_columns)
	# as we agreed, we remove plcg and her2 from the validation data:
	validation_data = validation_data %>% select(required_columns)
	
	# checking for any missing conditions
	required_conditions <- validation_data %>% select(cell_line,treatment,time) %>% unique()
	predicted_conditions <- prediction_data %>% select(cell_line,treatment,time)  %>% unique()
	
	missing_conditions = anti_join(required_conditions,predicted_conditions,by = c("cell_line", "treatment", "time"))
	
	if(nrow(missing_conditions)>0){
		error_status$state = -2
		error_status$message = 
			paste0("missing predictions detected for the following conditions (cell_line_treatment_time): ",
				   missing_conditions %>%
				   	unite(condition_id) %>%
				   	pull(condition_id) %>%
				   	paste(.,collapse = ", ")
			)
		return(error_status)
	} 
	
	# checking for NAs
	contains_NA_prediction <- prediction_data %>%  is.na() %>% any()
	if(contains_NA_prediction){
		error_status$state = -3
		error_status$message = paste0("NA value detected in the predictions.")
		return(error_status)
		
	}
	
	return(error_status)
}






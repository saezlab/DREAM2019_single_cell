library(readr)
library(tidyr)
library(dplyr)

#' validate_sc1
#' @param prediction_data_file path to prediction data file (.csv)
#' @param validation_data_file path to validation data file (.csv)
#' @description validates the predictions for subchallenge 1
#' checking:
#'  - missing columns
#'  - missing cells
#'  - NAs in the prediction
#'  @return error_status list(), with fields state and message. 
#'  state = 0 if validation was successful
validate_sc1 <- function(prediction_data_file, validation_data_file) {

  # stores the error message and state
  error_status = list(state=0,message="successful validation")
	
  # load validation data
  validation_data <- read_csv(validation_data_file)
  prediction_data <- read_csv(prediction_data_file)

  # Check for delim issues
  if (length(prediction_data)==1){
  	error_status$state = -1
  	error_status$message = paste0("Only 1 column detected. Make sure to submit comma-separated table. ",
  								  "We recommend readr::write_csv() if you are an R user.")
  	return(error_status)
  }
  
  ### Checking inputs -------------------
  # checking columns of input data table
  required_columns <- c("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID", "p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")

  if (!all(required_columns %in% names(prediction_data))) {
  	error_status$state = -1
  	error_status$message = paste0("missing columns detected. Required columns: ", paste(required_columns, collapse = ", "))
  	return(error_status)
  }
  
  prediction_data <- prediction_data %>% select(required_columns)

  # checking for any missing cell-predictions
  missing_prediction_data <- anti_join(
  	validation_data %>% select("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID"),
  	prediction_data %>% select("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID"), by = c("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID"))
  if (nrow(missing_prediction_data) > 0) {
  	
  	error_status$state = -2
  	
  	missing_ids = missing_prediction_data  %>% unite("id", glob_cellID,cell_line, treatment, time, cellID, fileID, sep = "_")
  	
  	error_status$message = paste0("missing predictions for cells (glob_cellID_cell_line_treatment_time_cellID_fileID): ", paste(missing_ids$id,collapse = ", "))
  	return(error_status)
  	
  }
  
  # checking for NAs
  contains_NA_prediction <- prediction_data %>% select("p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6") %>% is.na() %>% any()
  if(contains_NA_prediction){
  	error_status$state = -3
  	error_status$message = paste0("NA value detected in the predictions.")
  	return(error_status)
  	
  }
  
  return(error_status)
  
}

# utilities to process results

#' fix_leaderboard_raw
#' 
#' @param LB leaderboard table, imported from excel
#' @description  fix the problem of copy-pasting from leaderboard to excel: 
#' some names gets into 2 consequitive cells, but the score is only in the first. 
#' also look for NaNs (text) and converts to NA (numeric)

fix_leaderboard_raw <- function(LB){
	rm_rows = c()
	for(i in 2:nrow(LB)){
		
		# running name: name for prev row is in the this row, (other elements are NA)
		NA_name = slice(LB,i) %>% pull(submitterId) %>% is.na()
		is_name_running = !NA_name & (slice(LB,i) %>% select(-submitterId) %>% is.na() %>% all())
		if(is_name_running){
			LB[i-1,"submitterId"] = LB[i,"submitterId"]
			rm_rows = c(rm_rows,i)
		}
	}
	
	LB = LB[-rm_rows,]
	
	# fix NAN
	LB$score <- as.numeric(as.character(LB$score))
	return(LB)
}


#'  fix the dates of the submission
#'
#' dates comes with a free format that should be converted
#' @param date_string  vector containing the dates in the leaderboard
#' @return date as a posixct class object

fix_submission_dates <- function(date_string){
	
	as.POSIXct(date_string, tryFormats=c("%m/%d/%Y %I:%M %p","%m/%d/%Y %H:%M"))
}

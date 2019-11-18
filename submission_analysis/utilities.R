# utilities to process results

#' fix_leaderboard_raw
#' input: 
#'	@param LB leaderboard table, imported from excel
#' @description  fix the problem of copy-pasting from leaderboard to excel: 
#' some names gets into 2 consequitive cells, but the score is only in the first. 
#' also look for NaNs (text) and converts to NA (numeric)

fix_leaderboard_raw <- function(LB){
	rm_rows = c()
	for(i in 2:nrow(LB)){
		
		if(slice(LB,i) %>% select(-submitterId) %>% is.na() %>% all()){
			
			LB[i-1,"submitterId"] = LB[i,"submitterId"]
			rm_rows = c(rm_rows,i)
		}
	}
	
	LB = LB[-rm_rows,]
	
	# fix NAN
	LB$score <- as.numeric(as.character(LB$score))
	return(LB)
}
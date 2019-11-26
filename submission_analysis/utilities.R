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
#' dates comes in various format because the copy-paste excel import messup. 
#' we convert to a unified format:
#' - fix the '19 --> 2019 , then the date and month is also changed... 
#' - handle the AM/PM and the 24 hours base notation
#' @param date_string  vector containing the dates in the leaderboard
#' @return date as a posixct class object

fix_submission_dates <- function(date_string){
	
	# handles the date in a day.month.year AM/PM format, where year is only by 2 digits
	# e.g "9.10.19 6:28 PM"   
	# returns date in POSIXct format
	DMY_fix <- function(str_in){
		
		gsub(x = str_in, pattern = ".19 ",replacement = ".2019 ",fixed = T) %>%
			as.POSIXct(tryFormats=c("%d.%m.%Y %I:%M %p","%m.%d.%Y %I:%M %p"))
	}
	# handles the date in a month.day.year (AM/PM) format
	# e.g. "10.01.2019  07:31:00"   or "08.12.2019  21:48:00 PM"  "09.27.2019 10:20 PM"
	# returns date in POSIXct format
	MDY_convert <- function(str_in){
		
		as.POSIXct(str_in, tryFormats=c("%m.%d.%Y %I:%M %p","%m.%d.%Y %H:%M"))
	}
	
	timetable = tibble(date = date_string, order = 1:length(date_string))
	
	
	tt1 <- timetable %>% 
		filter(grepl(".19 ",date,fixed = TRUE)) %>%
		mutate(fixed_date = DMY_fix(date))
	tt2 <- timetable %>% 
		filter(!grepl(".19 ",date,fixed = TRUE)) %>%
		mutate(fixed_date = MDY_convert(date))
	
	
	results <- bind_rows(tt1,tt2) %>%
		arrange(order) %>%
		pull(fixed_date)
	return(results)
}

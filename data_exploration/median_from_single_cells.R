# Marco's median data does not have the Her2 and plcg, better to reproduce it 
# from single cell data

library(tidyverse)
library(DBI)
library(RSQLite)
library(progress)

con <- dbConnect(RSQLite::SQLite(), "./data/cleaned_single_cell_data/single_cell_dream_cls.sqlite")

tables <-  dbListTables(con)
cols <- dbListFields(con,tables[[1]])
timepoints <- dbGetQuery(con,paste0("SELECT DISTINCT time FROM  ",dbQuoteIdentifier(con,tables[[1]])))
treatments <- dbGetQuery(con,paste0("SELECT DISTINCT treatment FROM  ",dbQuoteIdentifier(con,tables[[1]])))


bar = progress::progress_bar$new(format = "  Processing [:bar] :percent eta: :eta",
								 total = length(tables))

median_data = sapply(tables, function(cl){
	#cl = tables[[1]]
	bar$tick()
	
	# load cell_line and limit the number of cells
	tmp = dbReadTable(con,  dbQuoteIdentifier(con,cl)) %>% 
		as_tibble()
	tmp %>% select(-cellID) %>%
		group_by(cell_line,treatment, time, fileID) %>% 
		summarise_all(median,na.rm=TRUE) %>% ungroup()
	
},simplify = F)


median_data = do.call("rbind",median_data) %>% as_tibble()
write_rds(median_data,"./data/median_data/median_all_reporters_mine.rds")



### Quality check --------------------------------------------------------------
# compare the new and the old median
# there are differences for sure:
# - single cell dataset: no info on timecourse A vs B
# - single cell data was normalised/batch-corrected differently

# I expected differences, because 
# (1) the batch correction was done differently in the 2 cases (on single cell level vs on the median level); 
# (2) the single cell files dont have information about if they belong to time course  A or B and finally 
# (3) I was averaging the median data on technical replicates (over fileID) in the new median, and taking the mean of A/B in the old. 

old_median_raw = read_rds("./data/median_data/Median_allsamples_nocontrols_withcellcount.rds") %>% as_tibble()

old_median <- old_median_raw %>% select(-`dmt$cellcount`) %>%
	gather(reporter,value,-cell_line,-treatment,-time,-time_course) %>%
	mutate(reporter = make.names(reporter)) %>%
	mutate(time = as.character(time),
		   time = ifelse(time=="0short","0",time),
		   time = ifelse(is.na(time),"0",time),
		   time = as.numeric(time)
	) %>%
	group_by(cell_line,treatment,time,reporter) %>% 
	summarise(value = mean(value, na.rm = TRUE) ) %>% # mean over time_course
	ungroup() %>%
	rename(old_median=value)

new_median = median_data %>% 
	gather(reporter,value,-cell_line,-treatment,-time,-fileID) %>%
	group_by(cell_line,treatment,time,reporter) %>% 
	summarise(value = mean(value, na.rm = TRUE) ) %>% # mean over fileID
	ungroup() %>%
	rename(new_median=value)
	
	
	
medians = full_join(old_median,new_median, by = c("cell_line", "treatment", "time", "reporter"))


library(ggplot2)

medians %>% ggplot(aes(old_median,new_median)) + geom_point() + geom_abline(slope = 1,intercept = 0) + 
	coord_equal() + facet_wrap(~reporter)


###### MEDIAN INTERPOLATION -----------------------------------------------------
## Median data interpolation
# the code will skip this section and data is read from the disk. 

# Load the median data to find suitable cell-lines

median_data <- read_rds("./data/median_data/median_all_reporters_mine.rds") %>% as_tibble()

# preparing time-courses:
# 	we interpolate the data to the same timepoints:


# uncomment if we use Marco's data
# median_data <- median_data %>% mutate(time = as.character(time)) %>%
# 	mutate(time = ifelse(time=="0short","0",time)) %>%
# 	mutate(time = ifelse(is.na(time),"0",time)) %>%
# 	mutate(time = as.numeric(as.character(time))) 

median_data%>%pull(time) %>% table()
base_time = c(0,5.5,7,9,13,17,23,30,40,60)

# In case of Marco's median data interpolation in 2 steps:
# 1. interpolate for the basal time per time_course
# - first we filter out treatment full, because that has only time 0 so we dont want to interpolate.
# - for each reporter time_course we interpolate for the base_timepoints, defined above
# 2. interpolate (average) between the time_courses


median_interpolated_data = median_data %>% select(-`dmt$cellcount`) %>% filter(treatment != "full") %>%
	gather(reporter,value,-cell_line,-treatment,-time,-time_course) %>%
	group_by(cell_line,treatment,time_course,reporter) %>% # interpolating per time_course!
	nest() %>%
	mutate(interp_value = map(data, function(data){
		
		out = approx(x = data$time,y = data$value,xout = base_time,method ="linear",rule = 2)
		data.frame(time=out$x,value=out$y)
	}
	))%>% unnest(interp_value) %>%
	# bind the EGF treatment: this has only time 0
	bind_rows( median_data %>% select(-`dmt$cellcount`) %>%
			   	filter(treatment == "full") %>%
			   	gather(reporter,value,-cell_line,-treatment,-time,-time_course)  
	) %>%
	
	group_by(cell_line,treatment,time,reporter) %>% # averaging over time_course!
	summarise(value = mean(value))

if(FALSE) saveRDS(median_interpolated_data,"./data/median_data/interpolated_median_allsamples_nocontrols_withcellcount.rds")

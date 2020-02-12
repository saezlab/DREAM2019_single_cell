# Marco's median data does not have the Her2 and plcg, better to reproduce it 
# from single cell data

# compute the median response for each reporter in each condition and
# in each replica (identified by the fileID)
#
# in the next step we map the fileIDs to time-course A and time-course B
# and compute the weighted mean response over multiple fileIDs if they belong
# to the same condition and time_course (there are some technical replica) 
#
# finally we interpolate for the general time courses 
# 
library(tidyverse)
library(DBI)
library(RSQLite)
library(progress)

con <- dbConnect(RSQLite::SQLite(), "./data/cleaned_single_cell_data/single_cell_dream_cls.sqlite")

cell_lines <-  dbListTables(con)

cell_lines = cell_lines[cell_lines!="HCC70_2"]

bar = progress::progress_bar$new(format = "  Processing [:bar] :percent eta: :eta",
								 total = length(cell_lines))


median_data = sapply(cell_lines, function(cl){
	#cl = cell_lines[[1]]
	bar$tick()
	
	# load cell_line and limit the number of cells
	tmp = dbReadTable(con,  dbQuoteIdentifier(con,cl)) %>% 
		as_tibble()
	cell_count_table <- tmp %>% select(-cellID) %>%
		group_by(cell_line,treatment, time, fileID) %>% 
		summarise(ncells = n()) 
	
	median_table <- tmp %>% select(-cellID) %>%
		group_by(cell_line,treatment, time, fileID) %>% 
		summarise_all(median,na.rm=TRUE) %>% ungroup()
	
	full_join(median_table,cell_count_table,by = c("cell_line", "treatment", "time", "fileID"))
	
},simplify = F)


median_data = do.call("rbind",median_data) %>% as_tibble()

fileID_tbl <- read_csv("./data/FileID_table.csv") %>% select(cell_line,treatment, time,time_course,fileID) %>%
	mutate(time = ifelse(is.na(time),0,time))

# Testing:
# full join fills missing info with NA. 
test_table <- full_join(fileID_tbl, median_data %>% select(cell_line, treatment,  time,fileID))
# we dont find any NA after merging => it is safe to merge
summarise_all(test_table,~any(is.na(.)))

# merge the median values with the information on the time-course:
median_data <- median_data %>% 
	left_join(fileID_tbl,by = c("cell_line", "treatment", "time", "fileID")) %>%
	select(cell_line,treatment,time,fileID,time_course,ncells,everything())  # just ordering columns


# 

# old version:
# this median not have info on number of cells: 
# write_rds(median_data,"./data/median_data/median_all_reporters.rds")
# new version:
write_rds(median_data,"./data/median_data/median_all_reporters_mine_timecourse_fileID_AB.rds")


median_data <- read_rds("./data/median_data/median_all_reporters_mine_timecourse_fileID_AB.rds")
# There are conditions where there are more than one fileID, but they belong to the same condition, they are technical replicates: 
tmp = median_data %>% group_by(cell_line,treatment, time, time_course) %>% summarise(n = n()) %>% filter(n>1)

# for each condition and time-course we compute the weighted mean of the medians. 
averaged_median_data <- median_data %>% gather(reporter,median_value,-cell_line, -treatment,-time, -fileID, -time_course, -ncells) %>%
	group_by(cell_line,treatment, time, time_course,reporter) %>%
	summarise(weighted_mean_median = weighted.mean(median_value,ncells), 
			  ncells_ave = mean(ncells)) %>% spread(reporter,weighted_mean_median)



write_rds(averaged_median_data,"./data/median_data/median_all_reporters_mine_timecourse_AB.rds")

### Quality check --------------------------------------------------------------
# compare the new and the old median
# there are differences for sure:
# - single cell data was normalised/batch-corrected differently

# I expected differences, because 
# (1) the batch correction was done differently in the 2 cases (on single cell level vs on the median level); 


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
	#summarise(value = mean(value, na.rm = TRUE) ) %>% # mean over time_course
	ungroup() %>%
	rename(old_median=value)

median_data = read_rds("./data/median_data/median_all_reporters_mine_timecourse_AB.rds") %>% as_tibble()

new_median = median_data %>% 
	gather(reporter,value,-cell_line,-treatment,-time,-time_course) %>%
	group_by(cell_line,treatment,time,reporter) %>% 
	#summarise(value = mean(value, na.rm = TRUE) ) %>% # mean over fileID
	ungroup() %>%
	rename(new_median=value)
	
	
	
medians = full_join(old_median,new_median, by = c("cell_line", "treatment", "time_course","time", "reporter"))


library(ggplot2)

medians %>% ggplot(aes(old_median,new_median)) + geom_point() + geom_abline(slope = 1,intercept = 0) + 
	coord_equal() + facet_wrap(~reporter)
# Stat 3
medians %>% filter(reporter=="p.STAT3") %>%  ggplot(aes(old_median,new_median)) + geom_point(aes(col=treatment)) + geom_abline(slope = 1,intercept = 0) + 
	coord_equal() + facet_wrap(~time) + ggtitle("STAT3")

medians %>% filter(reporter=="p.STAT3", time %in% c(0,30,40), treatment=="EGF", abs(old_median-new_median)>1  ) %>%  ggplot(aes(old_median,new_median)) + geom_point(aes(col=cell_line)) + geom_abline(slope = 1,intercept = 0) + 
	coord_equal() + facet_wrap(~time) + ggtitle("STAT3")
# outliers for stat3: HCC2157, MCF10F; treatment EGF, time 30 and 40
medians %>% filter(reporter=="p.STAT3", time %in% c(30,40), treatment=="EGF",cell_line %in% c("HCC2157", "MCF10F")  )

###### MEDIAN INTERPOLATION -----------------------------------------------------
## Median data interpolation
# the code will skip this section and data is read from the disk. 

# Load the median data to find suitable cell-lines

median_data <- read_rds("./data/median_data/median_all_reporters_mine_timecourse_AB.rds") %>% as_tibble()

# preparing time-courses:
# 	we interpolate the data to the same timepoints:


# uncomment if we use Marco's data
# median_data <- median_data %>% mutate(time = as.character(time)) %>%
# 	mutate(time = ifelse(time=="0short","0",time)) %>%
# 	mutate(time = ifelse(is.na(time),"0",time)) %>%
# 	mutate(time = as.numeric(as.character(time))) 

median_data%>%pull(time) %>% table()
base_time_EGF = c(0,5.5,7,9,13,17,23,30,40,60)
base_time_inhib = c(0,7,9,13,17,40,60)

# In case of Marco's median data interpolation in 2 steps:
# 1. interpolate for the basal time per time_course
# - first we filter out treatment full, because that has only time 0 so we dont want to interpolate.
# - for each reporter time_course we interpolate for the base_timepoints, defined above
# 2. interpolate (average) between the time_courses


median_interpolated_data = median_data %>% select(-ncells_ave) %>% filter(treatment != "full") %>%
	gather(reporter,value,-cell_line,-treatment,-time,-time_course) %>%
	filter(!is.na(value)) %>% # some her2 and plcg are NA which ruins the interpolation
	group_by(cell_line,treatment,time_course,reporter) %>% # interpolating per time_course!
	nest() %>%
	mutate(interp_value = map2(data,treatment, function(d,tr){
		
		if(tr=="EGF") base_time = base_time_EGF else base_time = base_time_inhib
		
		out = approx(x = d$time,y = d$value,xout =base_time,method ="linear",rule = 2)
		data.frame(time=out$x,value=out$y)
	}
	))%>% unnest(interp_value) %>%
	# bind the EGF treatment: this has only time 0
	bind_rows( median_data %>% select(-ncells_ave) %>%
			   	filter(treatment == "full") %>%
			   	gather(reporter,value,-cell_line,-treatment,-time,-time_course)  
	) %>%
	
	group_by(cell_line,treatment,time,reporter) %>% # averaging over time_course A/B
	summarise(value = mean(value))
 
if(FALSE) saveRDS(median_interpolated_data,"./data/median_data/interpolated_median_allsamples_correct_times.rds")

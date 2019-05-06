
library(tidyverse)
library(DBI)
library(RSQLite)
library(progress)


con <- dbConnect(RSQLite::SQLite(), "./data/cleaned_single_cell_data/single_cell_dream.sqlite")


# find files with single cell data: 
single_cell_files <- list.files("./data/cleaned_single_cell_data/SingleCells/",full.names = T)


sapply(single_cell_files,function(sc_file){
	
	print("Processing file:")
	print(sc_file)
	# import	
	single_cell_raw <- readRDS(sc_file) %>% as_tibble()
	
	# create condition ID
	single_cell_raw$CID = paste(single_cell_raw$cell_line,single_cell_raw$treatment,single_cell_raw$time,sep = "_")
	
	# create a condition matrix from existing cell_lines, treatment and time
	print(">> condition matrix")
	conds <- unique(single_cell_raw$CID)
	conditions = str_split_fixed(conds,"_",3)
	colnames(conditions) <- c("cell_line","treatment","time")
	conditions <- as_tibble(conditions)
	conditions$time <- as.numeric(conditions$time)
	conditions <- bind_cols(CID=conds,conditions)
	
	# write the conditions table to SQL
	dbWriteTable(con, "conditions", conditions, append=TRUE)
	
	# channels will be converted to columns, take care of the naming
	channels = unique(single_cell_raw$channel)
	names(channels) = make.names(channels)	
	
	bar = progress::progress_bar$new(format = "  Processing [:bar] :percent eta: :eta",
									 total = length(conds))
	
	# write the data table to SQL by each conditions
	sapply(conds,function(cid){
		
		bar$tick()
		single_cell_raw %>% filter(CID == cid) %>%   # cid is a nice small amount of cells
			select(-cell_line,-treatment,-time) %>%  # discard extra columns
			spread(.,channel,value) %>% # wide formar
			rename(!!channels) %>%  # take care of bad naming for colnames
			dbWriteTable(con, "measured_data", value = .,append = TRUE)   # write to DBase
		return(cid)
		
	})
})

dbDisconnect(con)

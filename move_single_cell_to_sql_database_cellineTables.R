# created: 7.May 2019
# last modified: 13.May 2019
# exports the single cell data to SQL database
# - each cell-line goes in a separate table. 

library(tidyverse)
library(DBI)
library(RSQLite)
library(progress)
# parallel version of readRDS:
#  -requires pigz: run in terminal> brew install pigz 
# install.packages("remotes")
# remotes::install_github("traversc/trqwe")
library(trqwe)   

con <- dbConnect(RSQLite::SQLite(), "./data/cleaned_single_cell_data/single_cell_dream_cls.sqlite")


# find files with single cell data: 
single_cell_files <- list.files("./data/cleaned_single_cell_data/SingleCellBatchCorrected",full.names = T)


# sc_file = single_cell_files[[8]]
sapply(single_cell_files,function(sc_file){
	
	print("Processing file:")
	print(sc_file)
	# import	
	single_cell_raw <- mcreadRDS(sc_file,mc.cores = 3) %>% as_tibble()
	#single_cell_raw <- readRDS(sc_file) %>% as_tibble()
	
	cell_lines <- unique(single_cell_raw$cell_line)
	# channels will be converted to columns, take care of the naming
	channels = unique(single_cell_raw$channel)
	names(channels) = make.names(channels)	
	
	bar = progress::progress_bar$new(format = "  Processing [:bar] :percent eta: :eta",
									 total = length(cell_lines))
	
	#cl = cell_lines[[1]]
	sapply(cell_lines,function(cl){
		
		bar$tick()
		tmp = single_cell_raw %>% filter(cell_line == cl) %>%   # find celllines
			spread(.,channel,value) %>% # wide formar
			rename(!!channels) 
		
		fieldTypes = unlist( lapply(tmp,class))
		
		dbWriteTable(con, cl, value = tmp,append = TRUE, field.type = fieldTypes)   # write to DBase
		return(cl)
		
	})
})

dbDisconnect(con)

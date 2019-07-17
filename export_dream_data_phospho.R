# export challenge data: single cell phospho

# export: csv-files per cell-lines
# read in single cell data per cell-lines from SQL, then process and then export to csv file
#
# we create a public_data --> for training
# we create a validation_data --> that we need to score the predictions
#
# we also export csv files, that indicate the values to be predicted by the participants.
# this is handled in export_prediction_conditions.R

library(tidyverse)
library(DBI)
library(RSQLite)
library(progress)

target_folder = "./challenge_data"
dir.create(file.path(target_folder,"single_cell_phospho","complete_cell_lines"))
dir.create(file.path(target_folder,"single_cell_phospho","aim_1_1"))
dir.create(file.path(target_folder,"single_cell_phospho","aim_1_2_1"))
dir.create(file.path(target_folder,"single_cell_phospho","aim_2"))

log_file = file("challenge_data/log.txt",open = "a")

# utility functions ------------------------------------------------------------

# end of utility functions -----------------------------------------------------

### Export single cell phospho data --------------------------------------------
# read in by cell-lines and create training and validation data

con <- dbConnect(RSQLite::SQLite(), "./data/cleaned_single_cell_data/single_cell_dream_cls.sqlite")
cell_lines <-  dbListTables(con)

cell_lines  <- cell_lines[cell_lines != "HCC70_2"]  # remove the duplicated

cell_line_sheet <- readxl::read_excel("./data/cell_line_distribution.xlsx",sheet = 1,range = "A1:I69")

bar = progress::progress_bar$new(format = "  Processing [:bar] :percent eta: :eta",
								 total = length(cell_lines))

current_cell_line = cell_lines[[1]]
current_cell_line = cell_lines[[22]]

for(current_cell_line in cell_lines){	
	bar$tick()
	
	print(paste("reading ",current_cell_line))
	cat(current_cell_line,"---\n",file = log_file)
	
	# load cell_line 
	sc_data = dbReadTable(con,  dbQuoteIdentifier(con,current_cell_line)) %>% 
		as_tibble() 
	
	stopifnot(length(unique(sc_data$cell_line))==1) # make sure there is exactly 1 cell-line there
	
	
	# process according to the purpose of the cell-line. 
	
	purpose = cell_line_sheet %>% filter(cell_line == current_cell_line)
	purpose[1,as.logical(is.na(purpose[1,]))] = ""
	
	public_data = sc_data
	validation_data = sc_data
	
	# 1. remove imTOR condition in the public_data from all cell-lines. 
	public_data <- public_data %>% filter(treatment!="imTOR")
	cat("imTOR condition removed from public data \n",file = log_file)
	
	if(purpose$AIM_1_1 =="test"){
		
		cat("AIM 1.1 test cell line \n",file = log_file)
		cat("removing targeted p-sites from public data \n",file = log_file)
		
		public_data[,c("p.ERK", "p.Akt.Ser473.","p.S6","p.HER2", "p.PLCg2")] = NA_real_
		
		cat("imTOR condition removed from validation data \n",file = log_file)
		validation_data = validation_data %>% filter(treatment!="imTOR") %>%
			select(cell_line, treatment,time,cellID,fileID,c("p.ERK", "p.Akt.Ser473.","p.S6","p.HER2", "p.PLCg2"))
		
		cat("writing public and validation datasets \n",file = log_file)
		
		# write out with 6 digit precision
		write_csv(public_data %>% mutate_if(is.double,format,digits=6), 
				  path = file.path(target_folder,'single_cell_phospho',"aim_1_1",paste0(current_cell_line,".csv")))
		
		# write out with 6 digit precision
		write_csv(validation_data %>% mutate_if(is.double,format,digits=6), 
				  path = file.path(target_folder,'validation_data',paste0("AIM_11_",current_cell_line,".csv")))
		
	}else if(purpose$AIM_1_2_1 == "test"){
		cat("AIM 1.2.1 test cell line \n",file = log_file)
		cat("removing iPKC condition from public data \n",file = log_file)
		
		# remove HER2 and PLCg2 from validation, because in some cell-lines they are misssing anyways
		required_columns <- c('cell_line','treatment', 'time', "cellID","fileID",
							  'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
							  'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
							  'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
							  'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
							  'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
							  'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
							  'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
							  'p.STAT3', 'p.STAT5') 
		
		if(current_cell_line %in% c("MDAMB468","MCF12A","BT483")){
			cat("iEGFR condition selected for validation data \n",file = log_file)
			public_data = public_data %>% filter(treatment!="iEGFR")
			validation_data = validation_data %>% filter(treatment=="iEGFR") %>% select(required_columns)
			
		}else if(current_cell_line %in% c("184B5","ZR751","HCC202")){
			cat("iMEK condition selected for validation data \n",file = log_file)
			public_data = public_data %>% filter(treatment!="iMEK")
			validation_data = validation_data %>% filter(treatment=="iMEK") %>% select(required_columns)
			
		}else if(current_cell_line  %in% c("UACC3199","SKBR3","MDAMB231")){
			cat("iPI3K condition selected for validation data \n",file = log_file)
			public_data = public_data %>% filter(treatment!="iPI3K")
			validation_data = validation_data %>% filter(treatment=="iPI3K") %>% select(required_columns)
			
		}else if(current_cell_line  %in% c("HCC1806","Hs578T","HCC1428")){
			cat("iPKC condition selected for validation data \n",file = log_file)
			public_data = public_data %>% filter(treatment!="iPKC")
			validation_data = validation_data %>% filter(treatment=="iPKC") %>% select(required_columns)
			
		}else stop("AIM_1_2_1 cell-lines changed!!!")
		
		
		cat("writing public and validation datasets \n",file = log_file)
		
		# write out with 6 digit precision
		write_csv(public_data %>% mutate_if(is.double,format,digits=6), 
				  path = file.path(target_folder,'single_cell_phospho',"aim_1_2_1",paste0(current_cell_line,".csv")))
		
		# write out with 6 digit precision
		write_csv(validation_data %>% mutate_if(is.double,format,digits=6), 
				  path = file.path(target_folder,'validation_data',paste0("AIM_121_",current_cell_line,".csv")))
		
	}else if(purpose$AIM_1_2_2 == "test"){
		cat("AIM 1.2.2 test cell line \n",file = log_file)
		cat("no training condition ! \n",file = log_file)
		
		validation_data = validation_data %>% filter(treatment=="imTOR") 
		
		cat("writing validation datasets \n",file = log_file)
		
		# write out with 6 digit precision
		write_csv(validation_data %>% mutate_if(is.double,format,digits=6), 
				  path = file.path(target_folder,'validation_data',paste0("AIM_122_",current_cell_line,".csv")))
		
	}else if(purpose$AIM2 == "test"){
		cat("AIM 2 test cell line \n",file = log_file)
		
		cat("keeping only Full condition from public data \n",file = log_file)
		
		public_data = public_data %>% filter(treatment=="full")
		
		cat("writing public datasets \n",file = log_file)
		
		# write out with 6 digit precision
		write_csv(public_data %>% mutate_if(is.double,format,digits=6), 
				  path = file.path(target_folder,'single_cell_phospho',"aim_2",paste0(current_cell_line,".csv")))
		
		# nothing to do with validation data
		# this challenge is handled in export_dream_median_phospho.R
		
	}else {
		cat("complete training data \n",file = log_file)
		
		# write out with 6 digit precision
		write_csv(public_data %>% mutate_if(is.double,format,digits=6), 
				  path = file.path(target_folder,'single_cell_phospho',"complete_cell_lines",paste0(current_cell_line,".csv")))
		
	}
}


close(log_file) 
dbDisconnect(con)


##### Aggregate files ----------------------------------------------------------
# previously we exported the validation per cell-line, now we import them and aggregate.

# AIM 1.1
prediction_data = list.files(file.path(target_folder,"validation_data"),pattern = "AIM_11_",full.names = T) %>% 
	map(read_csv) %>% bind_rows()

write_csv(prediction_data,"./challenge_data/validation_data/AIM_11_data.csv")


# AIM 1.2.1
prediction_data = list.files(file.path(target_folder,"validation_data"),pattern = "AIM_121_",full.names = T) %>% 
	map(read_csv) %>% bind_rows()

write_csv(prediction_data,"./challenge_data/validation_data/AIM_121_data.csv")


# AIM 1.2.2
prediction_data = list.files(file.path(target_folder,"validation_data"),pattern = "AIM_122_",full.names = T) %>% 
	map(read_csv) %>% bind_rows()

write_csv(prediction_data,"./challenge_data/validation_data/AIM_122_data.csv")



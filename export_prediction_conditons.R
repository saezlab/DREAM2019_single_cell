# export files that shows conditions to be predicted. 

library(tidyverse)
library(DBI)
library(RSQLite)
library(progress)

target_folder = "./challenge_data"


### Export single cell phospho data --------------------------------------------
con <- dbConnect(RSQLite::SQLite(), "./data/cleaned_single_cell_data/single_cell_dream_cls.sqlite")
median_data <- read_csv("./challenge_data/median_phospho/median_phospho_data.csv")	


cell_lines <-  dbListTables(con)

cell_lines = setdiff(cell_lines,"HCC70_2")


cell_line_sheet <- readxl::read_excel("./data/cell_line_distribution.xlsx",sheet = 1,range = "A1:I69")

bar = progress::progress_bar$new(format = "  Processing [:bar] :percent eta: :eta",
								 total = length(cell_lines))

current_cell_line = cell_lines[[1]]
current_cell_line = cell_lines[[3]]
current_cell_line = cell_lines[[2]]
current_cell_line = cell_lines[[4]]
current_cell_line = cell_lines[[8]]

for(current_cell_line in cell_lines){	
	bar$tick()
	
	print(paste("reading ",current_cell_line))
	
	purpose = cell_line_sheet %>% filter(cell_line == current_cell_line)
	purpose[1,as.logical(is.na(purpose[1,]))] = ""  # set NA to empty, to be used in logical evaluation
	
	
	
	# load cell_line 
	sc_data = dbReadTable(con,  dbQuoteIdentifier(con,current_cell_line)) %>% 
		as_tibble() 
	
	stopifnot(length(unique(sc_data$cell_line))==1) # make sure there is exactly 1 cell-line there
	
	
	# process according to the purpose of the cell-line. 
	public_data = sc_data
	
	
	if(purpose$AIM_1_1 =="test"){
		# 1. remove imTOR condition in the public_data from all cell-lines. 
		public_data <- public_data %>% filter(treatment!="imTOR")
		
		prediction_templates = public_data %>% select(cell_line,treatment,time,cellID,fileID,c("p.ERK", "p.Akt.Ser473.","p.S6","p.HER2", "p.PLCg2")) %>%
			mutate_at(c("p.ERK", "p.Akt.Ser473.","p.S6","p.HER2", "p.PLCg2"),~NA_real_)
		
		# write out with 6 digit precision
		write_csv(prediction_templates, 
				  path = file.path(target_folder,'prediction_templates',paste0("AIM_11_",current_cell_line,".csv")))
		
		
	}else if(purpose$AIM_1_2_1 == "test"){
		# Predict cells in iPKC conditions
		
		# 1. remove imTOR condition in the public_data from all cell-lines. 
		public_data <- public_data %>% filter(treatment!="imTOR")
		
		reporters = setdiff(colnames(public_data),c("treatment","cell_line","time","cellID","fileID","p.HER2","p.PLCg2"))
		
		# we find first the unique conditions
		# then we expand each conditions with cellID: 1-10k
		
		if(current_cell_line %in% c("MDAMB468","MCF12A","BT483")){
			public_data = public_data %>% filter(treatment=="iEGFR") 
	
		}else if(current_cell_line %in% c("184B5","ZR751","HCC202")){
			
			public_data = public_data %>% filter(treatment=="iMEK")
			
		}else if(current_cell_line  %in% c("UACC3199","SKBR3","MDAMB231")){
			
			public_data = public_data %>% filter(treatment=="iPI3K")
			
		}else if(current_cell_line  %in% c("HCC1806","Hs578T","HCC1428")){
			
			public_data = public_data %>% filter(treatment=="iPKC")
			
		}
		
		
		prediction_templates = public_data %>% 
			select(cell_line, treatment,time, -cellID,-fileID,reporters) %>%
			group_by(cell_line,treatment,time) %>% mutate_at(reporters,~NA_real_) %>%unique() %>%
			nest() %>% mutate(extended_data = map(data,function(data){
				data.frame(cellID = 1:10000,data)	
			})) %>% unnest(extended_data)
			
		# write out with 6 digit precision
		write_csv(prediction_templates, 
				  path = file.path(target_folder,'prediction_templates',paste0("AIM_121_", current_cell_line,".csv")))
		
	
	}else if(purpose$AIM_1_2_2 == "test"){
		
		reporters = setdiff(colnames(public_data),c("treatment","cell_line","time","cellID","fileID","p.HER2","p.PLCg2"))
		
		# we find first the unique conditions
		# then we expand each conditions with cellID: 1-10k
		prediction_templates = public_data %>% filter(treatment=="imTOR") %>% 
			select(cell_line, treatment,time, -cellID,-fileID,reporters) %>%
			group_by(cell_line,treatment,time) %>% mutate_at(reporters,~NA_real_) %>%unique() %>%
			nest() %>% mutate(extended_data = map(data,function(data){
				data.frame(cellID = 1:10000,data)	
			})) %>% unnest(extended_data)
		
		# write out with 6 digit precision
		write_csv(prediction_templates, 
				  path = file.path(target_folder,'prediction_templates',paste0("AIM_122_", current_cell_line,".csv")))
		
	}else if(purpose$AIM2 == "test"){
		# nothing to do.
		# this challenge is handled in export_dream_median_phospho.R
	}
}

dbDisconnect(con)



##### Aggregate files ----------------------------------------------------------
# previously we exported the validation per cell-line, now we import them and aggregate.

# AIM 1.1
temp_files =  list.files(file.path(target_folder,"prediction_templates"),pattern = "AIM_11_",full.names = T)
template_data = temp_files %>% 
	map(read_csv) %>% bind_rows()
# here we added the glob_cellID
write_csv(template_data,"./challenge_data/prediction_templates/AIM_11_template_data.csv")
file.remove(temp_files )

# AIM 1.2.1
temp_files =   list.files(file.path(target_folder,"prediction_templates"),pattern = "AIM_121_",full.names = T)
template_data = temp_files  %>% 
	map(read_csv) %>% bind_rows()

write_csv(template_data,"./challenge_data/prediction_templates/AIM_121_template_data.csv")
file.remove(temp_files )

# AIM 1.2.2
temp_files =  list.files(file.path(target_folder,"prediction_templates"),pattern = "AIM_122_",full.names = T)
template_data = temp_files %>% 
	map(read_csv) %>% bind_rows()

write_csv(template_data,"./challenge_data/prediction_templates/AIM_122_template_data.csv")

file.remove(temp_files )




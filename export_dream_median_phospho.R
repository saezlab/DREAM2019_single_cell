

# 3 parts:
# 1. export data for public use 
#	- remove the conditions used for the scoring
# 2. export prediction coinditions:
#	- export the conditions and reporters masked with NAs for aim2
# 3. export the validation data. 


median_interpolated_data <- read_rds("./data/median_data/interpolated_median_all_reporters_mine.rds")
median_data <- median_interpolated_data %>% spread(reporter,value)
# cell-line in aims
cell_line_sheet <- readxl::read_excel("./data/cell_line_distribution.xlsx",sheet = 1,range = "A1:I69")


# 1. export data for public use --------------------------------------------------
# export median phosphorylation data after removing the conditions that are to be predicted
# Export data for the public: all data except the test conditions   
# we remove imTOR from all cell-lines. 
# then depending which cell-line is used in which AIM, we remove a subset of the data. 


public_data <- median_data %>% filter(treatment != "imTOR")  %>%# remove imTOR condition from all
	group_by(cell_line) %>% nest(.key = "data") %>% 
	mutate(cleaned_data = map2(data,cell_line,function(data,cell_line){
		
		current_cell_line = cell_line
		
		purpose = cell_line_sheet %>% filter(cell_line == current_cell_line)
		purpose[1,as.logical(is.na(purpose[1,]))] = ""  
		
		
		if(purpose$AIM_1_1 =="test"){
			# remove psites to be predicted
			data[,c("p.ERK", "p.Akt.Ser473.","p.S6","p.HER2", "p.PLCg2")] = NA_real_
			
		}else if(purpose$AIM_1_2_1 == "test"){
			# remove all data in condition to be predicted
			
			data = data %>% filter(treatment!="iPKC")
			
		}else if(purpose$AIM_1_2_2 == "test"){
			
			# nothing to do here, imTOR condition already removed. 
			
		}else if(purpose$AIM2 == "test"){
			# providing only the full condition for the test
			data = data %>% filter(treatment=="full")
			
		}
		
		return(data)
	})) %>% unnest(cleaned_data)


write_csv(public_data,path = "./challenge_data/median_phospho/median_phospho_data.csv")	

### 2. export prediction conditions --------------------------------------------
# export the conditions and use NA fo the values that the participants ahve to predict

prediction_data <- median_data %>% filter(treatment != "imTOR")  %>%# remove imTOR condition from all
	group_by(cell_line) %>% nest(.key = "data") %>% 
	mutate(cleaned_data = map2(data,cell_line,function(data,cell_line){
		
		current_cell_line = cell_line
		
		purpose = cell_line_sheet %>% filter(cell_line == current_cell_line)
		purpose[1,as.logical(is.na(purpose[1,]))] = ""  
		
		reporters = colnames(data)[-1:-2]
		
		 if(purpose$AIM2 == "test"){
			# providing only the full condition for the test
			data = data %>% filter(treatment !="full") %>%
				mutate_at(reporters,~NA_real_)
			
		}else return(tibble())
		
		return(data)
		
	})) %>% unnest(cleaned_data)


write_csv(prediction_data,path = "./challenge_data/predict_conditions/AIM_2_median_data.csv")	


### # 3. export the validation data.  --------------------------------------------
# export the conditions for validaion

validation_data <- median_data %>% filter(treatment != "imTOR")  %>%# remove imTOR condition from all
	group_by(cell_line) %>% nest(.key = "data") %>% 
	mutate(cleaned_data = map2(data,cell_line,function(data,cell_line){
		
		current_cell_line = cell_line
		
		purpose = cell_line_sheet %>% filter(cell_line == current_cell_line)
		purpose[1,as.logical(is.na(purpose[1,]))] = ""  
		
		reporters = colnames(data)[-1:-2]
		
		if(purpose$AIM2 == "test"){
			# providing only the full condition for the test
			data = data %>% filter(treatment !="full") 
			
		}else return(tibble())
		
		return(data)
		
	})) %>% unnest(cleaned_data)


write_csv(validation_data,path = "./challenge_data/validation_data/AIM_2_median_data.csv")	

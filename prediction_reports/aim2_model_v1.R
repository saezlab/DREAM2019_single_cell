# aim2_model_v1: for each reporter we build a model condition-wise, independently
aim2_model_v1 <- function(median_interpolated_data,proteomics,training_cellines){
	library(ranger)
	library(progress)
	
	bar <- progress::progress_bar$new(total = 2257)
	
	model_results <- median_interpolated_data %>% filter(!is.na(value)) %>%
		mutate(data_purpose = ifelse(cell_line %in% training_cellines, "training","test")) %>%
		group_by(treatment,time,reporter) %>% nest() %>% 
		mutate(model = map(data,function(df,prot){
			#browser()
			bar$tick()
			df = df %>% inner_join(prot,by="cell_line")
			
			trainset <- df %>% filter(data_purpose == "training")
			testset <- df %>% filter(data_purpose == "test")
			
			rg.res <- ranger(value ~ ., data = trainset[,c(-1,-3)], mtry = 5, max.depth = 6,verbose = T )
			
			testset <- testset %>% select(cell_line,value,data_purpose) %>% 
				bind_cols(predicted_value = predictions(predict(rg.res, data = testset[,c(-1,-2,-3)])))
			
			trainset <- trainset %>% select(cell_line,value,data_purpose) %>% 
				bind_cols(predicted_value = rg.res$predictions)
			
			
			bind_rows(testset,trainset)
			
			
		},prot=proteomics)) 
	return(model_results)
}
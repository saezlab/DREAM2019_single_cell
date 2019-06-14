# aim2_model_v2: reduce the proteomics only to the measured proteins

aim2_model_v3 <- function(median_interpolated_data,training_cellines,proteomics){
	library(ranger)
	library(progress)
	
	
	
	AB_panel <- readxl::read_excel("./data/antibodyPanel_uniprot.xlsx",sheet = 1,skip = 1)
	nodes <- AB_panel$`UniProt Entry`
	nodes = nodes[!is.na(nodes)]
	
	
	
	
	# proteomics_neighbourhood_indicator: figure out which proteins are in the neighbourhood 
	# of the measured p-Sites. 
	proteomics_indicator = tibble(UID = names(proteomics)) %>% rowwise() %>%
		 mutate(included = any(strsplit(UID,".",fixed = T)[[1]] %in% nodes))
	
	# we keep only a subset of the proteins: 618
	reduced_proteomics = bind_cols(cell_line=proteomics$cell_line,
								   proteomics[,proteomics_indicator$included])
	
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
			
			
		},prot=reduced_proteomics)) 
	return(model_results)
}
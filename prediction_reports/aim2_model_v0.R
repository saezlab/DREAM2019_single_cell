# aim2_model_v0: for each reporter we build a mean prediction from the training

aim2_model_v0 <- function(median_interpolated_data,training_cellines){
	
	
	bar <- progress::progress_bar$new(total = 2257)
	
	model_results <- median_interpolated_data %>% filter(!is.na(value)) %>%
		mutate(data_purpose = ifelse(cell_line %in% training_cellines, "training","test")) %>%
		group_by(treatment,time,reporter) %>% nest() %>% 
		mutate(model = map(data,function(df){
			
			bar$tick()
			
			
			df$predicted_value = df %>% filter(data_purpose == "training") %>% 
				summarise(res = mean(.$value)) %>% pull(res)
			return(df)
			
		})) 
	return(model_results)
}

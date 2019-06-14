# aim2_model_v5: 
# comparing to v4 here we provide the "full" condition for all the cell-lines!
# then as in v4: 
# - compute the average response from the training data and remove it from the test, so 
# the random forest will focus on predicting the difference from the mean. 
# - build model for more than a single condition
# - reduce the proteomics only to the measured proteins

aim2_model_v5 <- function(median_data,training_cellines,proteomics,recompute = FALSE){
	library(ranger)
	library(progress)
	
	
	### !!!! BIG CHANGE:
	# convert the "full" condition to a feature!!!
	
	# convert the "full" condition to a feature!!!
	full_condition = median_data %>% filter(treatment=="full") %>% rename(full_cond = value) %>% select(cell_line,reporter,full_cond)
	# remove the data in full condition from the targets
	median_data <- median_data %>% filter(treatment != "full")
	# join back the feature
	median_data <- median_data %>% left_join(full_condition,by = c("cell_line", "reporter"))

	all_data <- median_data
	
	basic_features <-  median_data %>%  
		mutate(stimulated = time >0) %>% 
		model.matrix( ~0+treatment + time  + stimulated + full_cond ,data = .) %>% as_tibble()
	
	all_features = c(names(basic_features),names(proteomics)[-1])
	
	
	if(recompute){
		# we elimintate the values for the test, data so it is not used in any ways.
		# we keep it together with the training data, because the feature matrix (the contrast matrix)
		
		median_data <- median_data %>% mutate(value = ifelse(data_purpose =="test",NA_real_,value))
		
		
		# part 1.
		# compute the "mean" and "difference from the mean" over training data;
		# this is the first proxy of the test data
		median_data_with_mean <- median_data %>%
			group_by(treatment,time,reporter) %>% 
			mutate(mean_value = mean(value,na.rm = TRUE),
				   diff_value = value - mean(value,na.rm = TRUE)) %>% ungroup()
		
		## check: 
		# for the test data we already got the mean_value, but the diff_from mean is NA (because we dont know the value for the test)
		# median_data_with_mean %>% filter(data_purpose=="test")
		
		# combine training with proteomics
		median_data_with_proteomics = median_data_with_mean %>% left_join(proteomics,by="cell_line")
		
		
		# combine features in one table: conditons and proteomics
		rg_data = bind_cols(median_data_with_proteomics,basic_features) 
		
		
		# formula for the random forest
		f = as.formula( paste( "diff_value ~", paste( all_features, collapse = "+" ) ) )
		
		
		bar <- progress::progress_bar$new(total = length(unique(rg_data$reporter)))
		
		rg <- rg_data %>%
			group_by(reporter) %>% nest() %>% 
			mutate(model = map(data,function(df){
				#browser()
				bar$tick()
				
				
				df.training = df %>% filter(data_purpose == "train")
				# df.test = df %>% filter(data_purpose == "test")
				
				
				rg.res <- ranger(formula = f , data = df.training,
								 mtry = 1500, # make sure to use most of the 
								 max.depth = 6,
								 verbose = T ,
								 always.split.variables=names(basic_features) )
				
				return(rg.res)
				# 
				# df %>% bind_cols(pred_diff = rg.res$predictions) %>% 
				# 	filter(cell_line %in% training_cellines[6:21]) %>%
				# 	ggplot() + geom_line(aes(time,diff_value,col=treatment)) + 
				# 	geom_line(aes(time,pred_diff,col=treatment),linetype=2) + facet_grid(cell_line~treatment)
				# 
				# 
				# testset <- testset %>% select(cell_line,value,data_purpose) %>% 
				# 	bind_cols(predicted_value = predictions(predict(rg.res, data = testset[,c(-1,-2,-3)])))
				# 
				# trainset <- trainset %>% select(cell_line,value,data_purpose) %>% 
				# 	bind_cols(predicted_value = rg.res$predictions)
				# 
				# 
				# bind_rows(testset,trainset)
				
				
			}))
		
		
		
		if(FALSE) write_rds(rg,"./prediction_reports/aim2_model_v5.rds")
	}else{
		rg = read_rds("./prediction_reports/aim2_model_v5.rds")
	}
	
	# 
	# data = rg %>% pluck(2,1)
	# model = rg %>% pluck(3,1)
	rg_predictions <- rg %>%
		mutate(R2_OOB = map(model,function(model){
			model$r.squared
		})) %>%
		mutate(predictions = pmap(., function(data, model, ...) {
			
			pred.test = data %>% filter(data_purpose == "test") %>% 
				mutate(diff_predicted = predictions(predict(model,data = .)),
					   predicted_value = mean_value + diff_predicted) %>%
				select(-all_features,-time1,time)
			
			
			
			pred.training = data %>% filter(data_purpose == "train") %>% 
				mutate(diff_predicted = model$predictions,
					   predicted_value = mean_value + diff_predicted) %>%
				select(-all_features,-time1,time)
			
			
			bind_rows(pred.training,pred.test)
		}))
	
	#rg_predictions %>% unnest(R2_OOB) %>% print(n=37)
	
	
	# we compare the real data and the predictions
	model_v5_predictions <- rg_predictions %>% unnest(predictions) %>% select(-value) %>%
		left_join(all_data,
				  by = c("reporter", "cell_line", "treatment", "data_purpose", "time")) 
	
	
}



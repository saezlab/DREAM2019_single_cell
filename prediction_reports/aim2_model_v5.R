# aim2_model_v5: 
# comparing to v4 here we provide the "full" condition for all the cell-lines!
# then as in v4: 
# - compute the average response from the training data and remove it from the test, so 
# the random forest will focus on predicting the difference from the mean. 
# - build model for more than a single condition
# - reduce the proteomics only to the measured proteins

aim2_model_v5 <- function(median_data,test_data,test_celllines,training_cellines,proteomics,recompute = FALSE){
	library(ranger)
	library(progress)
	
	# we derive the model features from the input data. 
	# (1) take the median values measured at full condition
	# (2) other basic features (time, treatment)
	# this we have do for both the training and the predictive set. We do it at 
	# the same time, to make sure that the model matrix of the two data is compatible. 
	
	# bind the training and the test data. 
	all_data <- bind_rows(median_data,test_data)
	
	# (1) convert the "full" condition to a feature
	
	full_condition = all_data  %>% 
		filter(treatment=="full") %>%
		rename(full_cond = value) %>%
		select(cell_line,reporter,full_cond)
	
	
	all_data <- all_data %>% filter(treatment != "full") %>%
		left_join(full_condition,by = c("cell_line", "reporter"))
	
	# (2) basic features
	basic_features <-  all_data %>%  
		mutate(stimulated = time >0) %>% 
		model.matrix( ~0+treatment + time  + stimulated + full_cond ,data = ., na.action=na.fail()) %>% as_tibble()
	
	all_features = c(names(basic_features),setdiff(names(proteomics),"cell_line"))
	
	
	if(recompute){
		# we elimintate the values for the test, data so it is not used in any ways.
		# we keep it together with the training data, because the feature matrix (the contrast matrix)
		
		#median_data <- median_data %>% mutate(value = ifelse(purpose =="test",NA_real_,value))
		
		
		# part 1.
		# compute the "mean" and "difference from the mean" over training data;
		# this is the first proxy of the test data
		median_data_with_mean <- all_data %>%
			group_by(treatment,time,reporter) %>% 
			mutate(mean_value = mean(value,na.rm = TRUE),
				   diff_value = value - mean(value,na.rm = TRUE)) %>% ungroup()
		
		## check: 
		# for the test data we already got the mean_value, but the diff_from mean is NA (because we dont know the value for the test)
		# median_data_with_mean %>% filter(purpose=="test")
		
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
				
				
				df.training = df %>% filter(purpose == "train")
				# df.test = df %>% filter(data_purpose == "test")
				
				
				rg.res <- ranger(formula = f , data = df.training,
								 mtry = 1500, # make sure to use most of the 
								 max.depth = 6,
								 verbose = T ,
								 always.split.variables=names(basic_features) )
				
				return(rg.res)
				
			}))
		
		
		
		if(FALSE) write_rds(rg,"./prediction_reports/aim2_model_v5.rds")
	}else{
		rg = read_rds("./prediction_reports/aim2_model_v5.rds")
	}
	
	# 
	# data = rg %>% pluck(2,1)
	# model = rg %>% pluck(3,1)
	bar <- progress::progress_bar$new(total = length(rg$reporter))
	rg_predictions <- rg %>%
		mutate(R2_OOB = map(model,function(model){
			model$r.squared
		})) %>%
		mutate(predictions = pmap(., function(data, model, ...) {
			
			bar$tick()
			pred.test = data %>% filter(purpose == "predict") %>% 
				mutate(diff_predicted = predictions(predict(model,data = .)),
					   predicted_value = mean_value + diff_predicted) %>%
				select(-all_features,-time1,time)
			
			
			
			pred.training = data %>% filter(purpose == "train") %>% 
				mutate(diff_predicted = model$predictions,
					   predicted_value = mean_value + diff_predicted) %>%
				select(-all_features,-time1,time)
			
			
			bind_rows(pred.training,pred.test)
		}))
	
	#rg_predictions %>% unnest(R2_OOB) %>% print(n=37)
	
	
	# we compare the real data and the predictions
	model_v5_predictions <- rg_predictions %>% unnest(predictions) %>% select(-value) %>%
		left_join(all_data,
				  by = c("reporter", "cell_line", "treatment", "purpose", "time")) 
	
	
}



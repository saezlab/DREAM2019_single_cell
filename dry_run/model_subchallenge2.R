# Subchallenge II
# we build a random forest model from some cell-lines using a subset of cells from each conditions.
# Then we import the test data, merge with the basic features and model and do the predictions.
# Finally we export the data in the given format. 

# code is written based on ./prediction_reports/condition_prediction_single_cell.Rmd

library(tidyverse)
library(ranger)

## 1.0 model building --------------------------------------------------------------
training_celllines = c("BT474","CAL148","HBL100","MCF7","MDAMB157","T47D","ZR7530")

	
training_data <-  training_celllines %>% paste0("challenge_data/single_cell_phospho/complete_cell_lines/",., ".csv") %>%
	map(read_csv) %>%
 	bind_rows()

reporters = setdiff(names(training_data)[c(-1:-5)],c("value","target"))

# use subset of the data
set.seed(1235)
training_data <- training_data %>% group_by(cell_line, treatment, time) %>% sample_n(min(500,n()))


# data_to_stats: computes the mean and covariance matrix of the single cells, 
# reshapes the matrices for long-tables
# 
data_to_stats <- function(single_cell_data){
    
	# first we compute the mean and the covariance of the reporters
	single_cell_stats <-  single_cell_data %>% select(-fileID,-cellID) %>%
		group_by(cell_line,treatment,time) %>%
		nest() %>%
		mutate(mean_values = map(data,colMeans)) %>%
		mutate(cov_values = map(data,cov))
	
	flattenCovMatrix <- function(covmat) {
		ut <- upper.tri(covmat,diag = TRUE)
		tibble(
			stat_variable = paste0("cov_", rownames(covmat)[row(covmat)[ut]],"_",rownames(covmat)[col(covmat)[ut]]),
			stat_value  =(covmat)[ut]
		)
	}
	
	# we reshape the statistics to a column, these will be estimated. 
	# first the mean, then the cov matrix, finally we bind them 
	single_cell_stats_long <- single_cell_stats %>% 
		mutate(vec_mean = map(mean_values,function(x){
			# reshape the row vector to a column vector
			df = enframe(x,name = "stat_variable", value = "stat_value")
			df %>% mutate(stat_variable=paste0("mean_",stat_variable))
		}
		)) %>%
		mutate(vec_cov = map(cov_values,flattenCovMatrix)) %>% 
		mutate(all_stats = map2(vec_mean,vec_cov,function(mean,cov){
			rbind(mean,cov)
		})) %>% select(treatment, cell_line, time,all_stats) %>%
	    unnest(all_stats)
	

	stats_to_model =  single_cell_stats_long %>% spread(treatment,stat_value)

}
#training_data_prep <- training_data %>% gather(target,value,reporters) %>% spread(treatment,value)

stats_to_model <- data_to_stats(training_data)

stats_to_model <- stats_to_model %>% filter(!is.na(iPKC)) %>% select(-full) %>% mutate(iPI3K = ifelse(is.na(iPI3K),0,iPI3K))


# formula for the random forest
f_pkc = as.formula("iPKC ~ EGF + iEGFR + iMEK + iPI3K")
f_egfr = as.formula("iEGFR ~ EGF + iPKC + iMEK + iPI3K")
f_mek = as.formula("iMEK ~ EGF + iEGFR + iPKC + iPI3K")
f_pi3k = as.formula("iPI3K ~ EGF + iEGFR + iMEK + iPKC")


train_model <- function(df,model_formula){
	
	rg.res <- ranger(formula = model_formula , data = df,
					 mtry = 3, # make sure to use most of the 
					 max.depth = 3,
					 verbose = TRUE )
	return(rg.res)
}



recompute = FALSE
if(recompute){
	# training RF for each predictor
	
	training_res <- stats_to_model %>% group_by(stat_variable) %>% nest() %>%
		mutate(model_pkc  = map(data,train_model,f_pkc),
			   model_egfr = map(data,train_model,f_egfr),
			   model_mek  = map(data,train_model,f_mek),
			   model_pi3k = map(data,train_model,f_pi3k))
	
	# save model 
	if(FALSE) write_rds(training_res,"./dry_run/sc2_models.rds")
}else{
	training_res = read_rds("./dry_run/sc2_models.rds")
}


### predict the test data ------------------------------------------------------


# load the available data for the test cell-lines:
test_data <- list.files("./challenge_data/single_cell_phospho/subchallenge_2/", full.names = TRUE)  %>%
	map(read_csv) %>%
	bind_rows()

# load the conditions for predicting:
template_data <- read_csv("./challenge_data/prediction_templates/subchallenge_2_template_data.csv")  

# combine
prediction_data <- test_data %>% select(cell_line,treatment,time,cellID,-fileID,everything()) %>%
	bind_rows(template_data)

# we prepare the test data: compute the statistics (mean and correlation) and then we use the trained models
# to infer the statistics for the iPKC condition. 
prediction_data_stats <- data_to_stats(prediction_data)

what_to_predict <- template_data %>% select(cell_line,treatment) %>% unique() %>% rename(target = treatment)

prediction_data_stats <- prediction_data_stats %>% left_join(what_to_predict,by="cell_line")
###

# combine the test statistics with the trained models and predict the statistics in
# 
predict_missing_stats <- function(test_data, model_pkc, model_egfr, model_mek, model_pi3k,...){
	df = test_data
	#browser()
	# detect the target column and choose model
	target_col = which.max(colSums(is.na(df)))
	
	selected_model <-	switch (df$target[[1]],
		 iMEK = model_mek,
			iPKC = model_pkc,
			iEGFR = model_egfr,
			iPI3K = model_pi3k
	)	
		
	
	# replace NA values with 0. These line wont be used in the validations
	df <- df %>% mutate_each(function(x)ifelse(is.na(x),0,x))
	
	df[,target_col] <- predictions(predict(selected_model,data = df))
	return(df)
}


predicted_stats <-  prediction_data_stats %>% 
	select(-full) %>%  # full condition is not used in the prediction
	group_by(stat_variable,cell_line) %>% nest(.key = "test_data") %>%  # there is a model for each statistics
	full_join(training_res,by="stat_variable")  %>% # join the trained model to the test data
	# mutate(predicted_data = pmap(.,.f =predict_missing_stats(test_data, model_pkc, model_egfr, model_mek, model_pi3k)))
	mutate(predicted_data = pmap(.,.f =predict_missing_stats))


# reconstruct the mean vectors and the covariance matrix:

predicted_stats_matrix_form <-  predicted_stats %>% unnest(predicted_data) %>% arrange(cell_line,time) %>%
	mutate(is_mean_stat = grepl("mean_",stat_variable)) %>%
	group_by(cell_line,time) %>% nest() %>%
	# return the mean statistics for each variable in the targeted condition
	mutate(mean_vec = map(data,function(d){
		
		target_var = d$target[[1]]
		
		d_new = d %>% filter(is_mean_stat) %>% 
			mutate(variable=gsub("mean_","",stat_variable)) %>% 
			select(variable, target_var, target) %>% arrange(variable)
		names(d_new) = c("variable","prediction","target")
		return(d_new)
	})) %>% 
	# return the Covariance matrix for the iPKC condition
	mutate(cov_mat = map(data,function(d){
		#d = tmp$data[[1]]
		target_var = d$target[[1]]
		
		# find the components of the covariance matrix:
		Mcorr_upper = d %>% filter(!is_mean_stat) %>%
			mutate(variables=gsub("cov_","",stat_variable)) %>%
			separate(variables,into = c("var1","var2"),sep = "_") %>%
			select(var1,var2,target_var) %>% 
			arrange(var1,var2) %>%
			spread(var2,target_var) %>% column_to_rownames("var1") %>% as.matrix()
			
		# convert the upper triangular matrix to a symmtric , full matrix
		Mcorr_upper [lower.tri(Mcorr_upper )]  <- t(Mcorr_upper )[lower.tri(Mcorr_upper )]
		# check : any(Mcorr_upper - t(Mcorr_upper) > 0)
		return(Mcorr_upper)
	}))


predicted_cells <- predicted_stats_matrix_form %>%
	mutate(condid = 1:n()) %>% 
	mutate(single_cell_predictions = pmap(.l = ., .f = function(mean_vec,cov_mat,condid,...){
		print(condid)
		cells = MASS::mvrnorm(n = 10000,Sigma = cov_mat,mu = mean_vec$prediction, tol = 1)
		as_tibble(cells) %>% mutate(cellID=1:10000)
	})) %>% select(cell_line,time,single_cell_predictions) %>% unnest(single_cell_predictions)

predicted_cells %>% left_join(what_to_predict, by = "cell_line") %>% 
	rename(treatment=target) %>% write_csv("./dry_run/sc2_predictions.csv")



source("./scoring_scripts/score_sc2.R")
source("./scoring_scripts/validate_sc2.R")

validate_sc2("./dry_run/sc2_predictions.csv","./challenge_data/validation_data/sc2gold.csv")
score_sc2("./dry_run/sc2_predictions.csv","./challenge_data/validation_data/sc2gold.csv")

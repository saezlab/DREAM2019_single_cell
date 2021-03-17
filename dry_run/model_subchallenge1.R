# Subchallenge I
# we build a random forest model from some cell-lines using a subset of cells from each conditions.
# Then we import the test data, merge with the basic features and model and do the predictions.
# Finally we export the data in the given format. 

library(tidyverse)
library(ranger)

## 1.0 model building --------------------------------------------------------------
training_celllines = c("BT474","CAL148","HBL100","MCF7","MDAMB157","T47D","ZR7530")

prediction_targets = c("p.ERK", "p.Akt.Ser473.","p.S6","p.HER2", "p.PLCg2")


training_data <-  training_celllines %>% paste0("challenge_data/single_cell_phospho/complete_cell_lines/",., ".csv") %>%
	map(read_csv) %>%
 	bind_rows()


# use subset of the data
set.seed(1235)
training_data <- training_data %>% group_by(cell_line, treatment, time) %>% sample_n(min(500,n()))


training_data <- training_data %>% gather(target,value,prediction_targets)

known_reporters = setdiff(names(training_data)[c(-1:-5)],c("value","target"))

basic_features <-  training_data %>%  
	mutate(starved =  treatment!="full",
		   stimulated = treatment!="full" & time >0) %>% 
	model.matrix( ~0+treatment + time + starved + stimulated  ,data = .) %>% as_tibble()


rg_data <- bind_cols(training_data,basic_features) 

all_features = c(names(basic_features),known_reporters)

# formula for the random forest
f = as.formula( paste( "value ~", paste( all_features, collapse = "+" ) ) )

recompute = FALSE
if(recompute){
	# training RF for each predictor
	training_res <- rg_data %>% group_by(target) %>% nest() %>%
		mutate(model = map(data,function(df){
			
			rg.res <- ranger(formula = f , data = df,
							 mtry = 10, # make sure to use most of the 
							 max.depth = 6,
							 verbose = T ,
							 always.split.variables=names(basic_features) )
			
			return(rg.res)
		}))
	
	# save model 
	if(FALSE) write_rds(training_res,"./dry_run/subchallenge1_models.rds")
}else{
	training_res = read_rds("./dry_run/subchallenge1_models.rds")
}

### predict the test data ------------------------------------------------------

# load the data
test_data <- read_csv("./challenge_data/validation_data/sc1gold.csv")

# we have to prepare it with the features as in the training set
test_data <- test_data %>% gather(target,value,prediction_targets)

known_reporters = setdiff(names(test_data)[c(-1:-5)],c("value","target"))

basic_features_test <-  test_data %>%  
	mutate(starved =  treatment!="full",
		   stimulated = treatment!="full" & time >0) %>% 
	model.matrix( ~0+treatment + time + starved + stimulated, data = . , ) %>% as_tibble()

stopifnot(all(names(basic_features_test) == names(basic_features)))

rg_testdata <- bind_cols(test_data,basic_features_test) 


# training RF for each predictor
test_res <- rg_testdata %>% 
	group_by(target) %>%
	nest(.key = "test_data") %>%
	inner_join(training_res, by = "target")


test_res <- test_res %>% 
	mutate(test_data = pmap(.,function(test_data,model,...){
		#model = test_res$model[[1]]
		#test_data = test_res$test_data[[1]]
		model_preds = predict(model,data = test_data, verbose = TRUE)
		test_data$value = predictions(model_preds)
		
		return(test_data)
	}))

test_to_export <- test_res %>% unnest(test_data) 
test_to_export <- test_to_export %>% select(cell_line,treatment,time,cellID,fileID,target,value) %>% spread(target,value)

write_csv(test_to_export,"./dry_run/subchallenge1_predictions.csv")

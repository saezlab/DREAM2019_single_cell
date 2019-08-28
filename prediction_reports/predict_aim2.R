# we play around the AIM 2 idea. 


library(ranger)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)


### Data import

median_interpolated_data <- read_csv("./challenge_data/median_phospho/median_phospho_data.csv")
proteomics_raw <- read_csv("./challenge_data/proteomics/proteomics.csv") 
cell_lines_usage <- read_csv("./challenge_data/cell_line_usage.csv")
#rnaseq_raw <- readRDS("./data/genomics/dat_RNAseq_Marcotte.rds") %>% as_tibble()


# preprocessing data sets -----------------------------------------------------
proteomics <- proteomics_raw %>% 
	mutate(cell_line = gsub("normal_vs_","",Label)) %>%
	select(Protein,cell_line,log2FC) %>%
	mutate(log2FC_mod = ifelse(is.infinite(log2FC),NA,log2FC)) %>%
	select(-log2FC) %>%
	spread(cell_line,log2FC_mod) %>%
	 filter(complete.cases(.)) %>% 
	gather(cell_line,value,-Protein) %>% 
	spread(Protein,value)


# setup the model inputs ----------

# no need to build model for Her2 and plcgamma
median_interpolated_data$p.HER2 = NULL
median_interpolated_data$p.PLCg2 = NULL


# celllines with proteomics data:
cell_lines =  unique(median_interpolated_data$cell_line) 
cell_lines <- cell_lines[cell_lines %in% proteomics$cell_line]

colnames(proteomics) <- make.names(colnames(proteomics))

# Finding cells for training 

training_cellines = cell_lines_usage %>% filter(complete==TRUE) %>% pull(cell_line)

# test: 
test_data <- read_csv("challenge_data/prediction_templates/subchallenge_4_template_data.csv") %>% 
	gather(reporter,value,-cell_line, -treatment,  -time) %>% mutate(purpose = "predict")
test_celllines <- test_data %>% pull(cell_line) %>% unique()

median_data = median_interpolated_data %>% 
	filter(cell_line %in% c(training_cellines,test_celllines)) %>%
	mutate(purpose =ifelse(cell_line %in% training_cellines, "train","predict")) %>%
	gather(reporter,value,-cell_line, -treatment,  -time, -purpose) 


# model_v1: for each reporter we build a model condition-wise, independently
source("prediction_reports/aim2_model_v5.R")
model_v5_results = aim2_model_v5(median_data = median_data,
								 test_celllines = test_celllines,
								 test_data = test_data,
								 proteomics = proteomics,
								 training_cellines = training_cellines,
								 recompute = TRUE)


model_v5_predictions %>% write_csv("./dry_run/aim2_model_predictions_train_test.csv")
model_v5_predictions %>% filter(purpose == "predict") %>%
	select(cell_line,treatment,time,reporter,predicted_value) %>%
	spread(reporter,predicted_value) %>% write_csv("./dry_run/aim2_model_predictions_test.csv")



##### plot prediction results for model_v1 -------------------------------------
pdf("./prediction_reports/aim2_predictions.pdf")
sapply(cell_lines,function(cl){
	#cl = "BT20"
	data = 	model_results %>% unnest(model) %>% 
		filter(cell_line == cl)
	gg <- data %>% 
		ggplot() + 
		geom_line(aes(time,value,col=treatment)) +
		geom_line(aes(time,predicted_value,col=treatment),linetype=2) +
		facet_wrap(~reporter) + 
		theme_bw() + 
		ggtitle(paste(cl,",",data$data_purpose[[1]]))
	print(gg)
	return(NULL)
})
dev.off()

pdf("./prediction_reports/aim2_predictions_reporterwise.pdf")
sapply(unique(model_results$reporter),function(ireporter){
	#cl = "p.S6"
	data = 	model_results %>% unnest(model) %>% 
		filter(reporter == ireporter)
	gg <- data %>% 
		ggplot() + 
		geom_line(aes(time,value,col=treatment)) +
		geom_line(aes(time,predicted_value,col=treatment),linetype=2) +
		facet_wrap(~cell_line) + 
		theme_bw() + 
		ggtitle(paste(ireporter,",",data$data_purpose[[1]]))
	print(gg)
	return(NULL)
})
dev.off()


pdf("./prediction_reports/aim2_predictions_corrplot.pdf")
sapply(cell_lines,function(cl){
	#cl = "BT20"
	data = 	model_results %>% unnest(model) %>% 
		filter(cell_line == cl)
	gg <- data %>% 
		ggplot() + 
		geom_point(aes(predicted_value,value,col=treatment)) +
		facet_wrap(~reporter) + 
		theme_bw() + 
		ggtitle(paste(cl,",",data$data_purpose[[1]]))
	print(gg)
	return(NULL)
})
dev.off()

#### statistics of prediction results for model_v1 -----------------------------
model_condition_stats <- model_results %>% unnest(model) %>% 
	group_by(cell_line,treatment,reporter) %>%
	summarise(
		RMSE = sqrt(sum((value-predicted_value)^2)/n()),
		R2 = 1- sum((value-predicted_value)^2) / sum((value-mean(value))^2),
		data_sd = sd(value)
		) 
model_global_stats <- model_results %>% unnest(model) %>% ungroup() %>% 
	summarise(
		RMSE = sqrt(sum((value-predicted_value)^2)/n()),
		R2 = 1- sum((value-predicted_value)^2) / sum((value-mean(value))^2),
		data_sd = sd(value)
	) 

model_cell_line_stats <- model_results %>% unnest(model) %>% ungroup() %>% 
	group_by(cell_line) %>%
	summarise(
		RMSE = sqrt(sum((value-predicted_value)^2)/n()),
		R2 = 1- sum((value-predicted_value)^2) / sum((value-mean(value))^2),
		data_sd = sd(value)
	) 

model_reporter_stats <- model_results %>% unnest(model) %>% ungroup() %>% group_by(reporter) %>%
	summarise(
		RMSE = sqrt(sum((value-predicted_value)^2)/n()),
		R2 = 1- sum((value-predicted_value)^2) / sum((value-mean(value))^2),
		data_sd = sd(value)
	) 



# reduce the dimensionality of the proteomics by selecting the proteins close to the reporters. 

AB_panel <- readxl::read_excel("./data/antibodyPanel_uniprot.xlsx",sheet = 1,skip = 1)
nodes <- AB_panel$`UniProt Entry`

library(OmnipathR)

OP_interactions <- import_Omnipath_Interactions()
# find the nodes in the neighbourhood network:
neigbourhood <- OP_interactions %>% as_tibble() %>% 
	filter(source %in% nodes | target %in% nodes) %>% select(source,target) %>% as.list() %>% unlist() %>% unique()




















# we play around the AIM 2 idea. 


library(ranger)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)


### Data import

median_interpolated_data <- read_csv("./challenge_data/median_phospho/median_phospho_data.csv")
proteomics_raw <- read_csv("./challenge_data/proteomics/Proteomics_log2FC.csv") 
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


# run model
source("prediction_reports/aim2_model_v5.R")
model_v5_results = aim2_model_v5(median_data = median_data,
								 test_celllines = test_celllines,
								 test_data = test_data,
								 proteomics = proteomics,
								 training_cellines = training_cellines,
								 recompute = FALSE)


if(FALSE){ model_v5_results %>% write_csv("./dry_run/aim2_model_predictions_train_test.csv")
    
    model_v5_results %>% filter(purpose == "predict") %>%
        select(cell_line,treatment,time,reporter,predicted_value) %>%
        spread(reporter,predicted_value) %>% write_csv("./dry_run/subchallenge4_predictions.csv")
}



### Score predictions



source("./scoring_scripts/score_sc4.R")
source("./scoring_scripts/validate_sc4.R")

validate_sc4("./dry_run/subchallenge4_predictions.csv","./challenge_data/validation_data/sc4gold.csv")
score <- score_sc4("./dry_run/subchallenge4_predictions.csv","./challenge_data/validation_data/sc4gold.csv")

# score:  0.3408274















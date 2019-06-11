# we play around the AIM 2 idea. 



library(tidyverse)
library(pheatmap)
library(RColorBrewer)


### Data import

median_interpolated_data <- read_rds("./data/median_data/interpolated_median_all_reporters_mine.rds")
proteomics_raw <- readRDS("./data/proteomics/MSstat_groupComparison_selceted.rds") %>% as_tibble()
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


# cell_lines <- rnaseq_raw$Cellline
# rnaseq_raw$Cellline <- NULL
# rnaseq_raw <- type_convert(rnaseq_raw, col_types = cols(.default = col_double()))
# rnaseq_raw <- bind_cols(data.frame(cell_lines=cell_lines), rnaseq_raw)
# rnaseq <- rnaseq_raw %>% gather("gene","value",-cell_lines) %>% as_tibble()



### Proteomics only ------------------------------------------------------------



#### predict each phospho in each condition ----------------------------------
library(ranger)
# celllines with proteomics data:
cell_lines =  unique(median_interpolated_data$cell_line) 
cell_lines <- cell_lines[cell_lines %in% proteomics$cell_line]

colnames(proteomics) <- make.names(colnames(proteomics))

# Spliting training and testing dataset
training_cellines = sample(cell_lines , length(cell_lines) * 0.8, replace = FALSE ) 
# test: 
# test_cellines = c("BT20", "BT474", "CAL148", "EFM19", "EFM192A", "HCC1937", "HCC1954",
#  "HCC2157", "MDAkb2", "MDAMB361", "MDAMB436", "MDAMB453", "OCUBM")   

# genes = colnames(proteomics)
# f = as.formula( paste( "value ~", paste( n[!n %in% c("cell_line")], collapse = "+" ) ) )


# model_v1: for each reporter we build a model condition-wise, independently
source("prediction_reports/aim2_model_v1.R")
model_v1_results = aim2_model_v1(median_interpolated_data,proteomics,training_cellines)

if(FALSE) saveRDS(model_v1_results,"./prediction_reports/aim2_model_predictions.RDS")

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





















# reduce the dimensionality of the proteomics via PCA
#############
pca_trainset = trainset %>% select( -value,-cell_line )
pca_testset = testset
pca = prcomp( pca_trainset )

# variance
pr_var = ( pca$sdev )^2 

# % of variance
prop_varex = pr_var / sum( pr_var )

# Plot
plot( prop_varex, xlab = "Principal Component", 
	  ylab = "Proportion of Variance Explained", type = "b" )

summary(pca)
#plot(cumsum(proteomics_pca$sdev^2/sum(proteomics_pca$sdev^2)))

# Creating a new dataset
train = data.frame( value = trainset$value, pca$x )
t = as.data.frame( predict( pca, newdata = pca_testset ) )

new_trainset = train[, 1:36]
new_testset =  t[, 1:35]

n = names( new_trainset )
f = as.formula( paste( "value ~", paste( n[!n %in% c("value")], collapse = "+" ) ) )

library(ranger)

rg.res1 <- ranger(value ~ ., data = new_trainset, )
pred.training <- predict(rg.res1, data = new_trainset)
pred.test <- predict(rg.res1, data = new_testset)

plot(pred.value$predictions, trainset$value)
points(pred.test$predictions, test$value, col="red")

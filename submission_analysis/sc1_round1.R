# show some results on Subchallenge 1. 
# downloaded the data from the submissinos
# compared it with the golden standard
# - correlation an variance explained. 
# - also tried showing the cells on UMAP, no success; we only get a cloud. 

library(tidyverse)
library(ggplot2)
# best model from SC 1, round 1
# 9693326.csv -- round 1, sc 1 best
prediction <- read_csv("./submission_data/round1/SC1/9693326.csv")
goldStandard <- read_csv("./challenge_data/validation_data/sc1gold.csv")


data = full_join(
	prediction %>% gather(reporter, prediction,7:11),
	goldStandard %>% gather(reporter, measured,7:11),
	by = c("glob_cellID", "cell_line", "treatment", "time", "cellID", "fileID", "reporter"))


data %>% sample_n(1e5) %>% ggplot(aes(measured,prediction,col=time)) + 
	geom_point() + 
	facet_grid(treatment~reporter) + 
	geom_smooth(method="lm",col="red") + 
	theme_bw()


data %>% filter(treatment=="EGF") %>% sample_n(1e5) %>% ggplot(aes(measured,prediction)) + 
	geom_point() + 
	facet_grid(time~reporter) + 
	geom_smooth(method="lm",col="red") + 
	theme_bw() + coord_fixed(ratio = 1) 


# correlation:
corr_sum <- data %>% group_by(cell_line, treatment,reporter) %>% summarise(r = cor(prediction,measured))

hist(corr_sum$r)
mean(corr_sum$r)
sd(corr_sum$r)
RColorBrewer::display.brewer.all()
gg <- corr_sum %>% mutate(condition = paste(cell_line,treatment,sep = "_"))%>%
	ungroup()%>% select(condition,reporter,r) %>% spread(reporter,r) %>% column_to_rownames("condition") %>%
	pheatmap::pheatmap(cluster_cols = F, 
					   cluster_rows = F,
					   breaks = seq(0,1,length.out = 10),
					   colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(11))

ggsave(plot = gg$gtable,filename = "./submission_analysis/figures/round1_sc1_best_correlation_heatmap.pdf", width = 3,height = 6)
# variance explained
R2 <- data  %>% summarise(R2 = 1- sum((prediction-measured)^2)/sum((measured-mean(measured))^2))


data %>% filter(treatment=="EGF", cell_line =="LY2", reporter=="p.ERK" ) %>% 
	gather(source,value,prediction, measured) %>%
	ggplot() + 
	geom_violin(aes(as.factor(time),value,group=paste(time,source),fill=source),
				position = "identity", scale = "width", alpha=0.5) + 
	xlab("time [min]") + ylab("reporter value")  + 
	theme_bw()

data %>% filter(treatment=="iEGFR",cell_line =="LY2", reporter=="p.Akt.Ser473.")  %>% 
	ggplot(aes(measured,prediction)) + 
	geom_point() + 
	theme_bw() + coord_fixed(ratio = 1) 


#install.packages("uwot")
library(uwot)
scdata = bind_rows(
	prediction %>% mutate(source = "prediction"),
	goldStandard %>% mutate(source = "measured"))



sc_sample <- scdata %>% group_by(cell_line,treatment, time,source) %>% sample_frac(0.1) %>% ungroup() 
umap_eu_all <- sc_sample %>% select( 7:11) %>%
	umap(metric="euclidean", init="PCA", n_components = 2,
					fast_sgd=TRUE, min_dist=0.1, verbose = T)


sc_sample %>% bind_cols(umap_eu_all %>% as_tibble()) %>% ggplot() + geom_point(aes(V1,V2,col=source))
sc_sample %>% bind_cols(umap_eu_all %>% as_tibble()) %>% ggplot() + geom_point(aes(V1,V2,col=source)) + facet_grid(cell_line ~ treatment)
sc_sample %>% bind_cols(umap_eu_all %>% as_tibble()) %>% ggplot() + geom_point(aes(V1,V2,col=time)) 

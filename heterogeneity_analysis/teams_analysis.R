# Analysing heterogeneity in submissions in SC2
# This analysis is done by selecting the CC_markers, reducing the expressions values of these 7 markers
# to two dimensions using UMAP and then perform hierarchical clustering on these two dimensions
# The clusters are then visualised by plotting the first two UMAP dimensions and colouring points by the found clusters
# These steps are performed per model, and then per condition or all conditions together
# Also, sumsampling is performed per model and condition to keep the analysis fast.
# Created by Alice Driessen
# Created on: 04-02-20202

library(tidyverse)
library(uwot)
library(fastcluster)
library(ComplexHeatmap)

setwd("~/Desktop/BQ internship/DREAM2019_single_cell/heterogeneity_analysis/SC2")

markers <- c('b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU', 'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 
             'p.AKT.Thr308.', 'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b', 'p.H3', 'p.JNK', 
             'p.MAP2K3', 'p.MAPKAPK2', 'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38', 'p.p53', 
             'p.p90RSK', 'p.PDPK1', 'p.RB', 'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1', 'p.STAT3', 
             'p.STAT5')

nested_predictions <- readRDS("~/Desktop/BQ internship/DREAM2019_single_cell/submission_data/intermediate_data/sc2_all_predictions_nested.rds")

CC_markers <- c("Ki.67", "CyclinB", "GAPDH", "IdU", "p.RB", "cleavedCas", "p.H3")

pred_sub <- nested_predictions %>%
  group_by(team, cell_line, treatment, time) %>%
  mutate(sample = map(data, ~sample_frac(., 0.1))) %>%
  ungroup() %>%
  select(-data) %>%
  unnest(sample)


conditions <- nested_predictions %>%
  select(cell_line, treatment) %>%
  distinct()
for (j in 1:length(unique(nested_predictions$team))) {
  print( unique(nested_predictions$team)[j])
  team_sub <- pred_sub %>%
    filter(team ==  unique(nested_predictions$team)[j])
  
  
  # Clustering and analysis per cell line-treatment condition for selected team
  # For this step team_sub is sampled 0.1 per team and condition
  for (i in 1:dim(conditions)[1]){
    CL <- conditions$cell_line[i]
    TR <- conditions$treatment[i]
    print(paste0("iteration: ", i, ". Cell line: ", CL, " and treatment: ", TR))
    
    # To do UMAP dimensionality reduction and clustering and save result
    if (FALSE) {
     dir.create(paste0(TR, "_", CL))
      
      cond_data <- team_sub  %>%
        filter(cell_line == CL & treatment == TR)
      
      umap_eu_lin <- umap(as.matrix(select(cond_data, CC_markers)), metric="euclidean", init="PCA", 
                          n_components = 2, fast_sgd=TRUE, min_dist=0.1)
      
      hc_umap <- umap_eu_lin %>%
        dist(method = "euclidean") %>%
        fastcluster::hclust(method = "single") %>%
        cutree(5)
      
      cond_data <- cond_data %>% 
        add_column(V1 = umap_eu_lin[,1], 
                   V2 = umap_eu_lin[,2],
                   umap_hc_clust = as.factor(hc_umap))
      
      saveRDS(cond_data, paste0(TR, "_", CL, "/", unique(nested_predictions$team)[j], "_sub_clustered.rds"))}
    
    cond_data <- readRDS(paste0(TR, "_", CL, "/", unique(nested_predictions$team)[j], "_sub_clustered.rds"))
    
    pdf(paste0(TR, "_", CL, "/", unique(nested_predictions$team)[j], "_UMAP_HC_clust.pdf"), height = 5, width = 5)
    print(cond_data %>%
            ggplot(aes(x=V1, y=V2, colour=umap_hc_clust)) +
            geom_point(size=0.1) +
            labs(x = "UMAP V1", y = "UMAP V2", colour = "cluster",
                 title= paste0("UMAP colored by HC on first 2 UMAP comps ", CL, ", ", TR)) +
            guides(colour = guide_legend(override.aes = list(size=4.5))) +
            theme_bw()
    )
    dev.off()
    
    hm <- cond_data %>% 
      gather(marker, marker_value, markers) %>%
      group_by(umap_hc_clust, marker) %>%
      summarise(median_value = median(marker_value)) %>%
      spread(umap_hc_clust, median_value)
    
    pdf(paste0(TR, "_", CL, "/", unique(nested_predictions$team)[j], "_UMAP_HC_clust_HM.pdf"), height = 8, width = 4)
    print(Heatmap(as.matrix(select(hm, -marker)), right_annotation = rowAnnotation(marker = anno_text(hm$marker)), 
                  col = RColorBrewer::brewer.pal(9, "YlOrRd"), column_title = paste0("HC-UMAP clustering, ", CL, ", ", TR), 
                  bottom_annotation = HeatmapAnnotation(cluster_size = anno_barplot(count(cond_data, umap_hc_clust)$n)))
    )
    dev.off()
  }
}

# Analyse all confition together
setwd("~/Desktop/BQ internship/DREAM2019_single_cell/heterogeneity_analysis/SC2")
dir.create("all_conditions")
setwd("./all_conditions")

pred_sub <- nested_predictions %>%
  group_by(team, cell_line, treatment, time) %>%
  mutate(sample = map(data, ~sample_frac(., 0.01))) %>%
  ungroup() %>%
  select(-data) %>%
  unnest(sample)

# Do UMAP dimensionslity reduction and apply hierarchical clustering per team
# Clustering over all cell lines, treatments and timepoints
# For this step team_sub is samples with 0.01 per team and condition
if (TRUE) {
  for (i in 1:length(unique(nested_predictions$team))) {
    print(unique(nested_predictions$team)[i])
    team_sub <- pred_sub %>%
      filter(team == unique(nested_predictions$team)[i])
    
    # To do UMAP dimensionality reduction and clustering and save result
    if (FALSE) {
      umap_eu_lin <- umap(as.matrix(select(team_sub, CC_markers)), metric="euclidean", init="PCA",
                          n_components = 2, fast_sgd=TRUE, min_dist=0.1)
      
      hc_umap <- umap_eu_lin %>%
        dist(method = "euclidean") %>%
        fastcluster::hclust(method = "single") %>%
        cutree(5)
      
      team_clusters <- team_sub %>%
        add_column(V1 = umap_eu_lin[,1],
                   V2 = umap_eu_lin[,2],
                   umap_hc_clust = as.factor(hc_umap))
      saveRDS(team_clusters, paste0(unique(nested_predictions$team)[i], "_sub_clustered.rds"))}
    
    team_clusters <- readRDS( paste0(unique(nested_predictions$team)[i], "_sub_clustered.rds"))
    
    pdf(paste0(unique(nested_predictions$team)[i], "_UMAP_HC_clust.pdf"), height = 5, width = 5)
    print(team_clusters %>%
            ggplot(aes(x=V1, y=V2, colour=umap_hc_clust)) +
            geom_point(size=0.1) +
            labs(x = "UMAP V1", y = "UMAP V2", colour = "cluster",
                 title= "UMAP colored by HC on first 2 UMAP comps") +
            guides(colour = guide_legend(override.aes = list(size=4.5))) +
            theme_bw()
    )
    dev.off()
    
    hm <- team_clusters %>% 
      gather(marker, marker_value, markers) %>%
      group_by(umap_hc_clust, marker) %>%
      summarise(median_value = median(marker_value)) %>%
      spread(umap_hc_clust, median_value)
    
    pdf(paste0(unique(nested_predictions$team)[i], "_UMAP_HC_clust_HM.pdf"), height = 8, width = 4)
    print(Heatmap(as.matrix(select(hm, -marker)), right_annotation = rowAnnotation(marker = anno_text(hm$marker)), 
                  col = RColorBrewer::brewer.pal(9, "YlOrRd"), column_title = "HC-UMAP clustering", 
                  bottom_annotation = HeatmapAnnotation(cluster_size = anno_barplot(count(team_clusters, umap_hc_clust)$n)))
    )
    dev.off()
  }
}



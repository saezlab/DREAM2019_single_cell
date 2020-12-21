

library(tidyverse)

library(uwot)


# Participating team names
submissions <-  readRDS("./submission_analysis/intermediate_data/sc1_ranked_teams.rds") %>% as.character()
sub_data_all <- readRDS("./submission_analysis/intermediate_data/sc1_all_NP_predictions.rds")  



GS_data <- sub_data_all %>% select(1:8) %>%  spread(marker,standard)

team1 <- sub_data_all %>% select(1:7,9) %>%  spread(marker,icx_bxai)

set.seed(124869)
GS_subset_data = GS_data %>% group_by(cell_line,treatment,time) %>% sample_n(min(200,n()))

um_res_data <- uwot::umap(GS_subset_data[,-1:-6], 
                          metric="euclidean", 
                          init="pca", 
                          scale = "none",
                          n_components = 2,
                          fast_sgd=TRUE,
                          min_dist=0.01,
                          n_neighbors = 15,
                          verbose = T)

GS_subset_data %>% bind_cols(tibble(umap.x = um_res_data[,1],
                                   umap.y = um_res_data[,2])) %>% 
    #sample_frac(0.1) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=`p.ERK`),size = 0.2) +
    # coord_fixed() +
    theme_bw() + 
    ggtitle("data") +
    labs(x="",y="") +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) + facet_wrap(~cell_line)




set.seed(124869)
team1_subset_data = team1 %>% group_by(cell_line,treatment,time) %>% sample_n(min(500,n()))

um_res_team1_data <- uwot::umap(team1_subset_data[,-1:-6], 
                          metric="euclidean", 
                          init="pca", 
                          scale = "none",
                          n_components = 2,
                          fast_sgd=TRUE,
                          min_dist=0.01,
                          n_neighbors = 15,
                          verbose = T)



team1_subset_data %>% bind_cols(tibble(umap.x = um_res_team1_data[,1],
                     umap.y = um_res_team1_data[,2])) %>% 
    #sample_frac(0.1) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line),size = 0.2) +
    # coord_fixed() +
    theme_bw() + 
    ggtitle("data") +
    labs(x="",y="") +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())

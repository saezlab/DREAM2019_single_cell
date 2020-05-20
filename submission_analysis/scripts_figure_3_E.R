# SC2 UMAP of participants and real cells. 
# 


library(tidyverse)
library(uwot)
library(ggplot2)
# read golden standard
gs <- read_csv("./challenge_data/validation_data/sc2gold.csv") 

ranked_teams <- read_rds("./submission_analysis/intermediate_data/sc2_ranked_teams.rds")

hide_names <- tibble(submitterId = c("data",as.character(ranked_teams)),
                     alt_name = c("data", as.character(1:length(ranked_teams))))

hide_names$alt_name = factor(hide_names$alt_name,levels =hide_names$alt_name)

submission_folder = "./submission_data/final_round/SC2/"

# read leaderboard
SC_leaderboard = read_csv(file.path(submission_folder,"leaderboard_final_sc2.csv")) %>%
    select(-writeUp, -createdOn) %>% 
    mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
    mutate(submitterId = make.names(submitterId))

# read team's predictions from csv files and compute the stats
required_columns <- c('cell_line','treatment', 'time',
                      'b.CATENIN', 'cleavedCas', 'CyclinB', 'GAPDH', 'IdU',
                      'Ki.67', 'p.4EBP1', 'p.Akt.Ser473.', 'p.AKT.Thr308.',
                      'p.AMPK', 'p.BTK', 'p.CREB', 'p.ERK', 'p.FAK', 'p.GSK3b',
                      'p.H3', 'p.JNK', 'p.MAP2K3', 'p.MAPKAPK2',
                      'p.MEK', 'p.MKK3.MKK6', 'p.MKK4', 'p.NFkB', 'p.p38',
                      'p.p53', 'p.p90RSK', 'p.PDPK1', 'p.RB', 
                      'p.S6', 'p.S6K', 'p.SMAD23', 'p.SRC', 'p.STAT1',
                      'p.STAT3', 'p.STAT5') 

read_sample <- function(file_name){
    read_csv(file_name) %>% select(required_columns) %>% 
        group_by(cell_line,treatment,time) %>% 
        sample_n(.,1000)
}

prediction_data <- SC_leaderboard  %>%
    top_n(-4,score) %>% # select top 4 teams
    mutate(predictions = map(file.path(submission_folder,submissions),read_sample))






#### Data alone ---------------------------------------------------------------
set.seed(12354)
measured_data <- gs %>% add_column(submitterId="data",.before = 1) %>% 
    group_by(cell_line,treatment,time) %>% 
    sample_n(.,min(1000,n())) %>%
    select(submitterId,required_columns)

um_res_data <- uwot::umap(measured_data[,-1:-4], 
                          metric="euclidean", 
                          init="pca", 
                          scale = "none",
                          n_components = 2,
                          fast_sgd=TRUE,
                          min_dist=0.1,
                          n_neighbors = 15,
                          verbose = T)


measured_data %>% bind_cols(tibble(umap.x = um_res_data[,1],
                                   umap.y = um_res_data[,2])) %>% 
    sample_frac(0.1) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line),size = 0.2) +
    # coord_fixed() +
    theme_bw() + 
    ggtitle("data") +
    labs(x="",y="") +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())

ggsave("./publication/figures/figure3/sc2_UMAP_independent_data_legend.pdf", height = 4, width = 5)
ggsave("./publication/figures/figure3/sc2_UMAP_independent_data_legend.jpeg", height = 4, width = 5)


measured_data %>% bind_cols(tibble(umap.x = um_res_data[,1],
                                   umap.y = um_res_data[,2])) %>% 
    sample_frac(0.1) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line),size = 0.2) +
    # coord_fixed() +
    theme_bw() + 
    ggtitle("data") +
    guides(colour = "none") +
    labs(x="",y="") +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())


ggsave("./publication/figures/figure3/sc2_UMAP_independent_data.pdf", height = 4, width = 4)
ggsave("./publication/figures/figure3/sc2_UMAP_independent_data.jpeg", height = 4, width = 4)



####  Prediction alone ---------------------------------------------------------


set.seed(1493)
predictions_umap <- prediction_data %>% 
    mutate(umap_teams = map(predictions, ~ uwot::umap(.x[,-1:-3], 
                                                      metric="euclidean", 
                                                      init="pca", 
                                                      scale = "none",
                                                      n_components = 2,
                                                      fast_sgd=TRUE,
                                                      min_dist=0.1,
                                                      n_neighbors = 15,
                                                      verbose = T)))

# bind the umap coordinates to single cell data: 
predictions_umap_coords <- predictions_umap %>% 
    mutate(tmp = map2(.x = predictions,.y = umap_teams,
                      ~bind_cols(.x,tibble(umap.x = .y[,1],umap.y = .y[,2])) ))


ggplots <- predictions_umap_coords %>%
    mutate(ggs = map(tmp, ~ ggplot(data= sample_frac(.x,0.1)) +
                         geom_point(aes(umap.x,umap.y,col=cell_line),size = 0.1) +
                         theme_bw() + 
                         guides(colour = "none") +
                         labs(x="",y="") +
                         theme(panel.grid = element_blank(),
                               axis.text = element_blank(),
                               axis.ticks = element_blank())
    ))


ggplots$ggs[[1]] + ggtitle("team 1")
ggsave("./publication/figures/figure3/sc2_UMAP_independent_team1.pdf", height = 4, width = 4)
ggsave("./publication/figures/figure3/sc2_UMAP_independent_team1.jpeg", height = 4, width = 4)

ggplots$ggs[[2]]+ ggtitle("team 2")
ggsave("./publication/figures/figure3/sc2_UMAP_independent_team2.pdf", height = 4, width = 4)
ggsave("./publication/figures/figure3/sc2_UMAP_independent_team2.jpeg", height = 4, width = 4)

ggplots$ggs[[3]]+ ggtitle("team 3")
ggsave("./publication/figures/figure3/sc2_UMAP_independent_team3.pdf", height = 4, width = 4)
ggsave("./publication/figures/figure3/sc2_UMAP_independent_team3.jpeg", height = 4, width = 4)

ggplots$ggs[[4]]+ ggtitle("team 4")
ggsave("./publication/figures/figure3/sc2_UMAP_independent_team4.pdf", height = 4, width = 4)
ggsave("./publication/figures/figure3/sc2_UMAP_independent_team4.jpeg", height = 4, width = 4)







#### Data together with best 1,4 teams ------------------------------

set.seed(1938)
measured_data_plus_team1_4 <- prediction_data %>% slice(1,4) %>%
    unnest(predictions) %>%
    select(submitterId,required_columns) %>% 
    bind_rows(
        gs %>% add_column(submitterId="data",.before = 1) %>% 
            group_by(cell_line,treatment,time) %>% 
            sample_n(.,min(1000,n()))
    ) %>%
    select(submitterId,required_columns)


um_res_data_team1_4 <- uwot::umap(measured_data_plus_team1_4[,-1:-4], 
                                  metric="euclidean", 
                                  init="pca", 
                                  scale = "none",
                                  n_components = 2,
                                  fast_sgd=TRUE,
                                  min_dist=0.1,
                                  n_neighbors = 15,
                                  verbose = T)

preped_data <- measured_data_plus_team1_4 %>% 
    bind_cols(tibble(umap.x = um_res_data_team1_4[,1],
                     umap.y = um_res_data_team1_4[,2])) %>% 
    sample_frac(0.1) %>% 
    left_join(hide_names,by = "submitterId") 



preped_data %>%
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line), size=.2) +
    #geom_density2d(aes(umap.x,umap.y,col=cell_line),size = 0.5) +
    # coord_equal() +
    guides(colour = "none") +
    labs(x="",y="") +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())  +
    facet_wrap(~submitterId)



create_plot <- function(preped_data){
    preped_data %>%
        ggplot() + 
        geom_point(aes(umap.x,umap.y,col=cell_line), size=.2) +
        #geom_density2d(aes(umap.x,umap.y,col=cell_line),size = 0.5) +
        # coord_equal() +
        guides(colour = "none") +
        labs(x="",y="") +
        theme_bw() + 
        theme(panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank()) 
}
preped_data%>% filter(submitterId == "data") %>% create_plot()
ggsave("./publication/figures/figure3/sc2_UMAP_data_teams1_4_data.pdf", height = 4, width = 4)
ggsave("./publication/figures/figure3/sc2_UMAP_data_teams1_4_data.jpeg", height = 4, width = 4)

preped_data%>% filter(submitterId == "icx_bxai") %>% create_plot()
ggsave("./publication/figures/figure3/sc2_UMAP_data_teams1_4_team1.pdf", height = 4, width = 4)
ggsave("./publication/figures/figure3/sc2_UMAP_data_teams1_4_team1.jpeg", height = 4, width = 4)

preped_data%>% filter(submitterId == "PaL") %>% create_plot()
ggsave("./publication/figures/figure3/sc2_UMAP_data_teams1_4_team4.pdf", height = 4, width = 4)
ggsave("./publication/figures/figure3/sc2_UMAP_data_teams1_4_team4.jpeg", height = 4, width = 4)


#### LAyout by team and then add points --------

measured_data <- gs %>% add_column(submitterId="data",.before = 1) %>% 
    group_by(cell_line,treatment,time) %>% 
    sample_n(.,min(1000,n())) %>% select(required_columns)


umap_data <-  uwot::umap(measured_data[,-1:-3], 
                            metric="euclidean", 
                            ret_model = TRUE,
                            init="pca", 
                            scale = "none",
                            n_components = 2,
                            min_dist=0.01,
                            n_neighbors = 15,
                            verbose = T)

prediction_1_4 <- prediction_data %>% slice(1,4) %>%
    unnest(predictions) %>%
    select(submitterId,required_columns) 

mnist_predictions_umap <- umap_transform(prediction_1_4[,-1:-4], umap_data, verbose = TRUE)


# show the data only


measured_data %>% bind_cols(tibble(umap.x = umap_data$embedding[,1],
                                              umap.y = umap_data$embedding[,2])) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line), size=.2) +
    guides(colour = "none") +
    labs(x="",y="") +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) 



# show data and projected predictions:


measured_data %>% bind_cols(tibble(umap.x = umap_data$embedding[,1],
                                   umap.y = umap_data$embedding[,2])) %>% 
    mutate(submitterId="data") %>%
    bind_rows(bind_cols(prediction_1_4, tibble(umap.x = mnist_predictions_umap[,1],
                                 umap.y = mnist_predictions_umap[,2]))   ) %>%
    left_join(hide_names,by="submitterId") %>%
    mutate(alt_name = ifelse(alt_name!="data",paste0("team #",alt_name),'data')) %>%
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line), size=.2) +
    # guides(colour = "none") +
    labs(x="",y="") +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=1))) +
    facet_wrap(~alt_name)

ggsave("./publication/figures/figure3/sc2_UMAP_data_teams1_4.pdf", height = 4, width = 14)
ggsave("./publication/figures/figure3/sc2_UMAP_data_teams1_4.jpeg", height = 4, width = 14)



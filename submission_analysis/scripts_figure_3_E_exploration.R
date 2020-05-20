# SC2 UMAP of participants and real cells. 
# exploration of ideas


library(tidyverse)
library(uwot)
library(ggplot2)
# read golden standard
gs <- read_csv("./challenge_data/validation_data/sc2gold.csv") 

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
        sample_n(.,100)
}



prediction_data <- SC_leaderboard  %>%
    mutate(predictions = map(file.path(submission_folder,submissions),read_sample))



#### Data alone: 
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
                          min_dist=0.01,
                          n_neighbors = 5,
                          verbose = T)

measured_data %>% bind_cols(tibble(umap.x = um_res_data[,1],umap.y = um_res_data[,2])) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=treatment),size = 0.5) + facet_wrap(~cell_line)



measured_data %>% bind_cols(tibble(umap.x = um_res_data[,1],umap.y = um_res_data[,2])) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line),size = 0.5) +
    #geom_density2d(aes(umap.x,umap.y,col=cell_line),size = 0.5) +
    theme_bw()


### Data and best team togehter


measured_data_plus_team1 <- prediction_data %>% filter(submitterId=="icx_bxai") %>%
    unnest(predictions) %>%
    select(submitterId,required_columns) %>% 
    bind_rows(
        gs %>% add_column(submitterId="data",.before = 1) %>% 
            group_by(cell_line,treatment,time) %>% 
            sample_n(.,min(100,n()))
    ) %>%
    select(submitterId,required_columns)


um_res_data_team1 <- uwot::umap(measured_data_plus_team1[,-1:-4], 
                          metric="euclidean", 
                          init="pca", 
                          scale = "none",
                          n_components = 2,
                          fast_sgd=TRUE,
                          min_dist=0.01,
                          n_neighbors = 5,
                          verbose = T)

measured_data_plus_team1 %>% bind_cols(tibble(umap.x = um_res_data_team1[,1],umap.y = um_res_data_team1[,2])) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line,shape=submitterId)) +
    #geom_density2d(aes(umap.x,umap.y,col=cell_line),size = 0.5) +
    theme_bw() + facet_wrap(~submitterId)

measured_data_plus_team1 %>% bind_cols(tibble(umap.x = um_res_data_team1[,1],umap.y = um_res_data_team1[,2])) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=submitterId,shape=submitterId),size=0.3) +
    #geom_density2d(aes(umap.x,umap.y,col=cell_line),size = 0.5) +
    theme_bw() + facet_wrap(~cell_line)


measured_data_plus_team1 %>% bind_cols(tibble(umap.x = um_res_data_team1[,1],umap.y = um_res_data_team1[,2])) %>% 
    filter(submitterId=="data") %>%
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line),size=0.3) +
    #geom_density2d(aes(umap.x,umap.y,col=cell_line),size = 0.5) +
    theme_bw() 

### Data and icx_bxai + PAl (sampler) team togehter -----------------------


measured_data_plus_team1_4 <- prediction_data %>% filter(submitterId %in% c("PaL","icx_bxai")) %>%
    unnest(predictions) %>%
    select(submitterId,required_columns) %>% 
    bind_rows(
        gs %>% add_column(submitterId="data",.before = 1) %>% 
            group_by(cell_line,treatment,time) %>% 
            sample_n(.,min(100,n()))
    ) %>%
    select(submitterId,required_columns)


um_res_data_team1_4 <- uwot::umap(measured_data_plus_team1_4[,-1:-4], 
                                metric="euclidean", 
                                init="pca", 
                                scale = "none",
                                n_components = 2,
                                fast_sgd=TRUE,
                                min_dist=0.01,
                                n_neighbors = 5,
                                verbose = T)

measured_data_plus_team1_4 %>% bind_cols(tibble(umap.x = um_res_data_team1_4[,1],umap.y = um_res_data_team1_4[,2])) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line,shape=submitterId)) +
    #geom_density2d(aes(umap.x,umap.y,col=cell_line),size = 0.5) +
    theme_bw() + facet_wrap(~submitterId)




### Best team alone     ------------


team1 <- prediction_data %>% filter(submitterId=="icx_bxai") %>%
    unnest(predictions) %>%
    select(submitterId,required_columns) 


um_res_team1 <- uwot::umap(team1[,-1:-4], 
                                metric="euclidean", 
                                init="pca", 
                                scale = "none",
                                n_components = 2,
                                fast_sgd=TRUE,
                                min_dist=0.01,
                                n_neighbors = 5,
                                verbose = T)

team1 %>% bind_cols(tibble(umap.x = um_res_team1[,1],umap.y = um_res_team1[,2])) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line,shape=submitterId)) +
    #geom_density2d(aes(umap.x,umap.y,col=cell_line),size = 0.5) +
    theme_bw() + facet_wrap(~submitterId)



#### Sampling team (team 4)
team4 <- prediction_data %>% filter(submitterId=="PaL") %>%
    unnest(predictions) %>%
    select(submitterId,required_columns) 


um_res_team4 <- uwot::umap(team4[,-1:-4], 
                           metric="euclidean", 
                           init="pca", 
                           scale = "none",
                           n_components = 2,
                           fast_sgd=TRUE,
                           min_dist=0.01,
                           n_neighbors = 5,
                           verbose = T)

team4 %>% bind_cols(tibble(umap.x = um_res_team4[,1],umap.y = um_res_team4[,2])) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line,shape=submitterId)) +
    #geom_density2d(aes(umap.x,umap.y,col=cell_line),size = 0.5) +
    theme_bw() + facet_wrap(~submitterId)


### Data and prediction together

all_data <- prediction_data %>% unnest(predictions) %>%
    select(submitterId,required_columns) %>% 
    bind_rows(
        gs %>% add_column(submitterId="data",.before = 1) %>% 
            group_by(cell_line,treatment,time) %>% 
            sample_n(.,min(100,n()))
    ) %>%
    select(submitterId,required_columns)



um_res <- uwot::umap(all_data[,-1:-4], 
                     metric="euclidean", 
                     init="pca", 
                     scale = "Z",
                     n_components = 2,
                     fast_sgd=TRUE,
                     min_dist=0.01,
                     n_neighbors = 5,
                     verbose = T)


all_data %>% bind_cols(tibble(umap.x = um_res[,1],umap.y = um_res[,2])) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=cell_line),size = 0.5) + facet_wrap(~submitterId)


all_data %>% bind_cols(tibble(umap.x = um_res[,1],umap.y = um_res[,2])) %>% 
    filter(! submitterId %in% c("KAUST_RSS","X.Huiyuan","Raghava_India_SCS")) %>% 
    ggplot() + 
    geom_point(aes(umap.x,umap.y,col=submitterId),size = 0.5) + facet_wrap(~cell_line)


# per cell-line 
# run umap per cell-line
all_data_umap <- all_data %>%  group_by(cell_line) %>% nest() %>%
    mutate(umap_CL = map(data, ~ uwot::umap(.x[,-1:-3], 
                                                      metric="euclidean", 
                                                      init="pca", 
                                                      scale = "Z",
                                                      n_components = 2,
                                                      fast_sgd=TRUE,
                                                      min_dist=0.01,
                                                      n_neighbors = 5,
                                                      verbose = T)))



# bind the umap coordinates to single cell data: 
umap_coords <- all_data_umap %>% 
    mutate(tmp = map2(.x = data,.y = umap_CL, ~bind_cols(.x,tibble(umap.x = .y[,1],umap.y = .y[,2])) ))
    

ggplots <- umap_coords %>% mutate(ggs = map(tmp, ~ ggplot(data=.x) +
                                                geom_point(aes(umap.x,umap.y),size = 0.5) +
                                                facet_wrap(~submitterId)  ))

        
ggplots$ggs[[4]]

ggplots <- umap_coords %>% mutate(ggs = map(tmp, ~ ggplot(data=.x) +
                                                geom_density_2d(aes(umap.x,umap.y),size = 0.5) +
                                                facet_wrap(~submitterId)  ))
ggplots$ggs[[4]]

# DREAM challenge post analysis
# boxplots cov, mean
# SC2: compare the covariance and mean between data and predictions
# A. Gabor


library(tidyverse)
# SC2

combined_statistics <- read_rds("./submission_analysis/intermediate_data/sc2_stats_conditions.rds")
sserr <- read_rds("./submission_analysis/intermediate_data/sc2_stats_sumSquared_conditions.rds")
source("./scoring_scripts/score_sc2.R")




stats_by_var <- combined_statistics %>% gather(teams,prediction,7:22) %>%
    mutate(prediction_error = abs(standard-prediction)) %>%
    group_by(cell_line,treatment,time,stat_variable) %>%
    summarise(median_err  = median(prediction_error),
              median_prediction = median(prediction),
              top75_prediction = quantile(prediction,.75),
              top25_prediction = quantile(prediction,.25),
              measurement_stat = median(standard)) %>%
    arrange(desc(median_err)) 


# find some interesting variabels: the strongest
median_stat_var_magnitude <- combined_statistics %>%
    filter(!grepl("cov_",stat_variable)) %>%
    group_by(stat_variable) %>% 
    summarise(median_stat = median(standard)) %>%
    mutate(stat_variable = gsub("mean_","",stat_variable)) %>%
    arrange(desc(median_stat))

cov_stat_var_magnitude <- combined_statistics %>% 
    filter(grepl("cov_",stat_variable)) %>%
    group_by(stat_variable) %>% 
    summarise(median_stat = median(standard)) %>%
    arrange(desc(median_stat))

# pairs withaout variance (i.e. sd^2): 
cov_X_Y = cov_stat_var_magnitude %>% 
    separate(stat_variable,into = c(NA,"var1","var2"),sep = "_",remove = FALSE) %>%
    filter(var1!=var2) %>%
    arrange(desc(median_stat))

cov_stat_var_magnitude <- cov_stat_var_magnitude %>%
    mutate(stat_variable = gsub("cov_","",stat_variable)) %>%
    mutate(stat_variable = gsub("_"," - ",stat_variable)) %>%
    arrange(desc(median_stat))

# plot mean estimation
combined_statistics %>% gather(source,stat_value,6:13) %>%
    mutate(source = ifelse(source=="standard","measurement","prediction")) %>% 
    mutate(stat_variable = gsub("mean_","",stat_variable)) %>%
    filter(stat_variable %in% median_stat_var_magnitude$stat_variable) %>%
    mutate(stat_variable =factor(stat_variable,levels = median_stat_var_magnitude$stat_variable)) %>%
    ggplot() + 
    geom_boxplot(aes(stat_variable,stat_value,fill=source,color=source),outlier.size = .4) + 
    theme_bw() +
    xlab("") + ylab("median signal intensity") + 
    theme(axis.text.x = element_text(angle = 90,hjust = 1),
          legend.position = c(0.8,0.8),
          legend.title = element_blank(),
          panel.grid = element_blank())


ggsave("./publication/figures/figure3/sc2_mean_accuracy.pdf",width = 6,height = 4)

# plot variance estimation

combined_statistics %>% gather(source,stat_value,6:13) %>%
    mutate(source = ifelse(source=="standard","measurement","prediction")) %>% 
    filter(stat_variable %in% cov_X_Y$stat_variable[c(1:20,575:595)]) %>%
    mutate(stat_variable = gsub("cov_","",stat_variable)) %>%
    mutate(stat_variable = gsub("_"," - ",stat_variable)) %>% 
    mutate(stat_variable = factor(stat_variable,levels = cov_stat_var_magnitude$stat_variable)) %>%
    ggplot() + 
    geom_boxplot(aes(stat_variable,stat_value,fill=source,color=source)) + 
    theme(axis.text.x = element_text(angle = 90,hjust = 1),
          panel.grid = element_blank()) +
    xlab("") + ylab("Covariances") + guides(col="none",fill="none") 
# theme(legend.position = c(0.8,0.8),
#       legend.title = element_blank())

ggsave("./publication/figures/figure3/sc2_cov_accuracy.pdf",width = 6,height = 5)



### Checking for intersting stats, that change across conditions:
# 1. Build Anova model using CL, TR, Time -- unfortunately it is not possible for treatment

aov_ <- combined_statistics %>% 
    select(cell_line, treatment, time, stat_variable, standard) %>%
    group_by(stat_variable) %>% nest() %>%
    mutate(anova = map(data,function(d){
        aov_res <- aov(standard ~ cell_line + treatment + time,data = d)
        summary(aov_res)
    }))



combined_statistics %>% filter(grepl("cov",stat_variable)) %>%
    separate(stat_variable,into = c("stat","var1","var2"),sep = "_", remove = FALSE) %>%
    filter(var1 != var2) %>% group_by(stat_variable) %>%
    summarise(max_change = max(standard) - min(standard)) %>% arrange(desc(max_change))


combined_statistics %>% filter(stat_variable == "cov_p.ERK_p.MEK" ) %>% 
    ggplot()  +  geom_boxplot(aes(cell_line,standard))+ geom_point(aes(cell_line,standard,col=treatment)) +
    geom_boxplot(aes(cell_line,icx_bxai))+ geom_point(aes(cell_line,icx_bxai))

combined_statistics %>% filter(stat_variable == "cov_IdU_p.RB" ) %>% 
    ggplot(aes(cell_line,icx_bxai))  +  geom_boxplot()+ geom_point()


combined_statistics %>% filter(stat_variable == "cov_IdU_p.RB" ) %>% 
    ggplot(aes(time,standard)) +  geom_boxplot() + geom_point() 

combined_statistics %>% filter(stat_variable == "cov_p.ERK_p.MEK" ) %>% 
    ggplot(aes(treatment,standard)) +  geom_boxplot() + geom_point() 



####
# we would like to show how good is their estimates. 
# 1. which covariances are significant? 
# for this we generate a null-distribution by independently shuffling the nodes of single cells. 


# we shuffle the nodes to generate a null distribution for the covariance matrix.
# random Cov mat -------------------------------------------------------------
if(FALSE){
    gs <- read_csv("./challenge_data/validation_data/sc2gold.csv")
    
    nodes = names(gs)[-1:-5]
    tmp = gs %>% mutate_at(.vars = nodes, .funs = sample) %>% 
        select(-cellID,-fileID) %>% 
        data_to_stats()
    
    # flatten the upper triangular of a symmetric matrix with diagonal to a table (from Stackoverflow)
    flattenCovMatrix <- function(covmat) {
        ut <- upper.tri(covmat,diag = TRUE)
        tibble(
            stat_variable = paste0("cov_", rownames(covmat)[row(covmat)[ut]],"_",rownames(covmat)[col(covmat)[ut]]),
            stat_value  =(covmat)[ut]
        )
    }
    
    generate_rnd_cov <- function(i){
        pb$tick()
        gs %>% 
            select(-cellID,-fileID) %>% 
            group_by(cell_line,treatment,time) %>%
            nest(.key = "data") %>%
            mutate(data = map(data,~mutate_at(.tbl = .x,nodes,sample))) %>%
            mutate(cov_values = map(data,cov)) %>% 
            mutate(vec_cov = map(cov_values,flattenCovMatrix)) %>% 
            unnest(vec_cov)
    }
    
    pb <- progress::progress_bar$new(total = 100,format = "  downloading [:bar] :percent eta: :eta",)
    
    rnd_covs <- tibble(rep = 1:100) %>% rowwise() %>%
        mutate(rnd_cov = map(rep, generate_rnd_cov))
    
    pb <- progress::progress_bar$new(total = length(101:1000),format = " [:bar] :percent eta: :eta",clear = FALSE)
    rnd_covs_2 <- tibble(rep = 101:1000) %>% rowwise() %>%
        mutate(rnd_cov = map(rep, generate_rnd_cov))
    
    
    rnd_covs_all <- rnd_covs_2 %>% unnest(rnd_cov)
    rnd_covs_all <- bind_rows(rnd_covs_all,rnd_covs %>% unnest(rnd_cov))
    
    
    write_rds(rnd_covs_all,"./submission_analysis/intermediate_data/sc2_random_cov_stats.rds")
}else{
    rnd_covs_all = read_rds("./submission_analysis/intermediate_data/sc2_random_cov_stats.rds")
}
# end of random Cov mat -------------------------------------------------------------


# combine the real stats and the random to compute the pvalue:
all_stats <- rnd_covs_all %>% left_join(combined_statistics,by = c("cell_line", "treatment", "time", "stat_variable"))

alpha = 0.05

exceeding_rate <- all_stats %>% group_by(cell_line, treatment,time, stat_variable) %>%
    summarise(exc = ifelse(standard[[1]] > 0, sum(stat_value>standard)/1000,  # testing how many time Rnd is larger
                           sum(stat_value<standard)/1000)) # testing how many time Rnd is smaller


exceeding_rate[which(exceeding_rate$exc > alpha/2),] %>% pull(stat_variable) %>% table() %>% sort()

exceeding_rate[which(exceeding_rate$exc > alpha/2),] %>% pull(exc) %>% hist(.)



all_stats %>% separate(stat_variable,into = c("stat","var1","var2"),sep = "_", remove = FALSE) %>%
    filter(var1 !=var2)


exceeding_rate
# the problem is that the null distribution is really narrow. Many covariance is 
# significant, but the efect size is small. 


# check the score of random covariances :
# hypothesis, it is similar to the cov score we found for the combination of 
# participants

all_stats %>% group_by(rep) %>%
    filter(grepl("cov",stat_variable)) %>%
    summarise(squared_cov_score = sum((standard-stat_value)^2))



##### ANOVA across cell-lines -------------------------------------------------
# Run Anova to see, in which cases the cell-line makes a strong difference between covariances

aov_cell_line <- combined_statistics %>% 
    select(cell_line, treatment, time, stat_variable, standard) %>%
    group_by(stat_variable) %>% nest() %>%
    mutate(anova = map(data,function(d){
        
        aov_res <- aov(standard ~ cell_line,data = d)
        aov_res %>% broom::tidy()
    })) %>%  unnest(anova)

# see extremes (small adn large)
aov_cell_line %>% filter(term=="cell_line") %>% arrange(desc(statistic))
aov_cell_line %>% filter(term=="cell_line") %>% arrange(desc(p.value))


combined_statistics %>% filter(stat_variable =="cov_IdU_Ki.67") %>%
    ggplot() + geom_point(aes(cell_line,standard, col=treatment)) +
    geom_point(aes(cell_line,icx_bxai))  + 
    theme(axis.text.x = element_text(angle = 90))

# Let's check covariance for ERK
aov_cell_line %>% filter(term=="cell_line") %>% filter(grepl("ERK",stat_variable)) %>% arrange(desc(statistic))


# overall IdU and p-RB are the strongest, probably cell-cycle effect. 
cell_line_levels = c( "184B5","HCC202", "ZR751", # iMEK
                      "BT483","MCF12A","MDAMB468", # iEGFR
                      "MDAMB231", "SKBR3", "UACC3199",# iPI3K
                      "HCC1428", "HCC1806", "Hs578T") # iPKC

combined_statistics %>% filter(stat_variable =="cov_IdU_p.RB") %>%
    select(cell_line, treatment, time, stat_variable,standard,icx_bxai) %>%
    mutate(cell_line = factor(cell_line, levels = cell_line_levels)) %>%
    gather(source,stat_value,standard,icx_bxai) %>%
    ggplot() + 
    geom_boxplot(aes(cell_line,stat_value, fill=source),outlier.size = 1, size = 0.1) +
    theme(axis.text.x = element_text(angle = 90),
          panel.grid = element_blank(), legend.title = element_blank(), legend.position = c(0.8,0.2)) +
    scale_color_manual(values = c("standard" = RColorBrewer::brewer.pal(7,"Dark2")[1],
                                  "icx_bxai" = RColorBrewer::brewer.pal(7,"Dark2")[2])) +
    xlab("cell line") + 
    ylab("covariance(IdU, p-RB)") 
ggsave("./publication/figures/figure3/sc2_cov_example_IdU_RB.pdf",width = 3,height = 3)


# export boxplots showsn differenes between cell-lines -------------------------
cell_line_levels = c( "184B5","HCC202", "ZR751", # iMEK
                      "BT483","MCF12A","MDAMB468", # iEGFR
                      "MDAMB231", "SKBR3", "UACC3199",# iPI3K
                      "HCC1428", "HCC1806", "Hs578T") # iPKC

plot_cov_per_cellline <- function(data,target_dir,var){
    
    var1 = strsplit(var,"_")[[1]][[2]]
    var2 = strsplit(var,"_")[[1]][[3]]
    
    gg <- ggplot(data) + 
        geom_boxplot(aes(cell_line,stat_value, fill=source),outlier.size = 1, size = 0.1) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90),
              panel.grid = element_blank(), 
              legend.title = element_blank()) +
        scale_color_manual(values = c("standard" = RColorBrewer::brewer.pal(7,"Dark2")[1],
                                      "icx_bxai" = RColorBrewer::brewer.pal(7,"Dark2")[2])) +
        xlab("cell line") + 
        ylab(paste0( "covariance(",var1,", ", var2,")" ) ) 
    ggsave(filename = file.path(target_dir,paste0(var,".pdf")), plot = gg, width = 4, height = 3)
}

combined_statistics %>% group_by(stat_variable) %>% 
    filter(grepl("cov_",stat_variable)) %>% 
    select(cell_line, treatment, time, stat_variable,standard,icx_bxai) %>%
    mutate(cell_line = factor(cell_line, levels = cell_line_levels)) %>%
    gather(source,stat_value,standard,icx_bxai) %>%
    mutate(source=ifelse(source=="icx_bxai","team #1","measured")) %>%
    group_by(stat_variable) %>% 
    group_walk(.tbl = ., 
               ~ plot_cov_per_cellline(data = .,
               target_dir = "./submission_analysis/figures/sc2_cov_per_cellline/",
               var=.y$stat_variable)
    )

##### ANOVA across treatments --------------------------------------------------
# Run Anova to see, in which cases the cell-line makes a strong difference between covariances

aov_treatment <- combined_statistics %>% 
    select(cell_line, treatment, time, stat_variable, standard) %>%
    group_by(stat_variable) %>% nest() %>%
    mutate(anova = map(data,function(d){
        
        aov_res <- aov(standard ~ treatment,data = d)
        aov_res %>% broom::tidy()
    })) %>%  unnest(anova)

# see extremes (small adn large)
aov_treatment %>% filter(term=="treatment") %>% arrange(desc(statistic))
aov_treatment %>% filter(term=="treatment") %>% arrange(desc(p.value))


# export boxplots showsn between treatments------------------------------------
cell_line_levels = c( "184B5","HCC202", "ZR751", # iMEK
                      "BT483","MCF12A","MDAMB468", # iEGFR
                      "MDAMB231", "SKBR3", "UACC3199",# iPI3K
                      "HCC1428", "HCC1806", "Hs578T") # iPKC

plot_cov_per_treatment <- function(data,target_dir,var){
    
    var1 = strsplit(var,"_")[[1]][[2]]
    var2 = strsplit(var,"_")[[1]][[3]]
    
    gg <- ggplot(data) + 
        geom_boxplot(aes(treatment,stat_value, fill=source),outlier.size = 1, size = 0.1) +
        theme(axis.text.x = element_text(angle = 90),
              panel.grid = element_blank(),
              legend.title = element_blank()
              ) +
        scale_color_manual(values = c("measured" = RColorBrewer::brewer.pal(7,"Dark2")[1],
                                      "team #1" = RColorBrewer::brewer.pal(7,"Dark2")[2])) +
        xlab("treatment") + 
        ylab(paste0( "covariance(",var1,", ", var2,")" ) )
    ggsave(filename = file.path(target_dir,paste0(var,".pdf")), plot = gg, width = 4, height = 3)
}

combined_statistics %>% group_by(stat_variable) %>% 
    filter(grepl("cov_",stat_variable)) %>% 
    select(cell_line, treatment, time, stat_variable,standard,icx_bxai) %>%
    mutate(cell_line = factor(cell_line, levels = cell_line_levels)) %>%
    gather(source,stat_value,standard,icx_bxai) %>%
    mutate(source=ifelse(source=="icx_bxai","team #1","measured")) %>%
    group_by(stat_variable) %>% 
    group_walk(.tbl = ., 
               ~ plot_cov_per_treatment(data = .,
                                       target_dir = "./submission_analysis/figures/sc2_cov_per_treatment/",
                                       var=.y$stat_variable)
    )

# export median changes between treatments  -----------------------------
# cell_line_levels = c( "184B5","HCC202", "ZR751", # iMEK
#                       "BT483","MCF12A","MDAMB468", # iEGFR
#                       "MDAMB231", "SKBR3", "UACC3199",# iPI3K
#                       "HCC1428", "HCC1806", "Hs578T") # iPKC
# 
# 
# 
# 
# 
# 
# plot_median_per_treatment <- function(data,target_dir,var){
#     
#     var1 = strsplit(var,"_")[[1]][[2]]
#     var2 = strsplit(var,"_")[[1]][[3]]
#     
#     gg <- ggplot(data) + 
#         geom_boxplot(aes(treatment,stat_value, fill=source),outlier.size = 1, size = 0.1) +
#         theme(axis.text.x = element_text(angle = 90),
#               panel.grid = element_blank(),
#               legend.title = element_blank()
#         ) +
#         scale_color_manual(values = c("measured" = RColorBrewer::brewer.pal(7,"Dark2")[1],
#                                       "team #1" = RColorBrewer::brewer.pal(7,"Dark2")[2])) +
#         xlab("treatment") + 
#         ylab(paste0( "covariance(",var1,", ", var2,")" ) )
#     ggsave(filename = file.path(target_dir,paste0(var,".pdf")), plot = gg, width = 4, height = 3)
# }
# 
# combined_statistics %>% group_by(stat_variable) %>% 
#     filter(grepl("mean_",stat_variable)) %>% 
#     select(cell_line, treatment, time, stat_variable,standard,icx_bxai) %>%
#     mutate(cell_line = factor(cell_line, levels = cell_line_levels)) %>%
#     gather(source,stat_value,standard,icx_bxai) %>%
#     mutate(source=ifelse(source=="icx_bxai","team #1","measured")) %>%
#     group_by(stat_variable) %>% 
#     group_walk(.tbl = ., 
#                ~ plot_cov_per_treatment(data = .,
#                                         target_dir = "./submission_analysis/figures/sc2_cov_per_treatment/",
#                                         var=.y$stat_variable)
#     )
# 
# 
# tmp %>% filter(stat_variable == "mean_p.ERK") %>% 
#     ggplot() + 
#     geom_point(aes(time,stat_value, color=source)) +
#     geom_line(aes(time,stat_value, color=source)) +
#     theme(axis.text.x = element_text(angle = 90),
#           panel.grid = element_blank(),
#           legend.title = element_blank()
#     ) +
#     scale_color_manual(values = c("measured" = RColorBrewer::brewer.pal(7,"Dark2")[1],
#                                   "team #1" = RColorBrewer::brewer.pal(7,"Dark2")[2])) +
#     xlab("treatment") + 
#     facet_grid(treatment~cell_line)



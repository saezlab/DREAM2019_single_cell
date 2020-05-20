# DREAM challenge post analysis
# boxplots cov, mean
# SC3: compare the covariance and mean between data and predictions
# A. Gabor


library(tidyverse)
# SC3

combined_statistics <- read_rds("./submission_analysis/intermediate_data/sc3_stats_conditions.rds")
sserr <- read_rds("./submission_analysis/intermediate_data/sc3_stats_sumSquared_conditions.rds")
source("./scoring_scripts/score_sc2.R")



stats_by_var <- combined_statistics %>% gather(teams,prediction,7:20) %>%
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
# show only the top 7 teams
team_cutoff = 7
combined_statistics %>% gather(source,stat_value,6+0:team_cutoff) %>%
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


ggsave("./publication/figures/figure4/sc3_mean_accuracy.pdf",width = 6,height = 4)

# plot variance estimation

combined_statistics %>% gather(source,stat_value,6+0:team_cutoff) %>%
    mutate(source = ifelse(source=="standard","measurement","prediction")) %>% 
    filter(stat_variable %in% cov_X_Y$stat_variable[c(1:20,575:595)]) %>%
    mutate(stat_variable = gsub("cov_","",stat_variable)) %>%
    mutate(stat_variable = gsub("_"," - ",stat_variable)) %>% 
    mutate(stat_variable = factor(stat_variable,levels = cov_stat_var_magnitude$stat_variable)) %>%
    ggplot() + 
    geom_boxplot(aes(stat_variable,stat_value,fill=source,color=source)) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,hjust = 1),
          panel.grid = element_blank()) +
    xlab("") + ylab("Covariances") + guides(col="none",fill="none") 
# theme(legend.position = c(0.8,0.8),
#       legend.title = element_blank())

ggsave("./publication/figures/figure4/sc3_cov_accuracy.pdf",width = 6,height = 4.8)



### Checking for intersting stats, that change across conditions:
# 1. Build Anova model using CL, TR, Time -- unfortunately it is not possible for treatment

aov_ <- combined_statistics %>% 
    select(cell_line, treatment, time, stat_variable, standard) %>%
    group_by(stat_variable) %>% nest() %>%
    mutate(anova = map(data,function(d){
        aov_res <- aov(standard ~ cell_line + time,data = d)
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
    ggplot(aes(factor(time),standard)) +  geom_boxplot() + geom_point() 



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


combined_statistics %>% filter(stat_variable =="cov_IdU_p.RB") %>%
    select(cell_line, treatment, time, stat_variable,standard,icx_bxai) %>%
    gather(source,stat_value,standard,icx_bxai) %>%
    ggplot() + 
    geom_boxplot(aes(cell_line,stat_value, fill=source),outlier.size = 1, size = 0.1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          panel.grid = element_blank(), legend.title = element_blank(), legend.position = c(0.8,0.2)) +
    scale_color_manual(values = c("standard" = RColorBrewer::brewer.pal(7,"Dark2")[1],
                                  "icx_bxai" = RColorBrewer::brewer.pal(7,"Dark2")[2])) +
    xlab("cell line") + 
    ylab("covariance(IdU, p-RB)") 
ggsave("./publication/figures/figure4/sc3_cov_example_IdU_RB.pdf",width = 5,height = 3)


# export boxplots showsn differenes between cell-lines -------------------------

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
               target_dir = "./submission_analysis/figures/sc3_cov_per_cellline/",
               var=.y$stat_variable)
    )



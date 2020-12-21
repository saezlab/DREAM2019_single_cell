
library(tidyverse)
library(cowplot)

submission_folder <- "./submission_data/final_round/SC1"
rmse_data <- read_rds("./submission_analysis/intermediate_data/sc1_rmse_conditions.rds")
median_data <- read_rds("./submission_analysis/intermediate_data/sc1_median_conditions.rds")

SC_leaderboard = read_csv(
    file.path(submission_folder,"leaderboard_final_sc1.csv")) %>%
    select(-writeUp, -createdOn) %>% 
    mutate(submissions = paste0(objectId,".csv")) %>% arrange(score) %>% 
    mutate(submitterId = make.names(submitterId))


# General heatmap: 
# conditional RMSE for each team 
rmse_data %>% select(-cell_line,-treatment,-time,-marker) %>% 
    pheatmap::pheatmap(mat = .,
                       color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name ="YlOrRd"))(100),
                       show_rownames = FALSE)
# suggest to temove X.Huiyuan's and SCG's predictions
rmse_data %>% select(-X.Huiyuan,-SCG) %>%
    arrange(marker,cell_line,time,treatment) %>%
    select(-cell_line,-treatment,-time,-marker) %>% 
    pheatmap::pheatmap(mat = .,
                       color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name ="YlOrRd"))(10),
                       cluster_rows = FALSE)




# supp material 2:
# barplot that shows how many conditions were won by each team. 
#select(-rmse) %>% 
SC_leaderboard$submitterId = gsub("X.","",SC_leaderboard$submitterId)
    

rmse_data %>% 
    gather(team, rmse,-1:-4) %>%  
    mutate(team = gsub("X.","",team)) %>%
    group_by(cell_line, treatment, time, marker) %>%
    top_n(-1,wt = rmse) %>%
    ungroup() %>% 
    pull(team) %>% table() %>% enframe("team","cond_winner") %>%
    arrange(desc(cond_winner)) %>%
    mutate(team = factor(team, levels = SC_leaderboard$submitterId)) %>% 
    ggplot() + geom_col(aes(team, cond_winner,fill=team)) + theme_bw()+
    theme(axis.text.x = element_text(hjust=1,angle = 45)) + 
    ylab("conditional best") + guides(fill = "none") + xlab("ranked teams")

ggsave("./publication/figures/supp_figures/sc1_winners_by_conditions_bar.pdf",width = 5,height = 4.7)	


rmse_data %>% gather(team, rmse,-1:-4) %>% group_by(cell_line, treatment, time, marker) %>%
    top_n(-1,wt = rmse) %>%
    ungroup() %>% pull(team) %>% table() %>% sort()

#################
# what is the most difficult part to predict? 
# plot the RMSE from  different views

# Is there differences between cell-lines? (how difficult to predict them?)
# 

## plot tile 
# 1. for each condition we rank the participants (1:22)
# 2. we compute the 
RColorBrewer::display.brewer.all()
custom_colors = colorRampPalette(RColorBrewer::brewer.pal(9,"Greens"))

rmse_data %>% gather(team, rmse,-1:-4) %>%
    # fix time for plotting
    mutate(time = ifelse(time==18,17,time)) %>%
    mutate(time = ifelse(time==14,13,time)) %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(med_rmse = median(rmse),
              mad_rmse = mad(rmse)) %>% 
    #gather(grp,value,med_rmse,mad_rmse) %>% 
    unite(cond1, cell_line, marker,sep = " ",remove = F) %>% 
    arrange(cell_line,marker) %>% 
    mutate(cond1 = factor(cond1,levels = unique(cond1))) %>% 
    unite(cond2, treatment, time,sep = " ", remove = F) %>%
    arrange(treatment, time) %>% 
    mutate(cond2 = factor(cond2,levels = unique(cond2))) %>% 
    
    ggplot(aes(cond1,cond2)) +
    #geom_tile(aes(fill=med_rmse,width = mad_rmse/max(mad_rmse),height = mad_rmse/max(mad_rmse)),na.rm = T) +
    #geom_tile(aes(fill=mad_rmse,width = med_rmse/max(med_rmse),height = med_rmse/max(med_rmse)),na.rm = T) +
    geom_point(aes(fill=med_rmse, size =mad_rmse/max(mad_rmse)),na.rm = T,shape=21) +
    theme_bw() +
    theme(axis.text.x = element_text(hjust=0,angle = 90)) +
    xlab("") + ylab("") + 
    scale_fill_gradientn(colors = custom_colors(100)) +
    scale_size_area(max_size = 5)
#scale_color_gradientn(colors = custom_colors(100))


# no time
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

rmse_data %>% gather(team, rmse,-1:-4) %>%
    # fix time for plotting
    mutate(time = ifelse(time==18,17,time)) %>%
    mutate(time = ifelse(time==14,13,time)) %>%
    group_by(cell_line, treatment, marker) %>%
    summarise(med_rmse = median(rmse),
              QCV_rmse = (quantile(rmse,0.75)-quantile(rmse,0.25))/quantile(rmse,0.5)) %>% 
    #gather(grp,value,med_rmse,mad_rmse) %>% 
    
    ggplot(aes(marker,treatment)) +
    #geom_tile(aes(fill=med_rmse,width = mad_rmse/max(mad_rmse),height = mad_rmse/max(mad_rmse)),na.rm = T) +
    #geom_tile(aes(fill=mad_rmse,width = med_rmse/max(med_rmse),height = med_rmse/max(med_rmse)),na.rm = T) +
    geom_point(aes(size=med_rmse, fill =QCV_rmse),na.rm = T,shape=21) +
    theme_bw() +
    theme(axis.text.x = element_text(hjust=0,angle = 90)) +
    xlab("") + ylab("") + 
    scale_fill_gradientn(colors = custom_colors(100)) +
    scale_size_area(max_size = 5)+
    facet_grid(~cell_line) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    coord_equal()
ggsave("./publication/figures/figure2/conditional_prediction_accuracy.pdf",width = 8.3,height = 7.4)

#scale_color_gradientn(colors = custom_colors(100))
scale()

# use corrplot


col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue","#00007F"))

col3 <- colorRampPalette(RColorBrewer::brewer.pal(8,"Reds"))

rmse_med <- rmse_data %>% gather(team, rmse,-1:-4) %>%
    # fix time for plotting
    mutate(time = ifelse(time==18,17,time)) %>%
    mutate(time = ifelse(time==14,13,time)) %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(med_rmse = median(rmse),
              mad_rmse = mad(rmse)) %>% 
    #gather(grp,value,med_rmse,mad_rmse) %>% 
    unite(cond1, cell_line, marker,sep = " ",remove = F) %>% 
    arrange(cell_line,marker) %>% 
    mutate(cond1 = factor(cond1,levels = unique(cond1))) %>% 
    unite(cond2, treatment, time,sep = " ", remove = F) %>%
    arrange(treatment, time) %>% 
    mutate(cond2 = factor(cond2,levels = unique(cond2))) %>%
    ungroup() %>%
    select(cond1,cond2,med_rmse) %>%
    spread(cond2,med_rmse) %>% column_to_rownames("cond1") %>% 
    as.matrix() %>% t()

rmse_mad <- rmse_data %>% gather(team, rmse,-1:-4) %>%
    # fix time for plotting
    mutate(time = ifelse(time==18,17,time)) %>%
    mutate(time = ifelse(time==14,13,time)) %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(med_rmse = median(rmse),
              mad_rmse = mad(rmse)) %>% 
    #gather(grp,value,med_rmse,mad_rmse) %>% 
    unite(cond1, cell_line, marker,sep = " ",remove = F) %>% 
    arrange(cell_line,marker) %>% 
    mutate(cond1 = factor(cond1,levels = unique(cond1))) %>% 
    unite(cond2, treatment, time,sep = " ", remove = F) %>%
    arrange(treatment, time) %>% 
    mutate(cond2 = factor(cond2,levels = unique(cond2))) %>%
    ungroup() %>%
    select(cond1,cond2,mad_rmse) %>%
    spread(cond2,mad_rmse) %>% column_to_rownames("cond1") %>% 
    as.matrix() %>% t()


corrplot::corrplot(corr = rmse_med,
                   is.corr = FALSE,
                   col = col3(20),
                   text.col=1)







p1 <- rmse_data %>% select(cell_line, treatment,  time, marker, icx_bxai) %>%
    mutate(y_id = paste(treatment,time),
           x_id = paste(marker,cell_line)) %>%
    arrange(marker,treatment,time) %>% ungroup() %>% select(-cell_line, -marker,-treatment,-time) %>%
    spread(x_id,icx_bxai) %>%
    column_to_rownames("y_id") %>% 
    pheatmap::pheatmap(mat = .,cluster_cols = F,cluster_rows = F,
                       main = "RMSE of icx_bxai model",
                       gaps_col = c(6,12,18,24,30), 
                       gaps_row = c(12,13,22,31,40), 
                       breaks = seq(0.4,2,length.out = 50),
                       color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(49))

ggsave("./submission_analysis/figures/sc1_condition_error/icx_bxai_heatmap.pdf",height = 8,width = 6,plot = p1$gtable)

p2 <- rmse_data %>% select(cell_line, treatment,  time, marker, Sleeping) %>%
    mutate(y_id = paste(treatment,time),
           x_id = paste(marker,cell_line)) %>%
    arrange(marker,treatment,time) %>% ungroup() %>% select(-cell_line, -marker,-treatment,-time) %>%
    spread(x_id,Sleeping) %>%
    column_to_rownames("y_id") %>% 
    pheatmap::pheatmap(mat = .,cluster_cols = F,cluster_rows = F,
                       main = "RMSE of Sleeping model",
                       gaps_col = c(6,12,18,24,30), 
                       gaps_row = c(12,13,22,31,40),
                       breaks = seq(0.4,2,length.out = 50),
                       color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(49))
ggsave("./submission_analysis/figures/sc1_condition_error/Sleeping_heatmap.pdf",height = 8,width = 6,plot = p2$gtable)

# to see the effect of each conditions, for each condition we select the best
# predictions among the teams. 
team_independent_predictions <- 	rmse_data %>%
    mutate(time = ifelse(time==18,17,time)) %>%
    mutate(time = ifelse(time==14,13,time)) %>%
    # select(-X.Huiyuan,-SCG) %>%
    gather(team,RMSE,-1:-4) %>%
    mutate(time = as.factor(time)) %>% 
    group_by(cell_line,treatment,time, marker) %>%
    mutate(best_in_condition = RMSE == min(RMSE)) %>% 
    filter(best_in_condition) 

mean(team_independent_predictions$RMSE)


### RMSE violin plots --------------------------------------
fill_color <- RColorBrewer::brewer.pal(8,"Dark2")[[1]]
# Cell_line HCC2218 shows significantly larger prediction errors than other cell-lines
plot_rmse_cell_line <- team_independent_predictions %>% 
    ggplot(aes(cell_line,RMSE,fill=best_in_condition)) +
    #geom_hline(yintercept = median(rmse_data$icx_bxai), size=0.3) + 
    geom_violin(draw_quantiles = c(0.5),fill= fill_color) +
    # scale_fill_manual(values = c("TRUE"="grey","FALSE"=RColorBrewer::brewer.pal(8,"Dark2")[[1]])) + 
    theme_bw() + 
    guides(fill = "none") + 
    #ggtitle("RMSE by cell-lines") +
    labs( ylab= NULL) + xlab("Cell Lines")


plot_rmse_time <- team_independent_predictions %>%
    ggplot(aes(time,RMSE,fill=best_in_condition)) +
    #geom_hline(yintercept = median(rmse_data$icx_bxai), size=0.3) + 
    geom_violin(draw_quantiles = c(0.5),fill= fill_color) +
    # scale_fill_manual(values = c("TRUE"="grey","FALSE"=RColorBrewer::brewer.pal(8,"Dark2")[[1]])) + 
    theme_bw() + 
    guides(fill = "none") + 
    #ggtitle("RMSE by time-points") +
    labs( ylab= NULL) +
    xlab("Time points") 

plot_rmse_marker <- team_independent_predictions %>% 
    ggplot(aes(marker,RMSE,fill=best_in_condition)) +
    #geom_hline(yintercept = median(rmse_data$icx_bxai), size=0.3) + 
    geom_violin(draw_quantiles = c(0.5),fill= fill_color) +
    # scale_fill_manual(values = c("TRUE"="grey","FALSE"=)) + 
    theme_bw() + 
    guides(fill = "none") + 
    #ggtitle("RMSE by markers ") +
    labs( ylab= NULL) +
    xlab("Marker") 

plot_rmse_treatment <- team_independent_predictions %>% 
    ggplot(aes(treatment,RMSE)) +
    #geom_hline(yintercept = median(rmse_data$icx_bxai), size=0.3) + 
    geom_violin(draw_quantiles = c(0.5),fill= fill_color) +
    #scale_fill_manual(values = c("TRUE"="grey","TRUE"=RColorBrewer::brewer.pal(8,"Dark2")[[1]])) + 
    theme_bw() + 
    guides(fill = "none") + 
    #ggtitle("RMSE by treatment ") + 
    xlab("Treatment") +
    labs( ylab= NULL) #+ylab("RMSE per condition")

cowplot::plot_grid(plot_rmse_cell_line,plot_rmse_marker,plot_rmse_time, plot_rmse_treatment)

ggsave("./publication/figures/supp_figures/sc1_prediction_error_per_condition.pdf",width = 8,height = 6)

# HCC2218 and p.S6 was the most difficult to predict in SC1



###  Anova approach ----
aov_res <- team_independent_predictions %>% ungroup() %>% 
    mutate_at(c("cell_line","treatment","time","marker"),as.factor) %>%
    aov(data = ., RMSE ~ cell_line+treatment+time+marker ) %>% summary(.)

gt_data <- as.data.frame(aov_res[[1]]) %>% mutate_each(~format(.,digits = 2))
gt_table <- gridExtra::grid.table(gt_data)

ggsave(filename = "./publication/figures/supp_figures/sc1_prediction_error_annova_table.pdf",
       width = 4,height = 3)

team_independent_predictions %>% ungroup() %>% 
    mutate_at(c("cell_line","treatment","time","marker"),as.factor) %>%
    aov(data = ., RMSE ~ cell_line*marker + cell_line + marker) %>% coefficients(.)


team_independent_predictions %>% mutate()

rmse_data %>% select(cell_line, treatment,  time, marker, RMSE ) %>%
    mutate(y_id = paste(treatment,time),
           x_id = paste(marker,cell_line)) %>%
    arrange(marker,treatment,time) %>% ungroup() %>% select(-cell_line, -marker,-treatment,-time) %>%
    spread(x_id,RMSE) %>%
    column_to_rownames("y_id") %>% 
    pheatmap::pheatmap(mat = .,cluster_cols = F,cluster_rows = F)


# this looks really good: 
# easy to understand which combinations are the most difficult: 

team_independent_predictions %>% select(cell_line, treatment,  time, marker, RMSE ) %>%
    mutate(y_id = paste(treatment,time),
           x_id = paste(marker,cell_line)) %>%
    arrange(marker,treatment,time) %>% ungroup() %>% select(-cell_line, -marker,-treatment,-time) %>%
    spread(x_id,RMSE) %>%
    column_to_rownames("y_id") %>% 
    pheatmap::pheatmap(mat = .,cluster_cols = F,cluster_rows = F, )



pt <- team_independent_predictions %>% select(cell_line, treatment,  time, marker, RMSE ) %>%
    mutate(y_id = paste(treatment,time),
           x_id = paste(marker,cell_line)) %>%
    arrange(marker,treatment,time) %>% ungroup() %>% select(-cell_line, -marker,-treatment,-time) %>%
    spread(x_id,RMSE) %>%
    column_to_rownames("y_id") %>% 
    pheatmap::pheatmap(mat = .,cluster_cols = F,cluster_rows = F,
                       breaks = seq(0.4,2,length.out = 50),
                       color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(49),
                       main = "RMSE of best in condition",
                       gaps_col = c(6,12,18,24,30), 
                       gaps_row = c(12,13,22,31,40),
    )
ggsave("./submission_analysis/figures/sc1_condition_error/best_in_cond_heatmap.pdf",height = 8,width = 6,plot = pt$gtable)


#########################
## Plotting dyanamic (time-dependent responses) for different variables/ conditions

top_teams = c()


# 1. p-ERK in EGF treatment: palteau and peaking effect: -----------------------
median_data %>% filter(treatment=="EGF", marker=="p.ERK") %>%
    gather(team,prediction,6:27) %>%
    ggplot(aes(time,standard)) + geom_line() + 
    facet_wrap(~cell_line) +
    ylab("p.ERK") + 
    ggtitle("p.ERK signaling in EGF stimulation") + 
    theme_bw()


median_data %>% filter(treatment=="EGF", marker=="p.ERK") %>%
    gather(team,prediction,6:27) %>%
    filter(team %in% SC_leaderboard$submitterId[1:5]) %>%
    ggplot() + 
    geom_line(aes(time,standard)) +
    geom_line(aes(time,prediction,color=cell_line)) +
    ylab("p.ERK signaling") +
    ggtitle("p.ERK signaling in EGF stimulation", subtitle = "top 5 teams") + 
    facet_grid(team~cell_line) +
    theme_bw()


# 2. MEK independent p-ERK signaling:  -------------------------------
median_data %>% filter(treatment %in% c("EGF","iMEK"), marker=="p.ERK") %>%
    gather(team,prediction,6:27) %>%
    ggplot(aes(time,standard,color=treatment)) +
    geom_line() + 
    ylab("p.ERK signaling") + 
    ggtitle("p.ERK signaling in iMEK inhibition") +
    facet_wrap(~cell_line) + theme_bw()


median_data %>% filter(treatment %in% c("EGF","iMEK"), marker=="p.ERK") %>%
    gather(team,prediction,6:27) %>%
    filter(team %in% c(SC_leaderboard$submitterId[1:4],"X.msinkala","X.pqiu")) %>%
    ggplot() + 
    geom_line(aes(time,standard,color=treatment)) +
    geom_line(aes(time,prediction,color=treatment),linetype=2) +
    facet_grid(team~cell_line) + theme_bw() + 
    ylab("p.ERK activity") +
    ggtitle("MEK independet p.ERK", subtitle = "top 4 teams + 2 that predicts well HCC2218")

median_data %>% filter(treatment %in% c("EGF","iMEK"), marker=="p.ERK", cell_line == "HCC2218") %>%
    gather(team,prediction,6:27) %>%
    ggplot() + 
    geom_line(aes(time,standard,color=treatment)) +
    geom_line(aes(time,prediction,color=treatment),linetype=2) +
    facet_wrap(~team) + theme_bw() +
    ggtitle("MEK independet p.ERK in HCC2218", subtitle = "all teams")

# teams that predict almost no change to iMEK inhibition
median_data %>% filter(treatment %in% c("EGF","iMEK"), marker=="p.ERK") %>%
    select(1:5,hulab.SCS, X.pqiu, X.msinkala,X.Huiyuan) %>%
    gather(team,prediction,hulab.SCS, X.pqiu, X.msinkala,X.Huiyuan) %>%
    ggplot() + 
    geom_line(aes(time,standard,color=treatment)) +
    geom_line(aes(time,prediction,color=treatment),linetype=2) +
    facet_grid(team~cell_line) + theme_bw() +
    ggtitle("MEK independet p.ERK in HCC2218", subtitle = "all teams")



# 2.2. MEK independent p-ERK signaling: pS6 -------------------------------

median_data %>% filter(treatment %in% c("EGF","iMEK"), marker=="p.S6") %>%
    gather(team,prediction,6:27) %>%
    ggplot(aes(time,standard,color=treatment)) + geom_line() +  facet_wrap(~cell_line) + theme_bw()


median_data %>% filter(treatment %in% c("EGF","iMEK"), marker=="p.S6") %>%
    gather(team,prediction,6:27) %>%
    filter(team %in% SC_leaderboard$submitterId[1:5]) %>%
    ggplot() + 
    geom_line(aes(time,standard,color=treatment)) +
    geom_line(aes(time,prediction,color=treatment),linetype=2) +
    facet_grid(team~cell_line) + theme_bw() + 
    ggtitle("MEK independet p.S6", subtitle = "top 5 teams")

median_data %>% filter(treatment %in% c("EGF","iMEK"), marker=="p.ERK", cell_line == "HCC2218") %>%
    gather(team,prediction,6:27) %>%
    ggplot() + 
    geom_line(aes(time,standard,color=treatment)) +
    geom_line(aes(time,prediction,color=treatment),linetype=2) +
    facet_wrap(~team) + theme_bw() +
    ggtitle("MEK independet p.S6 in HCC2218", subtitle = "all teams")


# 3. PI3K independent AKT signaling -------------------------------------------


median_data %>% filter(treatment %in% c("EGF","iPI3K"), marker=="p.Akt.Ser473.") %>%
    gather(team,prediction,6:27) %>%
    ggplot(aes(time,standard,color=treatment)) + geom_line() +
    facet_wrap(~cell_line) +
    ylab("p.Akt.Ser473.") + 
    theme_bw()


median_data %>% filter(treatment %in% c("EGF","iPI3K"), marker=="p.Akt.Ser473.") %>%
    gather(team,prediction,6:27) %>%
    filter(team %in% SC_leaderboard$submitterId[1:5]) %>%
    ggplot() + 
    geom_line(aes(time,standard,color=treatment)) +
    geom_line(aes(time,prediction,color=treatment),linetype=2) +
    facet_grid(team~cell_line) + theme_bw() + 
    ylab("p.Akt.Ser473.") + 
    ggtitle("PI3K independet p.Akt.Ser473", subtitle = "top 5 teams")

median_data %>% filter(treatment %in% c("EGF","iPI3K"), marker=="p.Akt.Ser473.", cell_line == "AU565") %>%
    gather(team,prediction,6:27) %>%
    ggplot() + 
    geom_line(aes(time,standard,color=treatment)) +
    geom_line(aes(time,prediction,color=treatment),linetype=2) +
    facet_wrap(~team) + theme_bw() +
    ylab("p.Akt.Ser473.") + 
    ggtitle("PI3K dependet p.Akt.Ser473. in AU565", subtitle = "all teams")


# 4. EGFR independent signaling -------------------------------------------

median_data %>% filter(treatment %in% c("EGF","iEGFR")) %>%
    gather(team,prediction,6:27) %>%
    ggplot(aes(time,standard,color=marker)) + geom_line(aes(linetype=treatment)) +
    facet_grid(marker~cell_line) +
    ylab("value") + 
    theme_bw()


median_data %>% filter(treatment %in% c("EGF","iEGFR"), cell_line=="AU565") %>%
    gather(team,prediction,6:27) %>%
    filter(team %in% SC_leaderboard$submitterId[1:10]) %>%
    ggplot() + 
    geom_line(aes(time,standard,color=treatment)) +
    geom_line(aes(time,prediction,color=treatment),linetype=2) +
    facet_grid(marker~team,scales = "free_y") + theme_bw() + 
    ylab("markers") + 
    ggtitle("iEGFR dependet signaling in AU565", subtitle = "top 10 teams")

median_data %>% filter(treatment %in% c("EGF","iPI3K"), marker=="p.Akt.Ser473.", cell_line == "AU565") %>%
    gather(team,prediction,6:27) %>%
    ggplot() + 
    geom_line(aes(time,standard,color=treatment)) +
    geom_line(aes(time,prediction,color=treatment),linetype=2) +
    facet_wrap(~team) + theme_bw() +
    ylab("p.Akt.Ser473.") + 
    ggtitle("PI3K dependet p.Akt.Ser473. in AU565", subtitle = "all teams")


# 4. S6 modelling  -------------------------------------------

median_data %>% filter(marker %in% c('p.S6')) %>%
    gather(team,prediction,6:27) %>%
    ggplot(aes(time,standard,color=treatment)) + geom_line() +
    facet_grid(treatment~cell_line) +
    ylab("value") + 
    theme_bw() + ggtitle("S6 signaling in standards")


median_data %>% filter(treatment %in% c("EGF","iPKC"), marker =="p.S6") %>%
    gather(team,prediction,6:27) %>%
    filter(team %in% SC_leaderboard$submitterId[1:5]) %>%
    ggplot() + 
    geom_line(aes(time,standard,color=treatment)) +
    geom_line(aes(time,prediction,color=treatment),linetype=2) +
    facet_grid(team~cell_line) + theme_bw() + 
    ylab("markers") + 
    ggtitle("iPKC dependent S6 signaling", subtitle = "top 10 teams")

median_data %>% filter(treatment %in% c("EGF","iPI3K"), marker=="p.Akt.Ser473.", cell_line == "AU565") %>%
    gather(team,prediction,6:27) %>%
    ggplot() + 
    geom_line(aes(time,standard,color=treatment)) +
    geom_line(aes(time,prediction,color=treatment),linetype=2) +
    facet_wrap(~team) + theme_bw() +
    ylab("p.Akt.Ser473.") + 
    ggtitle("PI3K dependet p.Akt.Ser473. in AU565", subtitle = "all teams")

# 5. HER2 signaling  -----------------------------------------------------------
median_data %>% filter(cell_line %in% c("AU565","HCC2218"), marker=="p.ERK") %>%
    gather(team,prediction,6:27) %>%
    #mutate(time= as_factor(time)) %>% 
    ggplot(aes(time,standard,color=treatment)) +
    geom_line() +
    geom_boxplot(aes(time,prediction,group=paste(time,treatment)),position = "identity",outlier.alpha=0) +
    facet_grid(cell_line~treatment) + theme_bw()




### MEK independent ERK signaling 
# we compare 2 cell-lines ("AU565","HCC2218") 
plt_data <- median_data %>% 
    filter(cell_line %in% c("AU565","HCC2218")) %>%
    filter(marker=="p.ERK") %>%
    filter(treatment %in% c("iMEK","EGF")) %>%
    gather(team,prediction,6:27) 
#mutate(time= as_factor(time)) %>% 
plt_data %>% ggplot() +
    geom_boxplot(data = filter(plt_data,treatment %in% c("iMEK")), aes(time,prediction,group=time,fill="iMEK prediction"),
                 color=RColorBrewer::brewer.pal(8,"Dark2")[[2]],outlier.size = 0.5) +
    geom_line(data=filter(plt_data,treatment %in% c("iMEK"),team=="icx_bxai"), aes(time,prediction,col="iMEK prediction")) +
    geom_point(data=filter(plt_data,treatment %in% c("iMEK")), aes(time,standard,col=treatment)) +
    geom_line(aes(time,standard,col=treatment)) +
    facet_grid(~cell_line) + theme_bw() + 
    scale_color_manual(name = "", values = c("EGF" = "grey60", 
                                             "iMEK"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],
                                             "iMEK prediction"=RColorBrewer::brewer.pal(8,"Dark2")[[2]]), ) +
    scale_fill_manual(name = "",values = c("iMEK prediction"=RColorBrewer::brewer.pal(8,"Dark2")[[2]])) +
    ylab("p-ERK intensity") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggsave("./publication/figures/figure2/sc1_mek_indep_ERK.pdf",width = 5,height = 2.5)

### iEGFR independent ERK signaling -- skipped
# we compare 2 cell-lines ("AU565","HCC2218") 
plt_data2 <- median_data %>% 
    filter(cell_line %in% c("AU565","HCC2218")) %>%
    filter(marker=="p.ERK") %>%
    filter(treatment %in% c("iEGFR","EGF")) %>%
    gather(team,prediction,6:27) 
#mutate(time= as_factor(time)) %>% 
plt_data2 %>% ggplot() +
    geom_boxplot(data = filter(plt_data2,treatment %in% c("iEGFR")), aes(time,prediction,group=time),
                 fill=RColorBrewer::brewer.pal(8,"Dark2")[[2]],color=RColorBrewer::brewer.pal(8,"Dark2")[[2]]) +
    geom_point(aes(time,standard,col=treatment)) +
    geom_line(aes(time,standard,col=treatment)) +
    geom_line(data=filter(plt_data2,treatment %in% c("iEGFR"),team=="icx_bxai"), aes(time,prediction,col="iEGFR prediction")) +
    facet_grid(~cell_line) + theme_bw() + 
    scale_color_manual(name = "", values = c("EGF" = "grey60", 
                                             "iEGFR"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],
                                             "iEGFR prediction"=RColorBrewer::brewer.pal(8,"Dark2")[[2]])) +
    ylab("p-ERK intensity") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




### iEGFR independent AKT signaling  ----------
# we compare 2 cell-lines ("AU565","LY2") 
plt_data3 <- median_data %>% 
    filter(cell_line %in% c("AU565","LY2")) %>%
    filter(marker=="p.Akt.Ser473.") %>%
    filter(treatment %in% c("iEGFR","EGF")) %>%
    gather(team,prediction,6:27) 
#mutate(time= as_factor(time)) %>% 
plt_data3 %>% ggplot() +
    geom_boxplot(data = filter(plt_data3,treatment %in% c("iEGFR")), aes(time,prediction,group=time,fill="iEGFR prediction"),
                 color=RColorBrewer::brewer.pal(8,"Dark2")[[2]],outlier.size = 0.5) +
    geom_line(data=filter(plt_data3,treatment %in% c("iEGFR"),team=="icx_bxai"), aes(time,prediction,col="iEGFR prediction")) +
    geom_point(data=filter(plt_data3,treatment %in% c("iEGFR")),aes(time,standard,col=treatment)) +
    geom_line(aes(time,standard,col=treatment)) +
    facet_grid(~cell_line) + theme_bw() + 
    scale_color_manual(values = c("EGF" = "grey60",
                                  "iEGFR"=RColorBrewer::brewer.pal(8,"Dark2")[[1]],
                                  "iEGFR prediction"=RColorBrewer::brewer.pal(8,"Dark2")[[2]])) +
    scale_fill_manual(name = "",values = c("iEGFR prediction"=RColorBrewer::brewer.pal(8,"Dark2")[[2]])) +
    ylab(bquote('p-'*~Akt^Ser473*' intensity' ) ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggsave("./publication/figures/figure2/sc1_egfr_indep_Akt.pdf",width = 5,height = 2.5)


### iMEK independent AKT signaling 
# we compare 2 cell-lines ("AU565","LY2") 
plt_data3 <- median_data %>% 
    filter(cell_line %in% c("AU565","LY2")) %>%
    filter(marker=="p.Akt.Ser473.") %>%
    filter(treatment %in% c("iEGFR","EGF")) %>%
    gather(team,prediction,6:27) 
#mutate(time= as_factor(time)) %>% 
plt_data3 %>% ggplot() +
    geom_boxplot(data = filter(plt_data3,treatment %in% c("iEGFR")), aes(time,prediction,group=time),
                 fill=RColorBrewer::brewer.pal(8,"Dark2")[[2]],color=RColorBrewer::brewer.pal(8,"Dark2")[[2]]) +
    geom_line(data=filter(plt_data3,treatment %in% c("iEGFR"),team=="icx_bxai"), aes(time,prediction),color=RColorBrewer::brewer.pal(8,"Dark2")[[2]]) +
    geom_point(aes(time,standard,col=treatment)) +
    geom_line(aes(time,standard,col=treatment)) +
    facet_grid(~cell_line) + theme_bw() + 
    scale_color_manual(values = c("EGF" = "grey60","iEGFR"=RColorBrewer::brewer.pal(8,"Dark2")[[1]])) +
    ylab("p-Akt.Ser473. intensity")

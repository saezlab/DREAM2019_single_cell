# schow a good fitting example for SC2



# eg. pMEK and pERK



ics_bxai  <-  "./submission_data/final_round/SC2/9695439.csv"
gs_file <- "./challenge_data/validation_data/sc2gold.csv"

combined_statistics <- read_rds("./submission_analysis/intermediate_data/sc2_stats_conditions.rds")



combined_statistics %>% select(cell_line, treatment,time, stat_variable, standard, icx_bxai) %>%
    mutate(error = (standard-icx_bxai)^2) %>%
    filter(stat_variable == "cov_p.ERK_p.MEK") %>% arrange(desc(standard))
# 1 MDAMB468  iEGFR         7 cov_p.ERK_p.MEK 

prediction = read_csv(ics_bxai)
gs_data <- read_csv(gs_file)




selected_prediction <- prediction %>% filter(cell_line =="MDAMB468", treatment == "iEGFR", time ==7  ) 
selected_gs <- gs_data %>% filter(cell_line =="MDAMB468", treatment == "iEGFR", time ==7  ) 


RColorBrewer::display.brewer.all()
library(RColorBrewer)

d1 <- ggplot() + 
    geom_density(data=selected_prediction,aes(x=p.ERK),fill=brewer.pal(12,"Paired")[4],alpha=0.5) + 
    geom_density(data=selected_gs,aes(x=p.ERK),fill=brewer.pal(12,"Paired")[2],alpha=0.5) +
    theme_bw() + coord_flip() +scale_y_reverse() +
    theme( )

d2 <- ggplot() + 
    geom_density(data=selected_prediction,aes(x=p.MEK),fill=brewer.pal(12,"Paired")[4],alpha=0.5) + 
    geom_density(data=selected_gs,aes(x=p.MEK), fill=brewer.pal(12,"Paired")[2],alpha=0.5) +
    theme_bw()+ scale_y_reverse() + theme()

d3 <- ggplot() + 
    geom_point(data=selected_prediction,aes(x=p.MEK,y=p.ERK,col="prediction"),size=.5) + 
    geom_point(data=selected_gs,aes(x=p.MEK,y=p.ERK,col="data"),size=.5) +
    theme_bw() +xlab("") + ylab("") + theme(axis.text = element_blank()) + 
    scale_color_manual(values = c(prediction=brewer.pal(12,"Paired")[4],
                                  data=brewer.pal(12,"Paired")[2])) +
    theme(legend.position = c(0.8, 0.2),legend.title = element_blank())


cowplot::plot_grid(d1,d3,NULL,d2,ncol = 2,rel_widths = c(.3,.8),rel_heights = c(.8,.3),align = "hv",greedy=TRUE)
ggsave("./publication/figures/figure3/sc2_MEK_ERK_correlation_example.pdf",width = 4,height = 4)


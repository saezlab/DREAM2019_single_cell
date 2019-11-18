library(RColorBrewer)
library(tidyverse)
# TODO: download the leaderboard data to an excel file! 

#leader_board_raw = "./submission_data/round2/leaderboard_round2.xlsx"
leader_board_raw = "./submission_data/final_round/leaderboard_all_rounds.xlsx"

source("./submission_analysis/utilities.R")


LB_sc1 <- readxl::read_excel(leader_board_raw,sheet = "sc1") %>% fix_leaderboard_raw()
LB_sc2 <- readxl::read_excel(leader_board_raw,sheet = "sc2") %>% fix_leaderboard_raw()
LB_sc3 <- readxl::read_excel(leader_board_raw,sheet = "sc3") %>% fix_leaderboard_raw()
LB_sc4 <- readxl::read_excel(leader_board_raw,sheet = "sc4") %>% fix_leaderboard_raw()




#RColorBrewer::display.brewer.all()
colors  = brewer.pal(n = 8,"Dark2")

### Subchallange I -------------------------------------------------------------

LB_sc1 %>% ggplot() + 
	geom_rect(xmin=0.7, xmax=1.25, ymin=0, ymax=40, fill ="#F8F8CE",color="black")+
	geom_histogram(aes(score),bins = 50, fill = colors[[1]]) + theme_bw() + 
	ggtitle("Sub-challenge I: predict missing marker") +
	ylab("# of submissions") + 
	geom_vline(xintercept = 0.903103,color = colors[[2]],size=2) + 
	geom_vline(xintercept = 3.992752,color = colors[[3]],size=2) + coord_cartesian(xlim = c(0,6))

ggsave("./submission_analysis/figures/round2_sc1_score_dist.pdf",width = 6,height = 2)

LB_sc1 %>% filter(score<1.2) %>% ggplot() + geom_histogram(aes(score),bins = 40, fill = colors[[1]]) + theme_bw() + 
	ggtitle("Sub-challenge I: predict missing marker") +
	ylab("# of submissions") + 
	geom_vline(xintercept = 0.903103,color = colors[[2]],size = 2) + 
	xlim(0.7,1.25) + theme(panel.background = element_rect(fill = "#F8F8CE" ))  
	# geom_segment(x = 0.903103, y=20 ,xend =0.903103, yend = 0, 
	# 			 arrow.fill = colors[[2]],
	# 			 color = colors[[2]],
	# 			 arrow =  arrow(length = unit(1,"cm"),type = "closed",angle = 20),size=1)
	# #coord_cartesian(xlim = c(0,6))
	
ggsave("./submission_analysis/figures/round2_sc1_score_dist_top.pdf",width = 6,height = 2)

# Statistics: 
print("# Submissions:")
LB_sc1 %>% filter(is.finite(score)) %>% nrow()

print("# individuals:")
unique(LB_sc1$Submitter) %>% length()


### Subchallange II -------------------------------------------------------------

LB_sc2 %>% filter(score < 200) %>% ggplot() + geom_histogram(aes(score),bins = 20, fill = colors[[1]]) + theme_bw() + 
	ggtitle("Sub-challenge II: predict missing treatment") +
	ylab("# of submissions")  +
	xlab("score (mean + covariance)") + 
	geom_vline(xintercept = 44.6702,color = colors[[2]],size=2) + 
	geom_vline(xintercept = 179.613971,color = colors[[3]],size=2)

ggsave("./submission_analysis/figures/round2_sc2_score_dist.pdf",width = 6,height = 2)

# Statistics: 
print("# Submissions:")
LB_sc2 %>% filter(is.finite(score)) %>% nrow()

print("# individuals:")
unique(LB_sc2$Submitter) %>% length()


### Subchallange III -----------------------------------------------------------

LB_sc3 %>% filter(score < 300) %>% ggplot() + geom_histogram(aes(score),bins = 20, fill = colors[[1]]) + theme_bw() + 
	ggtitle("Sub-challenge III: predict new treatment") +
	ylab("# of submissions")  +
	xlab("score (mean + covariance)") + 
	geom_vline(xintercept = 336.76,color = colors[[3]],size=2)

ggsave("./submission_analysis/figures/round2_sc3_score_dist.pdf",width = 6,height = 2)

# Statistics: 
print("# Submissions:")
LB_sc3 %>% filter(is.finite(score)) %>% nrow()

print("# individuals:")
unique(LB_sc3$Submitter) %>% length()


### Subchallange IV -----------------------------------------------------------

LB_sc4 %>% ggplot() + 
	geom_rect(xmin=0.2, xmax=0.55, ymin=-1, ymax=27, fill ="#F8F8CE",color="black")+
	geom_histogram(aes(score),bins = 20, fill = colors[[1]]) + theme_bw() + 
	ggtitle("Sub-challenge IV: predict median time-course") +
	ylab("# of submissions")  +
	xlab("RMSE") + 
	geom_vline(xintercept = 0.340827,color = colors[[2]],size=2)+
	geom_vline(xintercept = 2.95,color = colors[[3]],size=2)

ggsave("./submission_analysis/figures/round2_sc4_score_dist.pdf",width = 6,height = 2)


LB_sc4 %>% filter(score<0.45) %>% ggplot() +
	geom_histogram(aes(score),bins = 40, fill = colors[[1]]) + theme_bw() + 
	ggtitle("Sub-challenge IV: predict median time-course") +
	ylab("# of submissions") + 
	xlab("RMSE") + 
	geom_vline(xintercept = 0.340827,color = colors[[2]],size=2)+ theme(panel.background = element_rect(fill = "#F8F8CE" ))  

ggsave("./submission_analysis/figures/round2_sc4_score_dist_top.pdf",width = 6,height = 2)

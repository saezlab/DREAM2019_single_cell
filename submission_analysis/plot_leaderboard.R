library(RColorBrewer)
library(tidyverse)

LB_sc1 <- readxl::read_excel("./submission_data/round1/leaderboard_round1.xlsx",sheet = "sc1")
LB_sc2 <- readxl::read_excel("./submission_data/round1/leaderboard_round1.xlsx",sheet = "sc2")
LB_sc3 <- readxl::read_excel("./submission_data/round1/leaderboard_round1.xlsx",sheet = "sc3")
LB_sc4 <- readxl::read_excel("./submission_data/round1/leaderboard_round1.xlsx",sheet = "sc4")

#RColorBrewer::display.brewer.all()
colors  = brewer.pal(n = 8,"Dark2")

### Subchallange I -------------------------------------------------------------

LB_sc1 %>% ggplot() + geom_histogram(aes(RMSE),bins = 20, fill = colors[[1]]) + theme_bw() + 
	ggtitle("Sub-challenge I: predict missing marker") +
	ylab("# of submissions") + 
	geom_vline(xintercept = 0.903103,color = colors[[2]]) + 
	geom_vline(xintercept = 3.992752,color = colors[[3]])

ggsave("./submission_analysis/figures/round1_sc1_score_dist.pdf",width = 6,height = 2)

LB_sc1 %>% filter(RMSE<1.2) %>% ggplot() + geom_histogram(aes(RMSE),bins = 40, fill = colors[[1]]) + theme_bw() + 
	ggtitle("Sub-challenge I: predict missing marker") +
	ylab("# of submissions") + 
	geom_vline(xintercept = 0.903103,color = colors[[2]])
	
ggsave("./submission_analysis/figures/round1_sc1_score_dist_top.pdf",width = 6,height = 2)


### Subchallange II -------------------------------------------------------------

LB_sc2 %>% filter(RMSE < 200) %>% ggplot() + geom_histogram(aes(RMSE),bins = 20, fill = colors[[1]]) + theme_bw() + 
	ggtitle("Sub-challenge II: predict missing treatment") +
	ylab("# of submissions")  +
	xlab("score (mean + covariance)") + 
	geom_vline(xintercept = 44.6702,color = colors[[2]],size=2) + 
	geom_vline(xintercept = 179.613971,color = colors[[3]],size=2)

ggsave("./submission_analysis/figures/round1_sc2_score_dist.pdf",width = 6,height = 2)


### Subchallange III -----------------------------------------------------------

LB_sc3 %>% filter(RMSE < 300) %>% ggplot() + geom_histogram(aes(RMSE),bins = 20, fill = colors[[1]]) + theme_bw() + 
	ggtitle("Sub-challenge III: predict new treatment") +
	ylab("# of submissions")  +
	xlab("score (mean + covariance)") + 
	geom_vline(xintercept = 336.76,color = colors[[3]],size=2)

ggsave("./submission_analysis/figures/round1_sc3_score_dist.pdf",width = 6,height = 2)

### Subchallange IV -----------------------------------------------------------

LB_sc4 %>% ggplot() + geom_histogram(aes(RMSE),bins = 20, fill = colors[[1]]) + theme_bw() + 
	ggtitle("Sub-challenge IV: predict median time-course") +
	ylab("# of submissions")  +
	xlab("RMSE") + 
	geom_vline(xintercept = 0.340827,color = colors[[2]],size=2)+
	geom_vline(xintercept = 2.95,color = colors[[3]],size=2)

ggsave("./submission_analysis/figures/round1_sc4_score_dist.pdf",width = 6,height = 2)


LB_sc4 %>% filter(RMSE<0.45) %>% ggplot() + geom_histogram(aes(RMSE),bins = 40, fill = colors[[1]]) + theme_bw() + 
	ggtitle("Sub-challenge IV: predict median time-course") +
	ylab("# of submissions") + 
	xlab("RMSE") + 
	geom_vline(xintercept = 0.340827,color = colors[[2]],size=2)

ggsave("./submission_analysis/figures/round1_sc4_score_dist_top.pdf",width = 6,height = 2)

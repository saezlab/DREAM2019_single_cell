# Round 1 fuinished
# check the first 3 teams for each challenge and download their submissions for 
# local analysis.


library(synapser)
library(tidyverse)


synLogin()

submission_folder = "./submission_data/round2/"
leader_board_raw = "./submission_data/round2/leaderboard_round2.xlsx"




# dir.create(path = file.path(submission_folder,"SC1"))
# dir.create(path = file.path(submission_folder,"SC2"))
# dir.create(path = file.path(submission_folder,"SC3"))
# dir.create(path = file.path(submission_folder,"SC4"))
		   		   

source("./submission_analysis/utilities.R")

LB_sc1 <- readxl::read_excel(leader_board_raw,sheet = "sc1") %>% fix_leaderboard_raw()
LB_sc2 <- readxl::read_excel(leader_board_raw,sheet = "sc2") %>% fix_leaderboard_raw()
LB_sc3 <- readxl::read_excel(leader_board_raw,sheet = "sc3") %>% fix_leaderboard_raw()
LB_sc4 <- readxl::read_excel(leader_board_raw,sheet = "sc4") %>% fix_leaderboard_raw()



# Download the top 3 submissions at Round 1 for each challenge


leaderboard = tibble(subchallenge = paste0(rep("SC",3),rep(1:4,each=3)), 
						 submission =  c(
						 	c(9694472, 9694390, 9693731), # SC1
							c(9694476, 9694448, 9694393), # SC2
							c(9694391, 9694466, 9694375), # SC3
							c(9694457, 9694465, 9694445)  # SC4
							))

# downloads the submission and moves to the designed folder:
download_submission <- function(submissionID,subchallenge_folder){
	sub = synGetSubmission(submissionID)
	file.copy(sub$filePath,to = file.path(subchallenge_folder,paste0(submissionID,".csv")))
	return(sub$filePath)
}


leaderboard %>% rowwise() %>%  do(download_submission( submissionID = .$submission,
										subchallenge_folder = file.path(submission_folder,.$subchallenge )))




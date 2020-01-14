
# check the first 3 teams for each challenge and download their submissions for 
# local analysis.


library(synapser)
library(tidyverse)


synLogin()

submission_folder = "./submission_data/final_round/"


# dir.create(path = file.path(submission_folder,"SC1"))
# dir.create(path = file.path(submission_folder,"SC2"))
# dir.create(path = file.path(submission_folder,"SC3"))
# dir.create(path = file.path(submission_folder,"SC4"))
		   		   

source("./submission_analysis/utilities.R")


# Download the top 3 submissions at Round 1 for each challenge


leaderboard = tibble(subchallenge = paste0(rep("SC",3),rep(1:4,each=3)), 
						 submission =  c(
						 	c(9695437, 9695091, 9695159), # SC1
							c(9695439, 9694578, 9695468), # SC2
							c(9695436, 9695469, 9695270), # SC3
							c(9695440, 9695158, 9695472)  # SC4
							))

# downloads the submission and moves to the designed folder:
download_submission <- function(submissionID,subchallenge_folder){
	sub = synGetSubmission(submissionID)
	file.copy(sub$filePath,to = file.path(subchallenge_folder,paste0(submissionID,".csv")))
	return(sub$filePath)
}


leaderboard %>% rowwise() %>%  do(download_submission( submissionID = .$submission,
										subchallenge_folder = file.path(submission_folder,.$subchallenge )))



submissionID = "9695178"

sub = synGetSubmission(submissionID)file.copy()
file.copy(sub$filePath,to = file.path(subchallenge_folder,paste0(submissionID,".csv")))
filesstrings::file.move(sub$filePath, file.path(subchallenge_folder,paste0(submissionID,".csv")))

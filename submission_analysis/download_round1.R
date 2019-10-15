# Round 1 fuinished
# check the first 3 teams for each challenge and download their submissions for 
# local analysis.


library(synapser)
library(tidyverse)


synLogin()

submission_folder = "./submission_data/round1/"


# dir.create(path = file.path(submission_folder,"SC1"))
# dir.create(path = file.path(submission_folder,"SC2"))
# dir.create(path = file.path(submission_folder,"SC3"))
# dir.create(path = file.path(submission_folder,"SC4"))
		   		   

# Download the top 3 submissions at Round 1 for each challenge


leaderboard = tibble(subchallenge = paste0(rep("SC",3),rep(1:4,each=3)), 
						 submission =  c(
						 	c(9693326, 9693414, 9693308), # SC1
							c(9693420, 9693418, 9693391), # SC2
							c(9693423, 9693317, 9693386), # SC3
							c(9693383, 9693336, 9693377)  # SC4
							))

# downloads the submission and moves to the designed folder:
download_submission <- function(submissionID,subchallenge_folder){
	sub = synGetSubmission(submissionID)
	file.copy(sub$filePath,to = file.path(subchallenge_folder,paste0(submissionID,".csv")))
	return(sub$filePath)
}


leaderboard %>% rowwise() %>%  do(download_submission( submissionID = .$submission,
										subchallenge_folder = file.path(submission_folder,.$subchallenge )))




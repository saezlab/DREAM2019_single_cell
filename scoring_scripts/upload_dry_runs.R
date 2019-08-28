# submit dry run results
library(tidyverse)
library(synapser)

####################
## subchallenge 1 ##
####################

# prepare the raw predictions and then submit to synapse
predictions_raw  =read_csv("./dry_run/aim1_predictions.csv")
template = read_csv("./challenge_data/prediction_templates/subchallenge_1_template_data.csv")


predictions_raw %>% full_join(template %>% 
							  	select(glob_cellID,cell_line, treatment, time, cellID, fileID),
							  by = c("cell_line", "treatment", "time", "cellID", "fileID")) %>% write_csv("dry_run/subchallenge_prediction_v1.csv")


library(synapser)

synLogin(email="attilagabor")

staging_challenge_project = "syn20366916"

base_folder_ent = Folder("test_submission", parentId=staging_challenge_project)
stored_base_folder_ent = synStore(base_folder_ent)


template_data_ent = File("dry_run/subchallenge_prediction_v1.csv", parentId=stored_base_folder_ent$properties$id)
synStore(template_data_ent)


####################
## subchallenge 2 ##
####################

# prepare the raw predictions and then submit to synapse
predictions_raw =read_csv("./dry_run/aim1_2_1_predictions.csv")
template = read_csv("./challenge_data/prediction_templates/subchallenge_2_template_data.csv")




####################
## subchallange 4 ##
####################

predictions_raw <- read_csv("./dry_run/aim2_model_predictions_test.csv")

synLogin(email="attilagabor", apiKey="/YwdwbcFPKfTSazwKChjXSXd/ZZ8BQ0DXOkh/JuGWvmGtnZCfUfRuZ5Ixid5RjGdvH0J8QlUImDtlBpsC4uRIA==",rememberMe=TRUE)

staging_challenge_project = "syn20366916"

base_folder_ent = Folder("test_submission", parentId=staging_challenge_project)
stored_base_folder_ent = synStore(base_folder_ent)


template_data_ent = File("./dry_run/aim2_model_predictions_test.csv", parentId=stored_base_folder_ent$properties$id)
synStore(template_data_ent)

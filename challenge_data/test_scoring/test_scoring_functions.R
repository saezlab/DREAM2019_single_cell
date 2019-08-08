# test scoring functions:
# we just run some random data to test for syntactical mistakes.

# use the prepare_data_4_test.R to generate dummy data. 
library(dplyr)
library(tidyr)
### AIM 1.1 --------------------------------------------------------------------
source("./scoring_scripts/score_sc1.R")
source("./scoring_scripts/validate_sc1.R")
test_aim11_file  = "./challenge_data/test_scoring/test_data/aim_11_test_data.csv"
validation_file_aim11 = "./challenge_data/validation_data/AIM_11_data.csv"

sc1_validation = validate_sc1(prediction_data_file = test_aim11_file,
						validation_data_file = validation_file_aim11)

sc1_score = score_sc1(prediction_data_file = test_aim11_file,
							validation_data_file = validation_file_aim11)
# > aim11_score
# [1] 3.99223


### AIM 1.2.1   ---------------------------------------------------------------
source("./scoring_scripts/score_sc2.R")
source("./scoring_scripts/validate_sc2.R")
test_aim121_file  = "./challenge_data/test_scoring/test_data/aim_121_test_data.csv"
validation_file_aim121 = "./challenge_data/validation_data/AIM_121_data.csv"

sc2_validation = validate_sc2(prediction_data_file = test_aim121_file,
					  validation_data_file = validation_file_aim121)

sc2_score = score_sc2(prediction_data_file = test_aim121_file,
							 validation_data_file = validation_file_aim121)
# > aim121_score
# test
# prediction 179.614


### AIM 1.2.2   ---------------------------------------------------------------

source("./scoring_scripts/score_sc3.R")
source("./scoring_scripts/validate_sc3.R")
test_aim122_file  = "./challenge_data/test_scoring/test_data/aim_122_test_data.csv"
validation_file_aim122 = "./challenge_data/validation_data/AIM_122_data.csv"

sc3_validation = validate_sc3(prediction_data_file = test_aim122_file,
						  validation_data_file = validation_file_aim122)


sc3_score = score_sc3(prediction_data_file = test_aim122_file,
							 validation_data_file = validation_file_aim122)
# > aim122_score
# test
# prediction 336.7389


### AIM 2   ---------------------------------------------------------------

source("./scoring_scripts/score_sc4.R")
source("./scoring_scripts/validate_sc4.R")
test_aim2_file  = "./challenge_data/test_scoring/test_data/aim_2_test_data.csv"
validation_file_aim2 = "./challenge_data/validation_data/AIM_2_median_data.csv"

sc4_validation = validate_sc4(prediction_data_file = test_aim2_file,
					   validation_data_file = validation_file_aim2)

sc4_score = score_sc4(prediction_data_file = test_aim2_file,
						 validation_data_file = validation_file_aim2)
# > aim122_score
# test
# prediction 336.7389


# test scoring functions:
# we just run some random data to test for syntactical mistakes.

# use the prepare_data_4_test.R to generate dummy data. 
library(dplyr)
library(tidyr)
### Subchallenge 1 --------------------------------------------------------------------
source("./scoring_scripts/score_sc1.R")
source("./scoring_scripts/validate_sc1.R")


test_sc1_file  = "./challenge_data/test_scoring/test_data/sc1_test_data_v1.csv"
validation_file_sc1 = "./challenge_data/validation_data/sc1gold.csv"

sc1_validation = validate_sc1(prediction_data_file = test_sc1_file,
						validation_data_file = validation_file_sc1)

sc1_score = score_sc1(prediction_data_file = test_sc1_file,
							validation_data_file = validation_file_sc1)
# > aim11_score
# [1] 3.992752



test_sc1_file  = "./challenge_data/test_scoring/test_data/sc1_test_data_v4.csv"
validation_file_sc1 = "./challenge_data/validation_data/sc1gold.csv"

sc1_validation = validate_sc1(prediction_data_file = test_sc1_file,
							  validation_data_file = validation_file_sc1)

sc1_score = score_sc1(prediction_data_file = test_sc1_file,
					  validation_data_file = validation_file_sc1)



### Subchallenge 2   ---------------------------------------------------------------
source("./scoring_scripts/score_sc2.R")
source("./scoring_scripts/validate_sc2.R")
test_sc2_file  = "./challenge_data/test_scoring/test_data/sc2_test_data_v1.csv"
validation_file_sc2 = "./challenge_data/validation_data/sc2gold.csv"

sc2_validation = validate_sc2(prediction_data_file = test_sc2_file,
					  validation_data_file = validation_file_sc2)

sc2_score = score_sc2(prediction_data_file = test_sc2_file,
							 validation_data_file = validation_file_sc2)
# > aim121_score
# test
# prediction 179.614


### Subchallenge 3   ---------------------------------------------------------------

source("./scoring_scripts/score_sc3.R")
source("./scoring_scripts/validate_sc3.R")
test_sc3_file  = "./challenge_data/test_scoring/test_data/sc3_test_data_v1.csv"
validation_file_sc3 = "./challenge_data/validation_data/sc3gold.csv"

sc3_validation = validate_sc3(prediction_data_file = test_sc3_file,
						  validation_data_file = validation_file_sc3)


sc3_score = score_sc3(prediction_data_file = test_sc3_file,
							 validation_data_file = validation_file_sc3)
# > aim122_score
# test
# prediction 336.7389


### Subchallenge 4   ---------------------------------------------------------------

source("./scoring_scripts/score_sc4.R")
source("./scoring_scripts/validate_sc4.R")
test_sc4_file  = "./challenge_data/test_scoring/test_data/sc4_test_data_v1.csv"
validation_file_sc4 = "./challenge_data/validation_data/sc4gold.csv"

sc4_validation = validate_sc4(prediction_data_file = test_sc4_file,
					   validation_data_file = validation_file_sc4)

sc4_score = score_sc4(prediction_data_file = test_sc4_file,
						 validation_data_file = validation_file_sc4)
# > aim122_score
# test
# prediction 2.952615


# test with wrong files
source("./scoring_scripts/validate_sc4.R")
validation_file_sc4 = "./challenge_data/validation_data/sc4gold.csv"
validate_sc4(prediction_data_file = "./challenge_data/test_scoring/test_data/sc4_test_data_v1.csv",
			 validation_data_file = validation_file_sc4)

validate_sc4(prediction_data_file = "./challenge_data/test_scoring/test_data/sc4_test_data_v2_wrong_cols.csv",
			 validation_data_file = validation_file_sc4)

validate_sc4(prediction_data_file = "./challenge_data/test_scoring/test_data/sc4_test_data_v4_wrong_delim.tsv",
			 validation_data_file = validation_file_sc4)

validate_sc4(prediction_data_file = "./challenge_data/test_scoring/test_data/sc4_test_data_v5_no_header.csv",
			 validation_data_file = validation_file_sc4)




library(dplyr)
library(tidyr)
library(readr)

# scores the subchallenge 1
#' @param prediction_data_file path to prediction data file (.csv)
#' @param validation_data_file path to validation data file (.csv)
#' @description computes root-mean square error by conditions, then averages these

score_sc1 <- function(prediction_data_file, validation_data_file) {

  # load validation data

  validation_data <- read_csv(validation_data_file)
  prediction_data <- read_csv(prediction_data_file)

  required_columns <- c("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID", "p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")

  prediction_data <- prediction_data %>% select(required_columns)

  ### Formating -------------------------
  # convert to long format
  reporters <- c("p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")
  validation_data_long <- validation_data %>% gather(key = "marker", value = "test", reporters)
  prediction_data_long <- prediction_data %>% gather(key = "marker", value = "prediction", reporters)

  combined_data <- validation_data_long %>%
    full_join(prediction_data_long, by = c("glob_cellID","cell_line", "treatment", "time", "cellID", "fileID", "marker"))

  ### Calculate score --------------------
  # calculate the RMSE for each condition
  RMSE_cond <- combined_data %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(RMSE = sqrt(sum((test - prediction)^2) / n()))
  final_score <- mean(RMSE_cond$RMSE)
  return(final_score)
}

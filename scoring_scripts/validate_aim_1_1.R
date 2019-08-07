require(tidyverse)

# scores the subchallenge aim 1.1
#' @param prediction_data_file path to prediction data file (.csv)
#' @param validation_data_file path to validation data file (.csv)
#' @description checks input for missing columns
#' check input for missing conditions (missing predicted cells)
#' computes root-mean square error by conditions, then averages these

score_aim_1_1 <- function(prediction_data_file, validation_data_file) {

  # load validation data

  validation_data <- read_csv(validation_data_file)
  prediction_data <- read_csv(prediction_data_file)

  ### Checking inputs -------------------
  # checking columns of input data table
  required_columns <- c("cell_line", "treatment", "time", "cellID", "fileID", "p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")

  if (!all(required_columns %in% names(prediction_data))) {
    stop(paste0("missing columns detected. Required columns: ", paste(required_columns, collapse = ", ")))
  }
  prediction_data <- prediction_data %>% select(required_columns)


  # checking for any missing cell-predictions
  missing_prediction_data <- anti_join(validation_data, prediction_data, by = c("cell_line", "treatment", "time", "cellID", "fileID"))
  if (nrow(missing_prediction_data) > 0) {
    print(missing_prediction_data %>% select(c("cell_line", "treatment", "time", "cellID", "fileID")))
    stop("missing predictions detected for above conditions")
  }

  ### Formating -------------------------
  # convert to long format
  reporters <- c("p.Akt.Ser473.", "p.ERK", "p.HER2", "p.PLCg2", "p.S6")
  validation_data_long <- validation_data %>% gather(key = "marker", value = "test", reporters)
  prediction_data_long <- prediction_data %>% gather(key = "marker", value = "prediction", reporters)

  combined_data <- validation_data_long %>%
    full_join(prediction_data_long, by = c("cell_line", "treatment", "time", "cellID", "fileID", "marker"))

  ### Calculate score --------------------
  # calculate the RMSE for each condition
  RMSE_cond <- combined_data %>%
    group_by(cell_line, treatment, time, marker) %>%
    summarise(RMSE = sqrt(sum((test - prediction)^2) / n()))
  final_score <- mean(RMSE_cond$RMSE)
}

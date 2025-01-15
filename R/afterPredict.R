#' afterPredict
#'
#' This function return predicts results of a caret trained model and a test dataset and
#' optionally saves the results and model to specified file paths.
#'
#' @param omic_test Data frame or matrix containing the test dataset.
#' @param omic_test_diagnosis Independent test dataset variable.
#' @param seed Seed for reproduciblity.
#' @param dir_path Directory where results should be saved.
#' @param file_suffix String suffix for naming result files.
#' @param model Caret trained model.
#' @param model_suffix String prefix for naming result files.
#' @param save_results Logical. Should the confusion matrix and model be saved? Default is TRUE.
#' @return A list containing the confusion matrix and predictions.
#' @export
afterPredict <- function(omic_test, omic_test_diagnosis, seed, dir_path, file_suffix, model, model_suffix, save_results = TRUE) {

  set.seed(seed)

  # Validate inputs (basic examples)
  if (!is.data.frame(omic_test) && !is.matrix(omic_test)) {
    stop("omic_test must be a data frame or matrix.")
  }
  if (!is.vector(omic_test_diagnosis)) {
    stop("omic_test_diagnosis must be a vector.")
  }

  # Predict on test set
  pred <- predict(model, omic_test)
  actual <- as.factor(omic_test_diagnosis)

  # Handle potential mismatches in levels
  if (!all(levels(pred) %in% levels(actual))) {
    warning("Mismatch in levels between predictions and actual outcomes.")
  }

  # Create confusion matrix
  conf_matrix <- confusionMatrix(data = pred, reference = actual)

  # Save results if required
  if (save_results) {
    tryCatch({
      saveRDS(conf_matrix, file.path(dir_path, paste0(model_suffix, "_", file_suffix, "_cf.rds")))
      saveRDS(model, file.path(dir_path, paste0(model_suffix, "_", file_suffix, "_model.rds")))
    }, error = function(e) {
      warning("Failed to save results: ", e$message)
    })
  }

  return(list(confusion_matrix = conf_matrix, predictions = pred))
}

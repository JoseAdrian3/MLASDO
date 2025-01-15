#' Train an SVM Model with Radial Kernel
#'
#' This function trains a Support Vector Machine (SVM) model with a radial kernel.
#' It evaluates the model using Balanced Accuracy (BA) and computes confidence intervals.
#' The trained model, confusion matrix, and confidence intervals are saved to specified files.
#'
#' @param omic_train Transcriptomic train dataframe with features (independent variables, x).
#' @param omic_train_diagnosis Vector of training data labels (dependent variables, y).
#' @param omic_test Transcriptomic test dataframe with features (independent variables, x).
#' @param omic_test_diagnosis Vector of test data labels (diagnoses).
#' @param cost Numeric value for the cost parameter of the SVM.
#' @param seed Seed for reproducibility.
#' @param dir_path Directory path where output files will be saved.
#' @param file_suffix Suffix to append to saved files.
#' @param n_cores Number of CPU cores to use for parallel processing.
#'
#' @return A list containing the trained SVM model, confusion matrix, balanced accuracy, and confidence intervals.
#' @export
train_svm_radial <- function(omic_train, omic_train_diagnosis, omic_test, omic_test_diagnosis, cost, seed, dir_path, file_suffix, n_cores) {

  # Seed for reproducibility
  set.seed(seed)

  # Custom function to calculate Balanced Accuracy (BA)
  calc_balanced_accuracy <- function(pred, actual) {
    sensitivity <- sum(pred == "Case" & actual == "Case") / sum(actual == "Case")
    specificity <- sum(pred == "Control" & actual == "Control") / sum(actual == "Control")
    return((sensitivity + specificity) / 2)
  }

  # Estimate the sigma parameter for the radial kernel
  sigma_est <- kernlab::sigest(as.matrix(omic_train))[[1]]

  # Train the SVM model with a radial kernel
  svm_model <- e1071::svm(
    x = omic_train,
    y = as.factor(omic_train_diagnosis),
    type = "C-classification",
    kernel = "radial",
    cost = cost,
    gamma = 1 / (2 * sigma_est^2),
    scale = TRUE
  )

  # Make predictions on the test data
  pred <- predict(svm_model, omic_test)
  actual <- as.factor(omic_test_diagnosis)

  # Create a confusion matrix
  confusion <- table(Predicted = pred, Actual = actual)

  # Calculate Balanced Accuracy
  ba <- calc_balanced_accuracy(pred, actual)

  # Generate confidence intervals for Balanced Accuracy using bootstrap
  resamples <- replicate(1000, {
    idx <- sample(seq_along(actual), replace = TRUE)
    calc_balanced_accuracy(pred[idx], actual[idx])
  })
  ci <- quantile(resamples, probs = c(0.025, 0.975))

  # Save the confusion matrix, trained model, and confidence intervals
  saveRDS(confusion, file.path(dir_path, paste0("svm_radial_", file_suffix, "_cf.rds")))
  saveRDS(svm_model, file.path(dir_path, paste0("svm_radial_", file_suffix, "_model.rds")))
  saveRDS(ci, file.path(dir_path, paste0("svm_radial_", file_suffix, "_ci.rds")))

  # Return the results as a list
  return(list(
    model = svm_model,
    confusion_matrix = confusion,
    balanced_accuracy = ba,
    confidence_interval = ci
  ))
}

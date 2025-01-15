#' Train a LASSO Regression Model
#'
#' This function trains a LASSO regression model using cross-validation and calculates performance metrics.
#' The trained model, confusion matrix, and confidence intervals are saved to specified files.
#'
#' @param omic_train Transcriptomic train dataframe with features (independent variables, x).
#' @param omic_train_diagnosis Vector of training data labels (dependent variables, y).
#' @param omic_test Transcriptomic test dataframe with features (independent variables, x).
#' @param omic_test_diagnosis Vector of test data labels (diagnoses).
#' @param lambdas Vector of lambda values for tuning the LASSO model.
#' @param seed Seed for reproducibility.
#' @param dir_path Directory path where output files will be saved.
#' @param file_suffix Suffix to append to saved files.
#' @param n_cores Number of CPU cores to use for parallel processing.
#'
#' @return A trained LASSO regression model.
#' @export
train_lasso <- function(omic_train, omic_train_diagnosis, omic_test, omic_test_diagnosis, lambdas, seed, dir_path, file_suffix, n_cores) {

  # Seed for reproducibility
  set.seed(seed)

  # Custom metric function to calculate Balanced Accuracy (BA)
  metrics <- function(data, lev = levels(as.factor(data$obs)), model = NULL) {
    c(
      ba = (sensitivity(data[, "pred"], data[, "obs"], positive = "Case", negative = "Control") +
              specificity(data[, "pred"], data[, "obs"], positive = "Case", negative = "Control")) / 2
    )
  }

  # Define training control parameters with cross-validation
  train_control <- trainControl(
    method = "cv",
    number = 10,
    savePredictions = "all",
    classProbs = TRUE,
    sampling = "down",
    returnResamp = "all",
    allowParallel = TRUE,
    summaryFunction = metrics
  )

  # Initialize a cluster for parallel processing
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # Train the LASSO model using glmnet
  lasso_model <- train(
    x = as.matrix(omic_train),
    y = as.factor(omic_train_diagnosis),
    method = "glmnet",
    trControl = train_control,
    tuneGrid = expand.grid(alpha = 1, lambda = lambdas),
    family = "binomial",
    seed = seed,
    metric = "ba"
  )

  # Stop the cluster after training
  stopCluster(cl)

  # Extract resampling results and calculate confidence intervals for Balanced Accuracy
  resamples <- lasso_model$resample
  ci <- t.test(resamples$ba)$conf.int

  # Predict on the test set using the best lambda
  model_prediction <- predict(lasso_model, omic_test)
  actual <- as.factor(omic_test_diagnosis)

  # Create a confusion matrix to evaluate test set predictions
  lasso_conf_matrix <- confusionMatrix(model_prediction, actual)

  # Save the confusion matrix, trained model, and confidence intervals
  saveRDS(lasso_conf_matrix, file.path(dir_path, paste0("lasso_", file_suffix, "_cf.rds")))
  saveRDS(lasso_model, file.path(dir_path, paste0("lasso_", file_suffix, "_model.rds")))
  saveRDS(ci, file.path(dir_path, paste0("lasso_", file_suffix, "_ci.rds")))

  # Return the trained model
  return(lasso_model)
}

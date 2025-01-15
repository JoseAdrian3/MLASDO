#' Train a Random Forest Model
#'
#' This function trains a Random Forest (RF) model using cross-validation and calculates performance metrics.
#' The trained model, confusion matrix, and confidence intervals are saved to specified files.
#'
#' @param omic_train Transcriptomic train dataframe with features (independent variables, x).
#' @param omic_train_diagnosis Vector of training data labels (dependent variables, y).
#' @param omic_test Transcriptomic test dataframe with features (independent variables, x).
#' @param omic_test_diagnosis Vector of test data labels (diagnoses).
#' @param num_trees Integer indicating the number of trees to grow in the random forest.
#' @param mtry Integer specifying the number of variables randomly sampled as candidates at each split.
#' @param seed Seed for reproducibility.
#' @param dir_path Directory path where output files will be saved.
#' @param file_suffix Suffix to append to saved files.
#' @param n_cores Number of CPU cores to use for parallel processing.
#'
#' @return A trained Random Forest model.
#' @export
train_rf <- function(omic_train, omic_train_diagnosis, omic_test, omic_test_diagnosis, num_trees, mtry, seed, dir_path, file_suffix, n_cores) {

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

  # Train the Random Forest model
  rf_model <- train(
    omic_train,
    as.factor(omic_train_diagnosis),
    method = "rf",
    trControl = train_control,
    tuneGrid = data.frame(mtry = mtry),
    ntree = num_trees,
    importance = TRUE,
    classification = TRUE,
    seed = seed,
    metric = "ba"
  )

  # Stop the cluster after training
  stopCluster(cl)

  # Extract resampling results and calculate confidence intervals for Balanced Accuracy
  resamples <- rf_model$resample
  ci <- t.test(resamples$ba)$conf.int

  # Predict on the test set
  pred <- predict(rf_model, omic_test)
  actual <- as.factor(omic_test_diagnosis)

  # Create a confusion matrix to evaluate test set predictions
  rf_conf_matrix <- confusionMatrix(data = pred, reference = actual)

  # Save the confusion matrix, trained model, and confidence intervals
  saveRDS(rf_conf_matrix, file.path(dir_path, paste0("rf_", file_suffix, "_cf.rds")))
  saveRDS(rf_model, file.path(dir_path, paste0("rf_", file_suffix, "_model.rds")))
  saveRDS(ci, file.path(dir_path, paste0("rf_", file_suffix, "_ci.rds")))

  # Return the trained model
  return(rf_model)
}

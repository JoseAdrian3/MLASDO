#' Execute Support Vector Machine (SVM) Model
#'
#' This function trains and evaluates an SVM model with hyperparameter tuning using cross-validation.
#' It identifies misclassified samples and selects some based on their misclassification distance.
#'
#' @param saving_name Directory path where results will be saved.
#' @param omic_data Data frame containing omics data.
#' @param class_variable Name of the column with class labels.
#' @param id_column Name of the column containing unique IDs for samples.
#' @param n_cores Number of CPU cores to use for parallel processing.
#' @param seed Random seed for reproducibility.
#' @param c List of cost values for the SVM model.
#' @param kernel List of kernels to test (e.g., "linear", "radial").
#' @param gamma List of gamma values for radial kernel.
#' @param diagnostic_change_probability Proportion of samples to select for diagnostic change.
#' @param max_iterations Maximum number of iterations for the optimization.
#' @param tolerance Tolerance level for the optimization.
#'
#' @return A list containing the best hyperparameters, balanced accuracy, and misclassified samples.
#' @export
execute_svm <- function(
    saving_name,
    omic_data,
    class_variable,
    id_column,
    n_cores,
    seed,
    c,
    kernel,
    gamma,
    diagnostic_change_probability,
    max_iterations,
    tolerance
) {
  # Set random seed for reproducibility
  set.seed(seed)

  # Prepare the feature matrix (x) and target vector (y)
  x <- omic_data %>% dplyr::select(-c(id_column, class_variable))
  y <- as.factor(omic_data[[class_variable]])

  # Initialize variables to store the best model and results
  best_model <- NULL
  best_balanced_acc <- -Inf
  best_c <- NULL
  best_kernel <- NULL
  best_gamma <- NULL

  # Create cross-validation folds
  folds <- createFolds(y, k = 10)

  # Data frames to store misclassified samples
  misclassified_data_best <- data.frame()
  misclassified_data_all <- data.frame()

  # Set up parallel processing
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # Loop over all combinations of hyperparameters
  for (c_value in c) {
    for (kernel_value in kernel) {
      for (gamma_value in gamma) {

        fold_balanced_acc <- c()
        misclassified_data_fold <- data.frame()

        for (i in seq_along(folds)) {
          cat(i)

          test_indices <- folds[[i]]
          train_indices <- setdiff(seq_len(nrow(x)), test_indices)

          x_train <- x[train_indices, ]
          x_test <- x[test_indices, ]
          y_train <- y[train_indices]
          y_test <- y[test_indices]

          # Balance the training data using down-sampling
          train_data <- downSample(x = x_train, y = y_train)
          x_train_balanced <- train_data[, -ncol(train_data)]
          y_train_balanced <- train_data$Class

          # Train the SVM model based on the kernel type
          if (kernel_value == "radial") {
            set.seed(seed)
            model <- svm(
              x = x_train_balanced,
              y = y_train_balanced,
              kernel = kernel_value,
              cost = c_value,
              gamma = gamma_value,
              probability = TRUE
            )
          } else if (kernel_value == "linear") {
            set.seed(seed)
            model <- svm(
              x = x_train_balanced,
              y = y_train_balanced,
              kernel = kernel_value,
              cost = c_value,
              probability = TRUE
            )
          }

          # Make predictions on the test set
          predictions <- predict(model, x_test)

          # Calculate balanced accuracy
          confusion_matrix <- confusionMatrix(predictions, y_test)
          balanced_acc <- confusion_matrix$byClass["Balanced Accuracy"]
          fold_balanced_acc <- c(fold_balanced_acc, balanced_acc)

          # Identify misclassified samples
          misclassified_indices <- which(predictions != y_test)

          distances <- rep(NA, length(misclassified_indices))

          if (length(misclassified_indices) > 0) {
            support_vectors <- model$SV
            coefs <- model$coefs
            intercept <- model$rho

            w <- t(coefs) %*% support_vectors

            if (kernel_value == "linear") {
              distances <- sapply(misclassified_indices, function(idx) {
                x_i <- as.numeric(x_test[idx, ])
                abs(sum(w * x_i) + intercept) / sqrt(sum(w^2))
              })
            } else if (kernel_value == "radial") {
              distances <- sapply(misclassified_indices, function(idx) {
                x_i <- as.numeric(x_test[idx, ])

                K <- function(x, y) exp(-gamma_value * sum((x - y)^2))

                decision_value <- sum(coefs * sapply(1:nrow(support_vectors), function(k) K(support_vectors[k, ], x_i)) + intercept)

                abs(decision_value)
              })
            }

            misclassified_fold <- data.frame(
              index = test_indices[misclassified_indices],
              distance = distances,
              c = c_value,
              kernel = kernel_value,
              gamma = gamma_value,
              balanced_accuracy = balanced_acc
            )

            misclassified_data_fold <- rbind(misclassified_data_fold, misclassified_fold)
          }
        }

        avg_balanced_acc <- mean(fold_balanced_acc)

        if (nrow(misclassified_data_fold) > 0) {
          misclassified_data_fold$mean_balanced_accuracy <- avg_balanced_acc
          misclassified_data_all <- rbind(misclassified_data_all, misclassified_data_fold)
        }

        if (avg_balanced_acc > best_balanced_acc) {
          best_balanced_acc <- avg_balanced_acc
          best_c <- c_value
          best_kernel <- kernel_value
          best_gamma <- gamma_value

          misclassified_data_best <- misclassified_data_fold
        }
      }
    }
  }

  stopCluster(cl)

  # Save the SVM solution
  dir_path <- paste(saving_name, "svm_data", sep = "/")

  # Save final results
  final_results <- list(
    best_c = best_c,
    best_kernel = best_kernel,
    best_gamma = best_gamma,
    best_balanced_acc = best_balanced_acc,
    selected_samples = misclassified_data_best
  )

  results_path <- paste(saving_name, "svm_data.rds", sep = "_")
  saveRDS(final_results, file = file.path(dir_path, results_path))

  results_path <- paste(saving_name, "svm_data_all.rds", sep = "_")
  saveRDS(misclassified_data_all, file = file.path(dir_path, results_path))
}

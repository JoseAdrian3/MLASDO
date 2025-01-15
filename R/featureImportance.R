#' Feature Importance Calculation
#'
#' This function calculates the feature importance using Random Forest models.
#' of control2case or case2control changes.
#'
#' @param omic_data Omics dataframe.
#' @param model_dir String of the directory to save the model.
#' @param solution Vector indicating the solution for transformation of control2case or case2control.
#' @param mtry Number of variables randomly sampled at each split in Random Forest models.
#' @param filename_suffix String of the suffix for saved file names.
#' @param saving_name Prefix for saved file names.
#' @param class_variable String of the class variable.
#' @param partition_percentage Percentage of data to partition for training.
#' @param id_column String of the IDs variable.
#' @param n_cores Number of cores for parallel processing.
#' @param seed Seed for reproducibility.
#' @param bool String indicating transformation type ("Control2Case" or "Case2Control").
#' @param num_trees Number of trees for Random Forest models.
#'
#' @return Dataframe with sorted feature importance.
#' @export
feature_importance <- function(omic_data, model_dir, solution, mtry, filename_suffix, saving_name, class_variable, partition_percentage, id_column, n_cores, seed, bool, num_trees) {

  # Seed for reproducibility
  set.seed(seed)

  # Create directory for saving results
  dir_path <- paste(saving_name, "importance", sep = "/")
  dir.create(dir_path)

  # Filter data based on the transformation type (Control2Case or Case2Control)
  if (bool == "Control2Case") {
    control_idx <- omic_data[, class_variable] == "Control"
    omic_data_filter <- omic_data[control_idx, ]
    solution_control <- solution[control_idx]
    omic_data_filter[, class_variable] <- solution_control
  } else if (bool == "Case2Control") {
    case_idx <- omic_data[, class_variable] == "Case"
    omic_data_filter <- omic_data[case_idx, ]
    solution_case <- solution[case_idx]
    omic_data_filter[, class_variable] <- solution_case
  }

  # Extract predictor (X) and response (Y) variables
  X <- omic_data_filter[, !names(omic_data_filter) %in% c(class_variable, id_column)]
  Y <- as.factor(omic_data_filter[, class_variable])

  # Proceed only if Y has two levels and is properly labeled
  if (length(levels(Y)) == 2 && sum(as.numeric(as.character(Y)) == 1) > 1) {

    # Relabel levels for easier interpretation
    levels(Y) <- c("maintained", "changed")

    # Custom metric for model evaluation: Balanced Accuracy
    metrics <- function(data, lev = levels(as.factor(data$obs)), model = NULL) {
      c(
        BA = (sensitivity(data[, "pred"], data[, "obs"], positive = "changed") +
                specificity(data[, "pred"], data[, "obs"], negative = "maintained")) / 2
      )
    }

    # Cross-validation setup
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

    # Hyperparameter grid for tuning
    my_grid <- expand.grid(
      mtry = c(mtry, mtry * 2, mtry * 3, mtry * 4),
      splitrule = c("gini", "hellinger", "extratrees"),
      min.node.size = c(1, 2, 4, 8)
    )

    # Define tree sizes to evaluate
    tree_sizes <- c(125, 250, 500)
    rf_cv_results <- list()

    # Enable parallel processing
    cluster <- makeCluster(n_cores)
    registerDoParallel(cluster)

    # Train Random Forest models with different tree sizes
    for (tree in tree_sizes) {
      rf_model <- train(
        X,
        Y,
        method = "ranger",
        trControl = train_control,
        tuneGrid = my_grid,
        num.trees = tree,
        classification = TRUE,
        seed = seed,
        metric = "BA",
        importance = "permutation"
      )
      rf_cv_results[[paste(tree, "trees")]] <- rf_model
    }

    # Stop parallel processing
    stopCluster(cluster)

    # Combine results from all tree sizes
    combined_results <- rbind(
      data.frame(mtry = rf_cv_results$`125 trees`$results$mtry, splitrule = rf_cv_results$`125 trees`$results$splitrule, min.node.size = rf_cv_results$`125 trees`$results$min.node.size, MetricValue = rf_cv_results$`125 trees`$results$BA, Trees = "125 trees", Metric = "Balanced Accuracy"),
      data.frame(mtry = rf_cv_results$`250 trees`$results$mtry, splitrule = rf_cv_results$`250 trees`$results$splitrule, min.node.size = rf_cv_results$`250 trees`$results$min.node.size, MetricValue = rf_cv_results$`250 trees`$results$BA, Trees = "250 trees", Metric = "Balanced Accuracy"),
      data.frame(mtry = rf_cv_results$`500 trees`$results$mtry, splitrule = rf_cv_results$`500 trees`$results$splitrule, min.node.size = rf_cv_results$`500 trees`$results$min.node.size, MetricValue = rf_cv_results$`500 trees`$results$BA, Trees = "500 trees", Metric = "Balanced Accuracy")
    )

    # Determine the best tree configuration
    best_model_tree <- aggregate(MetricValue ~ Trees, data = combined_results, FUN = max)
    best_config <- best_model_tree[which.max(best_model_tree$MetricValue), "Trees"]
    best_rf_model <- rf_cv_results[[best_config]]

    # Predict and calculate confusion matrix
    predictions <- predict(best_rf_model, X)
    confusionMatrix(predictions, Y, positive = "changed")

    # Calculate feature importance
    importance_df <- varImp(best_rf_model, scale = TRUE)$importance
    importance_df <- as.data.frame(importance_df)
    importance_df$variable <- rownames(importance_df)
    names(importance_df) <- c("importance", "variable")
    sorted_importance_df <- importance_df[order(-importance_df$importance), ]
  } else {
    # If criteria not met, return NULL for importance and model
    sorted_importance_df <- NULL
    best_rf_model <- NULL
  }

  # Save results to files
  saveRDS(sorted_importance_df, file.path(dir_path, paste0(saving_name, "_", filename_suffix, "_", bool, ".rds")))
  saveRDS(best_rf_model, file.path(dir_path, paste0(saving_name, "_", filename_suffix, "_", bool, "_model.rds")))

  return(sorted_importance_df)
}

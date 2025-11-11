#' Run the Random Forest feature selection step
.mlasdo_run_rf_feature_selection <- function(data,
                                             n_cores,
                                             clinical_data = NULL,
                                             id_column,
                                             target_column,
                                             always_split_variables,
                                             rf_mtry,
                                             rf_splitrule,
                                             rf_min_node_size,
                                             train_method,
                                             train_number,
                                             train_repeats,
                                             rf_num_trees,
                                             rf_model_path,
                                             seed,
                                             importance,
                                             verbose) {
  required_pkgs <- c("caret", "ranger", "doParallel")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' must be installed to run feature selection.", pkg))
    }
  }

  set.seed(seed)

  if (!target_column %in% names(data)) {
    stop(sprintf("Column '%s' not present in the omic dataset.", target_column))
  }

  if (!is.null(clinical_data)) {
    .mlasdo_require_columns(clinical_data, id_column, "clinical dataset")
    if (!identical(data[[id_column]], clinical_data[[id_column]])) {
      stop("Clinical data is not aligned with omic data; IDs differ.")
    }
  }

  if (is.null(always_split_variables)) {
    always_split_variables <- character()
  }

  if (!is.null(clinical_data) && length(always_split_variables) > 0) {
    clinical_vars <- intersect(always_split_variables, names(clinical_data))
    for (var in clinical_vars) {
      data[[var]] <- clinical_data[[var]]
    }
  }

  predictors <- setdiff(names(data), c(id_column, target_column))
  if (length(predictors) == 0) {
    stop("No predictors available for feature selection after removing ID and target columns.")
  }

  complete_idx <- !is.na(data[[target_column]])
  if (!all(complete_idx)) {
    data <- data[complete_idx, , drop = FALSE]
    if (isTRUE(verbose)) {
      message(sprintf("Removed %d rows with NA in '%s' before training.",
                      sum(!complete_idx), target_column))
    }
  }

  y <- factor(data[[target_column]])
  x <- data[, predictors, drop = FALSE]

  if (is.null(rf_mtry)) {
    mtry_values <- .mlasdo_default_mtry(ncol(x))
  } else {
    mtry_values <- unique(rf_mtry)
  }
  if (length(mtry_values) == 0) {
    stop("No valid mtry values provided.")
  }

  splitrule_values <- unique(rf_splitrule)
  if (length(splitrule_values) == 0) {
    stop("At least one splitrule value must be provided.")
  }
  min_node_values <- unique(rf_min_node_size)
  if (length(min_node_values) == 0) {
    stop("At least one min.node.size value must be provided.")
  }

  rf_grid <- expand.grid(
    mtry = mtry_values,
    splitrule = splitrule_values,
    min.node.size = min_node_values
  )

  train_method <- match.arg(train_method, choices = c("repeatedcv", "cv"))
  ctrl <- .mlasdo_build_train_control(
    method = train_method,
    number = train_number,
    repeats = train_repeats
  )

  available_split_vars <- intersect(always_split_variables, names(x))
  if (length(always_split_variables) > 0 && length(available_split_vars) == 0 &&
      isTRUE(verbose)) {
    message("None of the requested always.split.variables are present in the omic dataset.")
  }
  if (length(available_split_vars) == 0) {
    available_split_vars <- NULL
  }

  n_cores <- max(1L, min(.mlasdo_select_cores(n_cores), parallel::detectCores() - 1L))
  cl <- parallel::makeCluster(n_cores)
  if ("MLASDO" %in% loadedNamespaces()) {
    parallel::clusterExport(
      cl,
      varlist = c(".mlasdo_compute_balanced_accuracy"),
      envir = asNamespace("MLASDO")
    )
  }
  doParallel::registerDoParallel(cl)

  if (!is.null(importance)) {
    importance <- match.arg(
      importance,
      choices = c("none", "impurity", "impurity_corrected", "permutation")
    )
  }

  rf_model <- caret::train(
    x = x,
    y = y,
    method = "ranger",
    trControl = ctrl,
    tuneGrid = rf_grid,
    importance = importance,
    always.split.variables = available_split_vars,
    num.trees = rf_num_trees,
    metric = "Balanced_Accuracy",
  )

  parallel::stopCluster(cl)
  foreach::registerDoSEQ()

  dir.create(dirname(rf_model_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(rf_model, rf_model_path)

  meta <- list(
    model_path = normalizePath(rf_model_path, winslash = "/", mustWork = FALSE),
    best_tune = rf_model$bestTune,
    grid = rf_grid,
    train_method = train_method,
    train_number = train_number,
    train_repeats = if (train_method == "repeatedcv") train_repeats else NA,
    always_split_variables = available_split_vars,
    num_trees = rf_num_trees,
    class_levels = levels(y),
    seed = seed,
    importance = importance
  )

  list(model = rf_model, metadata = meta)
}

.mlasdo_default_mtry <- function(predictor_count) {
  if (predictor_count <= 1) {
    return(1)
  }
  unique(pmax(1, floor(seq(2, sqrt(predictor_count), length.out = 3))))
}

.mlasdo_build_train_control <- function(method, number, repeats) {
  ctrl_args <- list(
    method = method,
    number = number,
    verboseIter = TRUE,
    summaryFunction = .mlasdo_bal_acc_summary,
    classProbs = TRUE,
    savePredictions = "final",
    sampling = .mlasdo_downsample_wrapper
  )
  if (identical(method, "repeatedcv")) {
    ctrl_args$repeats <- repeats
  }
  do.call(caret::trainControl, ctrl_args)
}

.mlasdo_downsample_wrapper <- function(x, y) {
  predictor_names <- colnames(x)
  down <- caret::downSample(
    x = as.data.frame(x),
    y = y,
    yname = ".mlasdo_class"
  )
  if (is.null(predictor_names)) {
    predictor_names <- setdiff(colnames(down), ".mlasdo_class")
  }
  list(
    x = down[, predictor_names, drop = FALSE],
    y = down[[".mlasdo_class"]]
  )
}

.mlasdo_bal_acc_summary <- function(data, lev = NULL, model = NULL) {
  obs <- factor(data$obs, levels = lev)
  pred <- factor(data$pred, levels = levels(obs))
  ba <- .mlasdo_compute_balanced_accuracy(pred, obs)
  stats::setNames(ba, "Balanced_Accuracy")
}

.mlasdo_compute_balanced_accuracy <- function(pred, obs) {
  if (!is.factor(obs)) {
    obs <- factor(obs)
  }
  pred <- factor(pred, levels = levels(obs))
  conf <- table(obs, pred)
  if (length(conf) == 0) {
    return(NA_real_)
  }
  per_class_recall <- diag(conf) / rowSums(conf)
  mean(per_class_recall, na.rm = TRUE)
}

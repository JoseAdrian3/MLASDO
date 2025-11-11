#' Train SVM models on the feature-selected MLASDO dataset
#'
#' @description
#' Executes the post-feature-selection modelling stage described in the SVM
#' scripts: the same 10-fold split is reused across linear, radial and
#' polynomial kernels, the best hyperparameters per kernel are learned via
#' cross-validation on the full dataset, and each fold is then retrained with
#' those hyperparameters to obtain predictions, distances to the separating
#' hyperplane and Balanced Accuracy metrics. All artifacts are stored under the
#' provided `output_dir`, including a summary with the kernel that achieved the
#' highest average Balanced Accuracy.
#'
#' @param omic_data Data frame or path to an `.rds` file with the trimmed omic
#'   dataset. Must contain `id_column` and `target_column`.
#' @param id_column Identifier column (defaults to `"participant_id"`).
#' @param target_column Target variable column (defaults to `"diagnosis_type"`).
#' @param output_dir Directory where all kernel-specific folders and the summary
#'   will be written.
#' @param kernels Character vector with any subset of `"linear"`, `"radial"`,
#'   and `"polynomial"`. Defaults to all three.
#' @param linear_grid,radial_grid,polynomial_grid Optional `expand.grid()` data
#'   frames with the hyperparameter combinations evaluated for each kernel.
#'   `NULL` triggers the default grids.
#' @param linear_cost_values Candidate `C` values for the linear kernel when
#'   `linear_grid` is `NULL`.
#' @param radial_cost_values Candidate `C` values for the radial kernel when
#'   `radial_grid` is `NULL`.
#' @param radial_sigma_values Candidate `sigma` (gamma) values for the radial
#'   kernel when `radial_grid` is `NULL`.
#' @param polynomial_cost_values Candidate `C` values for the polynomial kernel
#'   when `polynomial_grid` is `NULL`.
#' @param polynomial_degree_values Candidate polynomial degrees evaluated when
#'   `polynomial_grid` is `NULL`.
#' @param polynomial_scale_values Candidate `scale` (`gamma`) values for the
#'   polynomial kernel when `polynomial_grid` is `NULL`.
#' @param num_folds Number of folds used for CV (defaults to 10).
#' @param folds_seed Seed applied before generating the shared folds.
#' @param n_cores Number of worker processes used by `caret::train`.
#' @param sampling Strategy passed to `caret::trainControl(sampling = ...)`.
#' @param positive_class Label treated as the positive class (used to order the
#'   target factor and for confusion matrices).
#' @param negative_class Label treated as the negative class.
#'
#' @return A list with the per-kernel results, the summary table, and the best
#'   kernel metadata. Side effects include all the saved RDS/CSV files.
#' @examples
#' \dontrun{
#' mlasdo_train_svm_models(
#'   omic_data = "data/mlasdo_runs/feature_selection_only/feature_selection/omic_selected_features.rds",
#'   output_dir = "data/svm/feature_importance_shared_cv"
#' )
#' }
#' @export
mlasdo_train_svm_models <- function(omic_data,
                                    id_column = "participant_id",
                                    target_column = "diagnosis_type",
                                    output_dir,
                                    kernels = c("linear", "radial", "polynomial"),
                                    linear_grid = NULL,
                                    radial_grid = NULL,
                                    polynomial_grid = NULL,
                                    linear_cost_values = 10^seq(-6, 6, length.out = 13),
                                    radial_cost_values = 10^seq(-4, 4, length.out = 9),
                                    radial_sigma_values = 10^seq(-6, 1, length.out = 8),
                                    polynomial_cost_values = 10^seq(-4, 4, length.out = 9),
                                    polynomial_degree_values = 2:4,
                                    polynomial_scale_values = 10^seq(-6, 0, length.out = 7),
                                    num_folds = 10,
                                    folds_seed = 2025,
                                    n_cores = max(1L, parallel::detectCores() - 1L),
                                    sampling = "down",
                                    positive_class = NULL,
                                    negative_class = NULL) {
  .mlasdo_stop_if_missing(omic_data, "omic_data")
  .mlasdo_stop_if_missing(output_dir, "output_dir")

  required_pkgs <- c("caret", "doParallel", "e1071", "kernlab", "LiblineaR")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' must be installed to train the SVM models.", pkg))
    }
  }

  kernels <- match.arg(kernels, choices = c("linear", "radial", "polynomial"), several.ok = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  omic_df <- .mlasdo_as_data_frame(omic_data, "omic_data")
  .mlasdo_require_columns(omic_df, c(id_column, target_column), "omic dataset")

  target_cfg <- .mlasdo_configure_target(
    omic_df[[target_column]],
    positive_class = positive_class,
    negative_class = negative_class
  )
  omic_df[[target_column]] <- target_cfg$values
  resolved_positive <- target_cfg$positive
  resolved_negative <- target_cfg$negative
  if (is.null(resolved_positive) || is.null(resolved_negative)) {
    lvl_order <- levels(target_cfg$values)
    if (is.null(resolved_negative) && length(lvl_order) >= 1) {
      resolved_negative <- lvl_order[1]
    }
    if (is.null(resolved_positive) && length(lvl_order) >= 2) {
      resolved_positive <- lvl_order[2]
    }
  }

  complete_idx <- !is.na(omic_df[[target_column]])
  if (!all(complete_idx)) {
    omic_df <- omic_df[complete_idx, , drop = FALSE]
  }
  omic_df[[target_column]] <- droplevels(omic_df[[target_column]])

  predictors <- setdiff(names(omic_df), c(id_column, target_column))
  if (length(predictors) == 0) {
    stop("No predictors available after removing ID and target columns.")
  }

  x <- as.matrix(omic_df[, predictors, drop = FALSE])
  y <- factor(omic_df[[target_column]])

  set.seed(folds_seed)
  shared_folds <- caret::createFolds(y, k = num_folds, returnTrain = TRUE)
  saveRDS(shared_folds, file.path(output_dir, "shared_folds_indices.rds"))

  default_linear_grid <- expand.grid(C = linear_cost_values)
  default_radial_grid <- expand.grid(
    C = radial_cost_values,
    sigma = radial_sigma_values
  )
  default_poly_grid <- expand.grid(
    C = polynomial_cost_values,
    degree = polynomial_degree_values,
    scale = polynomial_scale_values
  )

  kernel_configs <- list()
  if ("linear" %in% kernels) {
    kernel_configs$linear <- list(
      caret_method = "svmLinear",
      tune_grid = if (is.null(linear_grid)) default_linear_grid else linear_grid,
      extra_args = list()
    )
  }
  if ("radial" %in% kernels) {
    kernel_configs$radial <- list(
      caret_method = "svmRadial",
      tune_grid = if (is.null(radial_grid)) default_radial_grid else radial_grid,
      extra_args = list(prob.model = TRUE)
    )
  }
  if ("polynomial" %in% kernels) {
    kernel_configs$polynomial <- list(
      caret_method = "svmPoly",
      tune_grid = if (is.null(polynomial_grid)) default_poly_grid else polynomial_grid,
      extra_args = list()
    )
  }

  sampling_fn <- sampling
  if (is.character(sampling_fn)) {
    if (length(sampling_fn) != 1L) {
      stop("`sampling` must be length 1 when specified as a character string.")
    }
    if (identical(sampling_fn, "down")) {
      sampling_fn <- .mlasdo_downsample_wrapper
    }
  }

  kernel_results <- list()
  for (kernel_name in names(kernel_configs)) {
    cfg <- kernel_configs[[kernel_name]]
    kernel_results[[kernel_name]] <- .mlasdo_run_svm_kernel(
      kernel_name = kernel_name,
      caret_method = cfg$caret_method,
      tune_grid = cfg$tune_grid,
      extra_train_args = cfg$extra_args,
      shared_folds = shared_folds,
      x = x,
      y = y,
      id_column = id_column,
      omic_data = omic_df,
      base_output_dir = output_dir,
      n_cores = n_cores,
      sampling = sampling_fn,
      positive_class = resolved_positive
    )
  }

  summary_df <- data.frame(
    kernel = names(kernel_results),
    mean_balanced_accuracy = vapply(
      kernel_results,
      function(res) res$mean_balanced_accuracy,
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  summary_path_rds <- file.path(output_dir, "balanced_accuracy_summary.rds")
  summary_path_csv <- file.path(output_dir, "balanced_accuracy_summary.csv")
  saveRDS(summary_df, summary_path_rds)
  utils::write.csv(summary_df, summary_path_csv, row.names = FALSE)

  best_idx <- which.max(summary_df$mean_balanced_accuracy)
  best_kernel <- summary_df$kernel[best_idx]
  best_model_dir <- file.path(output_dir, "best_model")
  dir.create(best_model_dir, recursive = TRUE, showWarnings = FALSE)

  best_info <- list(
    best_kernel = best_kernel,
    mean_balanced_accuracy = summary_df$mean_balanced_accuracy[best_idx],
    source_dir = kernel_results[[best_kernel]]$output_dir,
    hyperparameters = kernel_results[[best_kernel]]$best_hyperparameters,
    summary_table = summary_df,
    positive_class = resolved_positive,
    negative_class = resolved_negative
  )
  saveRDS(best_info, file.path(best_model_dir, "best_model_summary.rds"))

  best_artifacts <- kernel_results[[best_kernel]]$artifacts
  for (artifact_path in best_artifacts) {
    file.copy(
      from = artifact_path,
      to = file.path(best_model_dir, basename(artifact_path)),
      overwrite = TRUE
    )
  }

  list(
    summary = summary_df,
    kernel_results = kernel_results,
    best_kernel = best_kernel,
    best_model_dir = normalizePath(best_model_dir, winslash = "/", mustWork = FALSE),
    shared_folds = shared_folds,
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    class_roles = list(
      positive = resolved_positive,
      negative = resolved_negative
    )
  )
}

.mlasdo_load_existing_svm_results <- function(svm_dir) {
  svm_dir <- normalizePath(svm_dir, winslash = "/", mustWork = TRUE)
  summary_path <- file.path(svm_dir, "balanced_accuracy_summary.rds")
  if (!file.exists(summary_path)) {
    stop(sprintf("Could not find 'balanced_accuracy_summary.rds' under '%s'.", svm_dir))
  }
  summary_df <- readRDS(summary_path)
  kernel_results <- list()
  for (i in seq_len(nrow(summary_df))) {
    kernel_name <- summary_df$kernel[i]
    kernel_dir <- file.path(svm_dir, kernel_name)
    svm_data_dir <- file.path(kernel_dir, "svm_data")
    artifacts <- list(
      fold_models = file.path(svm_data_dir, "fold_models_best_params.rds"),
      distances = file.path(svm_data_dir, "distances_best_params_cv.rds"),
      metrics = file.path(svm_data_dir, "cv_best_params_metrics.rds"),
      fold_train = file.path(svm_data_dir, "fold_train_sets_best_params.rds"),
      fold_test = file.path(svm_data_dir, "fold_test_sets_best_params.rds"),
      confmats = file.path(svm_data_dir, "fold_confmats_best_params.rds"),
      hyperparameters = file.path(svm_data_dir, "best_hyperparameters.rds")
    )
    best_params <- NULL
    if (file.exists(artifacts$hyperparameters)) {
      loaded_params <- readRDS(artifacts$hyperparameters)
      if (is.data.frame(loaded_params)) {
        best_params <- as.list(loaded_params[1, , drop = FALSE])
      } else {
        best_params <- loaded_params
      }
    }
    kernel_results[[kernel_name]] <- list(
      kernel = kernel_name,
      mean_balanced_accuracy = summary_df$mean_balanced_accuracy[i],
      best_hyperparameters = best_params,
      output_dir = kernel_dir,
      artifacts = artifacts
    )
  }
  shared_folds <- NULL
  shared_path <- file.path(svm_dir, "shared_folds_indices.rds")
  if (file.exists(shared_path)) {
    shared_folds <- readRDS(shared_path)
  }
  best_model_dir <- file.path(svm_dir, "best_model")
  best_kernel <- summary_df$kernel[which.max(summary_df$mean_balanced_accuracy)]
  class_roles <- list(positive = NULL, negative = NULL)
  best_summary_path <- file.path(best_model_dir, "best_model_summary.rds")
  if (file.exists(best_summary_path)) {
    best_info <- readRDS(best_summary_path)
    if (!is.null(best_info$best_kernel)) {
      best_kernel <- best_info$best_kernel
    }
    if (!is.null(best_info$hyperparameters)) {
      kernel_results[[best_kernel]]$best_hyperparameters <- best_info$hyperparameters
    }
    if (!is.null(best_info$positive_class)) {
      class_roles$positive <- best_info$positive_class
    }
    if (!is.null(best_info$negative_class)) {
      class_roles$negative <- best_info$negative_class
    }
  }
  list(
    summary = summary_df,
    kernel_results = kernel_results,
    best_kernel = best_kernel,
    best_model_dir = normalizePath(best_model_dir, winslash = "/", mustWork = FALSE),
    shared_folds = shared_folds,
    output_dir = normalizePath(svm_dir, winslash = "/", mustWork = FALSE),
    class_roles = class_roles
  )
}

.mlasdo_build_seeds_list <- function(num_folds, grid_size) {
  seeds <- vector("list", num_folds + 1)
  for (i in seq_len(num_folds)) {
    seeds[[i]] <- sample.int(100000, grid_size)
  }
  seeds[[num_folds + 1]] <- sample.int(100000, 1)
  seeds
}

.mlasdo_run_svm_kernel <- function(kernel_name,
                                   caret_method,
                                   tune_grid,
                                   extra_train_args,
                                   shared_folds,
                                   x,
                                   y,
                                   id_column,
                                   omic_data,
                                   base_output_dir,
                                   n_cores,
                                   sampling,
                                   positive_class) {
  num_folds <- length(shared_folds)
  if (nrow(tune_grid) == 0) {
    stop(sprintf("The tuning grid for kernel '%s' is empty.", kernel_name))
  }
  .mlasdo_log_step("  - Training kernel %s.", kernel_name)

  seeds_list <- .mlasdo_build_seeds_list(num_folds, nrow(tune_grid))
  ctrl_cv <- caret::trainControl(
    method = "cv",
    number = num_folds,
    index = shared_folds,
    summaryFunction = .mlasdo_bal_acc_summary,
    classProbs = TRUE,
    savePredictions = "final",
    allowParallel = TRUE,
    sampling = sampling,
    seeds = seeds_list
  )

  train_args <- c(
    list(
      x = x,
      y = y,
      method = caret_method,
      tuneGrid = tune_grid,
      trControl = ctrl_cv,
      metric = "Balanced_Accuracy"
    ),
    extra_train_args
  )

  workers <- .mlasdo_select_cores(n_cores)
  cl <- parallel::makeCluster(workers)
  if ("MLASDO" %in% loadedNamespaces()) {
    parallel::clusterExport(
      cl,
      varlist = c(".mlasdo_compute_balanced_accuracy"),
      envir = asNamespace("MLASDO")
    )
  }
  doParallel::registerDoParallel(cl)
  fit_global <- tryCatch(
    do.call(caret::train, train_args),
    finally = {
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }
  )
  best_params <- fit_global$bestTune
  best_params_list <- as.list(best_params[1, , drop = FALSE])

  fold_BA <- numeric(num_folds)
  fold_models <- vector("list", num_folds)
  fold_confmats <- vector("list", num_folds)
  fold_train_sets <- vector("list", num_folds)
  fold_test_sets <- vector("list", num_folds)
  margins_list <- vector("list", num_folds)

  for (fold_id in seq_len(num_folds)) {
    tr_idx <- shared_folds[[fold_id]]
    te_idx <- setdiff(seq_len(nrow(x)), tr_idx)

    fold_args <- c(
      list(
        x = x[tr_idx, , drop = FALSE],
        y = y[tr_idx],
        method = caret_method,
        tuneGrid = best_params,
        trControl = caret::trainControl(method = "none")
      ),
      extra_train_args
    )
    fit_fold <- do.call(caret::train, fold_args)
  
    preds <- predict(fit_fold$finalModel, x[te_idx, , drop = FALSE], type = "response")
    margin <- predict(fit_fold$finalModel, x[te_idx, , drop = FALSE], type = "decision")

    pred_factor <- factor(preds, levels = levels(y))
    if (!is.null(positive_class) && positive_class %in% levels(y)) {
      confmat <- caret::confusionMatrix(pred_factor, y[te_idx], positive = positive_class)
    } else {
      confmat <- caret::confusionMatrix(pred_factor, y[te_idx])
    }
    ba <- .mlasdo_compute_balanced_accuracy(pred_factor, y[te_idx])

    fold_models[[fold_id]] <- fit_fold
    fold_confmats[[fold_id]] <- confmat
    fold_train_sets[[fold_id]] <- omic_data[tr_idx, ]
    fold_test_sets[[fold_id]] <- omic_data[te_idx, ]
    fold_BA[fold_id] <- ba

    fold_margins <- data.frame(
      id = omic_data[[id_column]][te_idx],
      fold = fold_id,
      diagnosis_true = y[te_idx],
      diagnosis_pred = pred_factor,
      signed_margin = as.numeric(margin),
      abs_margin = abs(as.numeric(margin)),
      kernel = kernel_name,
      stringsAsFactors = FALSE
    )
    for (param in names(best_params_list)) {
      fold_margins[[param]] <- best_params_list[[param]]
    }
    margins_list[[fold_id]] <- fold_margins
  }

  kernel_dir <- file.path(base_output_dir, kernel_name)
  out_dir <- file.path(kernel_dir, "svm_data")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  artifacts <- list(
    fold_models = file.path(out_dir, "fold_models_best_params.rds"),
    distances = file.path(out_dir, "distances_best_params_cv.rds"),
    metrics = file.path(out_dir, "cv_best_params_metrics.rds"),
    fold_train = file.path(out_dir, "fold_train_sets_best_params.rds"),
    fold_test = file.path(out_dir, "fold_test_sets_best_params.rds"),
    confmats = file.path(out_dir, "fold_confmats_best_params.rds"),
    hyperparameters = file.path(out_dir, "best_hyperparameters.rds")
  )

  saveRDS(fold_models, artifacts$fold_models)
  saveRDS(do.call(rbind, margins_list), artifacts$distances)
  saveRDS(
    data.frame(
      fold = seq_len(num_folds),
      Balanced_Accuracy = fold_BA,
      stringsAsFactors = FALSE
    ),
    artifacts$metrics
  )
  saveRDS(fold_train_sets, artifacts$fold_train)
  saveRDS(fold_test_sets, artifacts$fold_test)
  saveRDS(fold_confmats, artifacts$confmats)
  saveRDS(best_params, artifacts$hyperparameters)

  .mlasdo_log_step("  - Kernel %s completed.", kernel_name)
  list(
    kernel = kernel_name,
    mean_balanced_accuracy = mean(fold_BA),
    best_hyperparameters = best_params_list,
    output_dir = kernel_dir,
    artifacts = artifacts
  )
}

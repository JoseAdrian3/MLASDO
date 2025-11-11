#' Select top omic features based on a trained MLASDO RF model
#'
#' This helper loads (or receives) an already-trained `caret::train` object,
#' retrieves permutation importances, estimates the elbow point using
#' `inflection::uik()`, and trims the omic dataset to the most informative
#' variables—mirroring the workflow in
#' `scripts/feature_selection/feature_selection.Rmd`.
#'
#' @param model A fitted `caret` model or a path to an `.rds` file containing
#'   one (usually the RF produced by `run_mlasdo()`).
#' @param omic_data Data frame or `.rds` path with the full omic dataset.
#' @param id_column Identifier column preserved in the output.
#' @param target_column Target variable preserved in the output.
#' @param exclude_variables Optional character vector with covariates that should
#'   be excluded from the trimmed dataset (e.g., `"age_at_baseline"` and `"sex"`).
#'   These variables still appear in the importance ranking but are removed before
#'   saving the reduced omic table.
#' @param elbow_method Currently only `"uik"` is implemented; reserved for future
#'   strategies.
#' @param output_dir Optional directory where the trimmed dataset will be
#'   written. When provided, a default filename (`omic_selected_features.rds`) is
#'   used unless `output_path` overrides it.
#' @param output_path Optional explicit file path for the trimmed dataset. Takes
#'   precedence over `output_dir`. Directories are created automatically.
#' @param verbose Print short progress messages.
#'
#' @return A list with the sorted importances, selected feature names, and the
#'   trimmed omic data frame.
#' @export
mlasdo_select_top_features <- function(model,
                                       omic_data,
                                       id_column = "participant_id",
                                       target_column = "diagnosis_type",
                                       exclude_variables = c("age_at_baseline", "sex"),
                                       elbow_method = c("uik"),
                                       output_dir = NULL,
                                       output_path = NULL,
                                       verbose = interactive()) {
  if (missing(model)) stop("Argument 'model' is required.")
  if (missing(omic_data)) stop("Argument 'omic_data' is required.")

  elbow_method <- match.arg(elbow_method)

  rf_model <- .mlasdo_as_model_object(model)
  omic_df <- .mlasdo_as_data_frame(omic_data, "omic_data")
  .mlasdo_require_columns(omic_df, c(id_column, target_column), "omic dataset")

  imp <- caret::varImp(rf_model)
  if (!"importance" %in% names(imp)) {
    stop("Unexpected varImp structure; missing 'importance' data frame.")
  }
  importance_df <- imp$importance

  sorted_idx <- order(importance_df$Overall, decreasing = TRUE)
  imp_sorted <- importance_df$Overall[sorted_idx]
  names(imp_sorted) <- rownames(importance_df)[sorted_idx]
  feature_names <- rownames(importance_df)[sorted_idx]

  knee_idx <- switch(
    elbow_method,
    uik = {
      x <- seq_along(imp_sorted)
      estimate <- inflection::uik(x, imp_sorted)
      if (length(estimate) > 1) estimate <- estimate[1]
      max(1, min(length(imp_sorted), as.integer(round(estimate))))
    }
  )

  selected_features <- feature_names[seq_len(knee_idx)]

  trimmed_features <- setdiff(selected_features, exclude_variables %||% character())
  missing_feats <- setdiff(trimmed_features, names(omic_df))
  if (length(missing_feats) > 0) {
    stop(
      sprintf(
        "Selected features not found in omic dataset: %s",
        paste(missing_feats, collapse = ", ")
      )
    )
  }

  trimmed_columns <- unique(c(id_column, target_column, trimmed_features))
  trimmed_df <- omic_df[, trimmed_columns, drop = FALSE]

  final_output_path <- output_path
  if (is.null(final_output_path) && !is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    final_output_path <- file.path(output_dir, "omic_selected_features.rds")
  }

  if (!is.null(final_output_path)) {
    dir.create(dirname(final_output_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(trimmed_df, final_output_path)
    # Verbose output omitted here because callers already report the saved data.
  }

  list(
    importance = importance_df,
    sorted_importance = imp_sorted,
    selected_features = trimmed_features,
    knee_index = knee_idx,
    trimmed_data = trimmed_df,
    output_path = final_output_path,
    importance_values = importance_df,
    sorted_importance = imp_sorted
  )
}

.mlasdo_as_model_object <- function(model) {
  if (inherits(model, "train")) {
    return(model)
  }
  if (is.character(model) && length(model) == 1) {
    path <- normalizePath(model, winslash = "/", mustWork = TRUE)
    obj <- readRDS(path)
    if (!inherits(obj, "train")) {
      stop(sprintf("Object loaded from '%s' is not a caret::train model.", path))
    }
    return(obj)
  }
  stop("`model` must be either a caret::train object or a path to an .rds file.")
}

.mlasdo_load_existing_feature_selection <- function(dir,
                                                    id_column,
                                                    target_column) {
  dir <- normalizePath(dir, winslash = "/", mustWork = TRUE)
  trimmed_path <- file.path(dir, "omic_selected_features.rds")
  if (!file.exists(trimmed_path)) {
    stop(sprintf("Existing feature-selection directory '%s' is missing 'omic_selected_features.rds'.", dir))
  }

  importance_path <- file.path(dir, "feature_importance_values.rds")
  summary_path <- file.path(dir, "feature_importance_summary.rds")
  if (!file.exists(importance_path) || !file.exists(summary_path)) {
    assets_dir <- file.path(dirname(dir), "report_assets")
    importance_candidate <- file.path(assets_dir, "feature_importance_values.rds")
    summary_candidate <- file.path(assets_dir, "feature_importance_summary.rds")
    if (file.exists(importance_candidate) && file.exists(summary_candidate)) {
      importance_path <- importance_candidate
      summary_path <- summary_candidate
    } else {
      stop(
        sprintf(
          "Existing feature-selection directory '%s' must contain feature importance files.",
          dir
        )
      )
    }
  }
  trimmed <- readRDS(trimmed_path)
  if (!all(c(id_column, target_column) %in% names(trimmed))) {
    stop("Trimmed dataset does not contain the required ID/target columns.")
  }
  importance_df <- readRDS(importance_path)
  if (!all(c("Feature", "Importance") %in% names(importance_df))) {
    stop("feature_importance_values.rds must contain 'Feature' and 'Importance' columns.")
  }
  summary_info <- readRDS(summary_path)
  knee_index <- summary_info$knee_index
  selected_features <- summary_info$selected_features
  if (is.null(knee_index) || !is.numeric(knee_index)) {
    knee_index <- nrow(importance_df)
  }
  if (is.null(selected_features)) {
    selected_features <- importance_df$Feature[seq_len(min(knee_index, nrow(importance_df)))]
  }
  sorted_importance <- importance_df$Importance
  names(sorted_importance) <- importance_df$Feature
  list(
    importance = importance_df,
    sorted_importance = sorted_importance,
    selected_features = selected_features,
    knee_index = knee_index,
    trimmed_data = trimmed,
    output_path = trimmed_path,
    importance_values = importance_df
  )
}

.mlasdo_build_rf_report_payload <- function(rf_result) {
  if (is.null(rf_result) || is.null(rf_result$model)) {
    return(NULL)
  }
  results_df <- rf_result$model$results
  if (is.null(results_df) || nrow(results_df) == 0) {
    return(NULL)
  }
  list(
    results = results_df,
    best_tune = rf_result$model$bestTune,
    metric = rf_result$model$metric %||% "Balanced_Accuracy",
    importance = rf_result$metadata$importance %||% NA_character_
  )
}

.mlasdo_build_svm_report_payload <- function(svm_results) {
  if (is.null(svm_results)) {
    return(NULL)
  }
  kernel_details <- lapply(svm_results$kernel_results, function(res) {
    params <- res$best_hyperparameters
    if (is.data.frame(params)) {
      params <- as.list(params[1, , drop = TRUE])
    }
    list(
      kernel = res$kernel,
      mean_balanced_accuracy = res$mean_balanced_accuracy,
      best_hyperparameters = params,
      metrics_path = normalizePath(res$artifacts$metrics, winslash = "/", mustWork = TRUE),
      distances_path = normalizePath(res$artifacts$distances, winslash = "/", mustWork = TRUE)
    )
  })
  list(
    summary = svm_results$summary,
    best_kernel = svm_results$best_kernel,
    num_folds = if (!is.null(svm_results$shared_folds)) length(svm_results$shared_folds) else NA_integer_,
    kernel_details = kernel_details,
    class_roles = svm_results$class_roles
  )
}

.mlasdo_build_anomaly_payload <- function(anomalies) {
  if (is.null(anomalies)) {
    return(NULL)
  }
  method_label <- if (identical(anomalies$method, "sd")) {
    sprintf("Mean ± %.2f SD", anomalies$sd_multiplier)
  } else {
    sprintf("Z-score threshold |z| > %.2f", anomalies$z_threshold)
  }
  list(
    method = anomalies$method,
    method_label = method_label,
    sd_multiplier = anomalies$sd_multiplier,
    z_threshold = anomalies$z_threshold,
    thresholds = anomalies$thresholds,
    counts = anomalies$counts,
    anomalies = anomalies$anomalies,
    plot_path = anomalies$plot_path,
    positive_class = anomalies$positive_class,
    negative_class = anomalies$negative_class,
    clinical_tests = anomalies$clinical_tests
  )
}

.mlasdo_render_report <- function(report_path) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' must be installed to render the importance report.")
  }
  report_dir <- dirname(report_path)
  report_file <- basename(report_path)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(report_dir)
  output_file <- rmarkdown::render(
    input = report_file,
    output_format = "html_document",
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )
  if (!grepl("^(/|[A-Za-z]:)", output_file)) {
    output_file <- file.path(report_dir, output_file)
  }
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

#' Run the MLASDO pipeline
#' 
#' @description
#' Loads the omic and clinical datasets, validates the ID/target columns, and
#' returns aligned data frames so the remaining analytical scripts can rely on
#' a single entry point. Once aligned, a feature-selection step based on a
#' random forest is executed to mirror the current MLASDO workflow and the
#' trimmed omic dataset with the most important features is stored. All
#' artifacts are organised in dedicated subfolders under the provided
#' `output_dir`.
#'
#' @param omic_data Path to an `.rds` file or a data frame with molecular
#'   measurements. Defaults to `data/datasets/omic_data_ppmi_sc_npc_nm.rds`
#'   inside the repository.
#' @param clinical_data Path or data frame with the clinical annotations. Defaults
#'   to `data/datasets/clinic_data_PPMI_SC_no_mutation.rds`.
#' @param id_column Name of the column that uniquely identifies each subject.
#' @param target_column Name of the response variable used across downstream
#'   analyses.
#' @param target_positive_class Label within `target_column` treated as the
#'   positive class (defaults to `"Case"`).
#' @param target_negative_class Label treated as the negative class (defaults to
#'   `"Control"`). Both values must exist in the dataset.
#' @param output_dir Directory where MLASDO will create subfolders for
#'   preprocessing and feature selection. The directory is created if it does
#'   not exist.
#' @param always_split_variables Character vector passed to
#'   `ranger::ranger(always.split.variables = ...)` to force covariates to be
#'   part of each split. Values are pulled from the clinical dataset (and added
#'   to the modeling table if needed). Defaults to `c("age_at_baseline", "sex")`.
#' @param feature_selection_exec Execute the RF-based feature-selection step.
#'   When `FALSE`, `feature_selection_existing_dir` must point to a folder with
#'   previously generated artifacts (trimmed dataset and importance files).
#' @param feature_selection_existing_dir Directory containing an existing
#'   feature-selection run (expects `omic_selected_features.rds`,
#'   `feature_importance_values.rds`, and `feature_importance_summary.rds`).
#' @param rf_mtry Candidate `mtry` values. If `NULL`, defaults to the sequence
#'   originally used in MLASDO (`floor(seq(2, sqrt(p), length.out = 3))`).
#' @param rf_splitrule Candidate split rules (defaults to
#'   `c("gini", "extratrees", "hellinger")`).
#' @param rf_min_node_size Candidate `min.node.size` values (defaults to
#'   `c(1, 5, 10)`).
#' @param rf_train_method Resampling strategy for `caret::trainControl`
#'   (`"repeatedcv"` or `"cv"`). Defaults to `"repeatedcv"`.
#' @param rf_train_number Number of folds used in resampling (defaults to 10).
#' @param rf_train_repeats Number of repeats when `rf_train_method = "repeatedcv"`
#'   (defaults to 5).
#' @param rf_num_trees Number of trees passed to `ranger` (defaults to 500).
#' @param rf_model_filename File name used to persist the feature-selection
#'   model inside the feature-selection folder. Defaults to
#'   `"rf_feature_selection.rds"`.
#' @param rf_importance Importance measure passed to `ranger`. Defaults to
#'   `"permutation"` but it can be set to any value supported by `ranger`
#'   (e.g., `"impurity"`).
#' @param n_cores Number of CPU cores allocated to MLASDO (feature selection and
#'   SVM training). Defaults to `parallel::detectCores() - 1`.
#' @param seed Random seed applied before starting the feature-selection phase
#'   and reused for the SVM folds.
#' @param svm_exec Run the post-selection SVM training workflow once the
#'   trimmed dataset is available. Defaults to `TRUE`.
#' @param svm_output_dir Directory used to store the SVM artifacts. Defaults to
#'   a `svm_models` subfolder inside `output_dir`.
#' @param svm_existing_dir Directory with previously trained SVM models. Used
#'   when `svm_exec = FALSE` but the artifacts should be reused.
#' @param svm_kernels Character vector with any subset of `"linear"`,
#'   `"radial"` or `"polynomial"` kernels passed to
#'   `mlasdo_train_svm_models()`.
#' @param svm_linear_cost_values Candidate `C` values for the linear SVM grid
#'   when `svm_kernels` includes `"linear"`.
#' @param svm_radial_cost_values Candidate `C` values for the radial SVM grid.
#' @param svm_radial_sigma_values Candidate `sigma`/`gamma` values for the
#'   radial kernel.
#' @param svm_polynomial_cost_values Candidate `C` values for the polynomial
#'   kernel.
#' @param svm_polynomial_degree_values Candidate degrees evaluated for the
#'   polynomial kernel.
#' @param svm_polynomial_scale_values Candidate `scale` (`gamma`) values
#'   evaluated for the polynomial kernel.
#' @param svm_num_folds Number of folds reused across the SVM kernels. Defaults
#'   to 10 and independent from the RF resampling settings.
#' @param svm_sampling Sampling strategy passed to the underlying
#'   `trainControl(sampling = ...)`.
#' @param anomalies_exec Run the anomalous-individuals selection after the SVM
#'   stage finishes (defaults to `TRUE`).
#' @param anomalies_method Strategy used to flag anomalies. Defaults to `"sd"`
#'   (mean ± SD) but `"zscore"` can also be used.
#' @param anomalies_sd_multiplier Multiplier applied to the per-class SD when
#'   `anomalies_method = "sd"`.
#' @param anomalies_z_threshold Absolute z-score threshold applied when
#'   `anomalies_method = "zscore"`.
#' @param anomalies_output_dir Directory where anomaly-related artifacts will be
#'   stored. Defaults to an `anomalies` subfolder under `output_dir`.
#' @param anomaly_numeric_vars Optional character vector with the clinical
#'   numeric variables evaluated when comparing anomalous individuals.
#' @param anomaly_categorical_vars Optional character vector with the clinical
#'   categorical variables evaluated during the anomaly analysis.
#' 
#' @return An invisible list with `omic`, `clinical`, `feature_selection_model`,
#'   the trimmed omic dataset, and `metadata` entries.
#' @export
run_mlasdo <- function(
    omic_data = "data/datasets/omic_data_ppmi_sc_npc_nm.rds",
    clinical_data = "data/datasets/clinic_data_PPMI_SC_no_mutation.rds",
    id_column = "participant_id",
    target_column = "diagnosis_type",
    target_positive_class = "Case",
    target_negative_class = "Control",
    output_dir = NULL,
    always_split_variables = c("age_at_baseline", "sex"),
    feature_selection_exec = TRUE,
    feature_selection_existing_dir = NULL,
    rf_mtry = NULL,
    rf_splitrule = c("gini", "extratrees", "hellinger"),
    rf_min_node_size = c(1, 5, 10),
    rf_train_method = "repeatedcv",
    rf_train_number = 10,
    rf_train_repeats = 5,
    rf_num_trees = 500,
    rf_model_filename = "rf_feature_selection.rds",
    rf_importance = "permutation",
    n_cores = max(1, parallel::detectCores() - 1),
    seed = 123,
    svm_exec = TRUE,
    svm_output_dir = NULL,
    svm_existing_dir = NULL,
    svm_kernels = c("linear", "radial", "polynomial"),
    svm_linear_cost_values = 10^seq(-6, 6, length.out = 13),
    svm_radial_cost_values = 10^seq(-4, 4, length.out = 9),
    svm_radial_sigma_values = 10^seq(-6, 1, length.out = 8),
    svm_polynomial_cost_values = 10^seq(-4, 4, length.out = 9),
    svm_polynomial_degree_values = 2:4,
    svm_polynomial_scale_values = 10^seq(-6, 0, length.out = 7),
    svm_num_folds = 10,
    svm_sampling = "down",
    anomalies_exec = TRUE,
    anomalies_method = "sd",
    anomalies_sd_multiplier = 2,
    anomalies_z_threshold = 1.65,
    anomalies_output_dir = NULL,
    anomaly_numeric_vars = NULL,
    anomaly_categorical_vars = NULL
  ) {
  rf_train_method <- match.arg(rf_train_method, c("repeatedcv", "cv"))
  rf_importance <- match.arg(
    rf_importance,
    c("permutation", "impurity", "impurity_corrected", "none")
  )
  anomalies_method <- match.arg(anomalies_method, c("sd", "zscore"))
  anomaly_sd_label <- if (identical(anomalies_method, "sd")) {
    sprintf("mean ± %.2f SD", anomalies_sd_multiplier)
  } else {
    sprintf("|z| > %.2f", anomalies_z_threshold)
  }

  if (is.null(output_dir)) {
    stop("output_dir must be provided so the feature-selection model can be saved there.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  preprocessing_dir <- .mlasdo_prepare_subdir(output_dir, "preprocessing")
  feature_dir <- NULL
  if (isTRUE(feature_selection_exec)) {
    feature_dir <- .mlasdo_prepare_subdir(output_dir, "feature_selection")
  } else {
    if (is.null(feature_selection_existing_dir)) {
      stop("feature_selection_existing_dir must be provided when feature_selection_exec = FALSE.")
    }
    feature_dir <- normalizePath(feature_selection_existing_dir, winslash = "/", mustWork = TRUE)
  }

  domain_expert_dir <- .mlasdo_prepare_subdir(output_dir, "domain_expert_knowledge")

  rf_model_path <- NULL
  if (isTRUE(feature_selection_exec)) {
    rf_model_path <- file.path(feature_dir, rf_model_filename)
    rf_model_path <- normalizePath(rf_model_path, winslash = "/", mustWork = FALSE)
  }

  .mlasdo_log_step("Starting data preparation.")
  prep <- prepare_mlasdo_data(
    omic_data = omic_data,
    clinical_data = clinical_data,
    id_column = id_column,
    target_column = target_column,
    verbose = TRUE
  )
  .mlasdo_log_step("Data preparation completed.")
  omic_aligned <- prep$omic
  clinical_aligned <- prep$clinical
  
  target_cfg <- .mlasdo_configure_target(
    omic_aligned[[target_column]],
    positive_class = target_positive_class,
    negative_class = target_negative_class
  )
  resolved_positive <- target_cfg$positive
  resolved_negative <- target_cfg$negative
  if (is.null(resolved_positive) || is.null(resolved_negative)) {
    lvls <- levels(target_cfg$values)
    if (is.null(resolved_negative) && length(lvls) >= 1) {
      resolved_negative <- lvls[1]
    }
    if (is.null(resolved_positive) && length(lvls) >= 2) {
      resolved_positive <- lvls[2]
    }
  }
  omic_aligned[[target_column]] <- target_cfg$values
  clinical_aligned[[target_column]] <- factor(
    clinical_aligned[[target_column]],
    levels = levels(target_cfg$values)
  )
  .mlasdo_save_preprocessing_data(
    preprocess_dir = preprocessing_dir,
    omic = omic_aligned,
    clinical = clinical_aligned,
    verbose = FALSE
  )
  .mlasdo_log_step("Aligned datasets stored.")

  .mlasdo_log_step("Building Domain Expert Knowledge table.")
  domain_expert_result <- .mlasdo_compute_domain_expert_table(
    clinical_data = clinical_aligned,
    target_column = target_column,
    positive_class = resolved_positive,
    negative_class = resolved_negative,
    output_dir = domain_expert_dir,
    verbose = TRUE
  )
  if (is.null(domain_expert_result)) {
    .mlasdo_log_step("Domain Expert Knowledge table not generated (no valid numeric covariates).")
  } else {
    .mlasdo_log_step("Domain Expert Knowledge table generated.")
  }

  rf_result <- NULL
  top_features <- NULL
  if (isTRUE(feature_selection_exec)) {
    .mlasdo_log_step("Starting Random Forest feature selection.")
    rf_result <- .mlasdo_run_rf_feature_selection(
      data = omic_aligned,
      n_cores = n_cores,
      clinical_data = clinical_aligned,
      id_column = id_column,
      target_column = target_column,
      always_split_variables = always_split_variables,
      rf_mtry = rf_mtry,
      rf_splitrule = rf_splitrule,
      rf_min_node_size = rf_min_node_size,
      train_method = rf_train_method,
      train_number = rf_train_number,
      train_repeats = rf_train_repeats,
      rf_num_trees = rf_num_trees,
      rf_model_path = rf_model_path,
      seed = seed,
      importance = rf_importance,
      verbose = TRUE
    )
    .mlasdo_log_step("Feature selection completed.")

    top_features <- mlasdo_select_top_features(
      model = rf_result$model,
      omic_data = omic_aligned,
      id_column = id_column,
      target_column = target_column,
      exclude_variables = always_split_variables,
      output_dir = feature_dir,
      verbose = TRUE
    )
    .mlasdo_log_step("Selected features stored.")
  } else {
    .mlasdo_log_step("Reusing existing feature selection from '%s'.", feature_dir)
    top_features <- .mlasdo_load_existing_feature_selection(
      dir = feature_dir,
      id_column = id_column,
      target_column = target_column
    )
  }
  
  svm_results <- NULL
  if (isTRUE(svm_exec)) {
    .mlasdo_log_step(
      "Training SVM models (kernels: %s).",
      paste(svm_kernels, collapse = ", ")
    )
    svm_dir <- svm_output_dir
    if (is.null(svm_dir)) {
      svm_dir <- file.path(output_dir, "svm_models")
    }
    svm_results <- mlasdo_train_svm_models(
      omic_data = top_features$trimmed_data,
      id_column = id_column,
      target_column = target_column,
      output_dir = svm_dir,
      kernels = svm_kernels,
      linear_cost_values = svm_linear_cost_values,
      radial_cost_values = svm_radial_cost_values,
      radial_sigma_values = svm_radial_sigma_values,
      polynomial_cost_values = svm_polynomial_cost_values,
      polynomial_degree_values = svm_polynomial_degree_values,
      polynomial_scale_values = svm_polynomial_scale_values,
      num_folds = svm_num_folds,
      folds_seed = seed,
      n_cores = n_cores,
      sampling = svm_sampling,
      positive_class = resolved_positive,
      negative_class = resolved_negative
    )
    .mlasdo_log_step("SVM training completed. Best kernel: %s.", svm_results$best_kernel %||% "not available")
  } else if (!is.null(svm_existing_dir)) {
    .mlasdo_log_step("Reusing existing SVM results from '%s'.", svm_existing_dir)
    svm_results <- .mlasdo_load_existing_svm_results(svm_existing_dir)
  } else {
    .mlasdo_log_step("SVM training skipped.")
  }
  if (!is.null(svm_results) && !is.null(svm_results$class_roles$positive)) {
    resolved_positive <- svm_results$class_roles$positive
  }
  if (!is.null(svm_results) && !is.null(svm_results$class_roles$negative)) {
    resolved_negative <- svm_results$class_roles$negative
  }
  
  anomalies_result <- NULL
  individual_interpretability <- NULL
  anomalies_dir <- NULL
  if (!is.null(svm_results) && isTRUE(anomalies_exec)) {
    if (is.null(anomalies_output_dir)) {
      anomalies_dir <- .mlasdo_prepare_subdir(output_dir, "anomalies")
    } else {
      anomalies_dir <- .mlasdo_prepare_directory(anomalies_output_dir)
    }
    best_kernel_name <- svm_results$best_kernel
    if (!is.null(best_kernel_name) &&
        best_kernel_name %in% names(svm_results$kernel_results)) {
      .mlasdo_log_step("Running anomaly detection (%s).", anomalies_method)
      distances_path <- svm_results$kernel_results[[best_kernel_name]]$artifacts$distances
      anomalies_result <- .mlasdo_detect_anomalies(
        distances = distances_path,
        method = anomalies_method,
        sd_multiplier = anomalies_sd_multiplier,
        z_threshold = anomalies_z_threshold,
        output_dir = anomalies_dir,
        positive_class = resolved_positive,
        negative_class = resolved_negative
      )
      if (!is.null(anomalies_result) && !is.null(anomalies_result$anomalies)) {
        .mlasdo_log_step(
          "Anomaly detection finished (n = %d).",
          nrow(anomalies_result$anomalies)
        )
      }
    } else {
      warning("Unable to locate best-kernel distances for anomaly detection.")
    }
  } else if (isTRUE(anomalies_exec)) {
    .mlasdo_log_step("Skipping anomaly detection: SVM results unavailable.")
  } else {
    .mlasdo_log_step("Anomaly detection disabled.")
  }
  if (!is.null(anomalies_result) && !is.null(anomalies_dir)) {
    .mlasdo_log_step("Running statistical tests for anomalous individuals.")
    clinical_tests <- .mlasdo_anomaly_clinical_analysis(
      anomalies = anomalies_result,
      clinical_data = clinical_aligned,
      positive_class = resolved_positive,
      negative_class = resolved_negative,
      numeric_vars = anomaly_numeric_vars,
      categorical_vars = anomaly_categorical_vars,
      sd_label = anomaly_sd_label,
      output_dir = anomalies_dir,
      id_column = id_column,
      target_column = target_column
    )
    if (!is.null(clinical_tests)) {
      anomalies_result$clinical_tests <- clinical_tests
      .mlasdo_log_step("Statistical tests generated.")
    } else {
      .mlasdo_log_step("No statistical tests generated (no valid covariates).")
    }
  }
  
  anomaly_payload <- NULL
  if (!is.null(anomalies_result)) {
    anomaly_payload <- .mlasdo_build_anomaly_payload(anomalies_result)
  }

  if (!is.null(anomalies_result) && !is.null(domain_expert_result)) {
    .mlasdo_log_step("Building individual explanations for anomalous cases.")
    individual_interpretability <- .mlasdo_build_individual_interpretability(
      anomalies = anomalies_result,
      clinical_data = clinical_aligned,
      domain_expert = domain_expert_result,
      id_column = id_column,
      target_column = target_column
    )
    if (!is.null(individual_interpretability)) {
      .mlasdo_log_step("Individual explanations ready.")
    } else {
      .mlasdo_log_step("Unable to build individual explanations (insufficient information).")
    }
  } else {
    .mlasdo_log_step("Individual explanations skipped (missing anomalies or expert table).")
  }
  
  rf_payload <- .mlasdo_build_rf_report_payload(rf_result)
  .mlasdo_log_step("Generating final report.")
  report_files <- .mlasdo_write_importance_report(
    report_dir = output_dir,
    sorted_importance = top_features$sorted_importance,
    knee_index = top_features$knee_index,
    selected_features = top_features$selected_features,
    rf_performance = rf_payload,
    svm_results = svm_results,
    anomalies = anomaly_payload,
    domain_expert = domain_expert_result,
    individual_interpretability = individual_interpretability
  )
  
  rendered_report <- .mlasdo_render_report(report_files$report_path)
  .mlasdo_log_step("Final report rendered: %s.", rendered_report %||% "not available")

  run_metadata <- list(
    id_column = id_column,
    target_column = target_column,
    classes = list(
      positive = resolved_positive,
      negative = resolved_negative
    ),
    timestamp = Sys.time(),
    directories = list(
      root = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
      preprocessing = preprocessing_dir,
      feature_selection = feature_dir,
      domain_expert = domain_expert_dir
    ),
    sources = prep$sources,
      feature_selection = list(
        model = if (!is.null(rf_result)) rf_result$metadata else list(source = "existing", directory = feature_dir),
        selected_features = list(
          count = length(top_features$selected_features),
          knee_index = top_features$knee_index,
          dataset = top_features$output_path,
          report_rmd = report_files$report_path,
          report_html = rendered_report
        )
    ),
    rf = list(
      importance = if (!is.null(rf_result)) rf_result$metadata$importance else NULL,
      metric = if (!is.null(rf_result) && !is.null(rf_result$model$metric)) rf_result$model$metric else NULL,
      best_tune = if (!is.null(rf_result)) rf_result$model$bestTune else NULL
    ),
    svm = list(
      enabled = !is.null(svm_results),
      source = if (isTRUE(svm_exec)) "train" else if (!is.null(svm_existing_dir)) "existing" else "disabled",
      output_dir = if (!is.null(svm_results)) svm_results$output_dir else NULL,
      best_kernel = if (!is.null(svm_results)) svm_results$best_kernel else NULL,
      summary = if (!is.null(svm_results)) svm_results$summary else NULL,
      num_folds = svm_num_folds,
      folds_seed = seed,
      grid = list(
        linear_cost = svm_linear_cost_values,
        radial_cost = svm_radial_cost_values,
        radial_sigma = svm_radial_sigma_values,
        polynomial_cost = svm_polynomial_cost_values,
        polynomial_degree = svm_polynomial_degree_values,
        polynomial_scale = svm_polynomial_scale_values
      )
    ),
    anomalies = list(
      enabled = !is.null(anomalies_result),
      method = anomalies_method,
      sd_multiplier = anomalies_sd_multiplier,
      z_threshold = anomalies_z_threshold,
      output_dir = anomalies_dir,
      anomalies_path = if (!is.null(anomalies_result)) anomalies_result$anomalies_path else NULL,
      thresholds_path = if (!is.null(anomalies_result)) anomalies_result$thresholds_path else NULL,
      counts_path = if (!is.null(anomalies_result)) anomalies_result$counts_path else NULL,
      clinical_tests = if (!is.null(anomalies_result)) anomalies_result$clinical_tests else NULL
    ),
    interpretability = list(
      domain_expert = if (!is.null(domain_expert_result)) {
        list(
          table_path = domain_expert_result$table_path,
          csv_path = domain_expert_result$csv_path,
          alpha = domain_expert_result$alpha,
          positive_class = domain_expert_result$positive_class,
          negative_class = domain_expert_result$negative_class
        )
      } else {
        list(
          table_path = NULL,
          csv_path = NULL,
          alpha = NULL,
          positive_class = resolved_positive,
          negative_class = resolved_negative
        )
      },
      individual = if (!is.null(individual_interpretability) && !is.null(individual_interpretability$data)) {
        list(
          count = nrow(individual_interpretability$data)
        )
      } else {
        NULL
      }
    )
  )

  metadata_path <- .mlasdo_save_metadata(output_dir, run_metadata)

  invisible(
    list(
      preprocessing_dir = preprocessing_dir,
      feature_selection_dir = feature_dir,
      omic = omic_aligned,
      clinical = clinical_aligned,
      feature_selection_model = if (!is.null(rf_result)) rf_result$model else NULL,
      omic_selected = top_features$trimmed_data,
      svm_results = svm_results,
      anomalies = anomalies_result,
      domain_expert = domain_expert_result,
      individual_interpretability = individual_interpretability,
      class_roles = list(
        positive = resolved_positive,
        negative = resolved_negative
      ),
      metadata_path = metadata_path,
      metadata = run_metadata
    )
  )
}

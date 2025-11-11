.mlasdo_compute_domain_expert_table <- function(clinical_data,
                                                target_column,
                                                positive_class,
                                                negative_class,
                                                output_dir,
                                                alpha = 0.05,
                                                verbose = interactive()) {
  if (is.null(clinical_data) || nrow(clinical_data) == 0) {
    if (isTRUE(verbose)) {
      message("Clinical dataset is empty; skipping domain expert table generation.")
    }
    return(NULL)
  }

  .mlasdo_require_columns(clinical_data, target_column, "clinical dataset")
  target_values <- clinical_data[[target_column]]
  if (!is.factor(target_values)) {
    target_values <- factor(target_values)
  }

  if (is.null(positive_class) || is.null(negative_class)) {
    lvls <- levels(target_values)
    if (length(lvls) < 2) {
      if (isTRUE(verbose)) {
        message("Target column does not contain at least two classes; domain expert table not generated.")
      }
      return(NULL)
    }
    if (is.null(positive_class)) positive_class <- lvls[1]
    if (is.null(negative_class)) {
      negative_class <- setdiff(lvls, positive_class)[1]
    }
  }

  positive_idx <- !is.na(target_values) & target_values == positive_class
  negative_idx <- !is.na(target_values) & target_values == negative_class
  if (sum(positive_idx) == 0 || sum(negative_idx) == 0) {
    if (isTRUE(verbose)) {
      message("Both positive and negative classes require at least one observation; domain expert table not generated.")
    }
    return(NULL)
  }

  numeric_cols <- names(clinical_data)[vapply(
    clinical_data,
    function(col) is.numeric(col) && length(stats::na.omit(unique(col))) > 1,
    logical(1)
  )]
  numeric_cols <- setdiff(numeric_cols, target_column)
  if (length(numeric_cols) == 0) {
    if (isTRUE(verbose)) {
      message("No numeric clinical covariates available for the domain expert table.")
    }
    return(NULL)
  }

  summary_rows <- lapply(numeric_cols, function(var) {
    pos_vals <- clinical_data[[var]][positive_idx]
    neg_vals <- clinical_data[[var]][negative_idx]
    pos_vals <- pos_vals[!is.na(pos_vals)]
    neg_vals <- neg_vals[!is.na(neg_vals)]
    if (length(pos_vals) < 3 || length(neg_vals) < 3) {
      return(NULL)
    }
    test_res <- tryCatch(
      stats::wilcox.test(pos_vals, neg_vals, exact = FALSE),
      error = function(e) NULL
    )
    if (is.null(test_res)) {
      return(NULL)
    }
    median_pos <- stats::median(pos_vals, na.rm = TRUE)
    median_neg <- stats::median(neg_vals, na.rm = TRUE)
    threshold <- mean(c(median_pos, median_neg))
    direction <- if (isTRUE(median_pos > median_neg)) positive_class else negative_class
    data.frame(
      variable = var,
      positive_median = median_pos,
      negative_median = median_neg,
      threshold = threshold,
      p_value = test_res$p.value,
      direction = direction,
      positive_n = length(pos_vals),
    negative_n = length(neg_vals),
    stringsAsFactors = FALSE
  )
})
summary_rows <- Filter(Negate(is.null), summary_rows)
  if (length(summary_rows) == 0) {
    if (isTRUE(verbose)) {
      message("Not enough numeric variables passed the requirements for the domain expert table.")
    }
    return(NULL)
  }

  table_df <- do.call(rbind, summary_rows)
  table_df <- table_df[order(table_df$p_value), , drop = FALSE]
  table_df$significant <- !is.na(table_df$p_value) & table_df$p_value < alpha
  table_df$p_value <- sprintf("%.3e", table_df$p_value)
  significant_df <- table_df[table_df$significant, , drop = FALSE]

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  table_path <- file.path(output_dir, "domain_expert_table.rds")
  csv_path <- file.path(output_dir, "domain_expert_table.csv")
  to_save <- significant_df
  if (nrow(to_save) == 0) {
    to_save <- table_df[0, , drop = FALSE]
  }
  saveRDS(to_save, table_path)
  utils::write.csv(to_save, csv_path, row.names = FALSE)

  list(
    table = to_save,
    table_path = normalizePath(table_path, winslash = "/", mustWork = FALSE),
    csv_path = normalizePath(csv_path, winslash = "/", mustWork = FALSE),
    positive_class = positive_class,
    negative_class = negative_class,
    alpha = alpha,
    total_candidates = length(numeric_cols)
  )
}

.mlasdo_value_supports_class_change <- function(value,
                                                sample_class,
                                                positive_class,
                                                negative_class,
                                                positive_median,
                                                negative_median,
                                                threshold) {
  if (is.na(value) || is.na(threshold) || is.na(sample_class)) {
    return(NA)
  }
  target_class <- NULL
  if (!is.null(positive_class) && !is.null(negative_class)) {
    if (identical(sample_class, positive_class)) {
      target_class <- negative_class
    } else if (identical(sample_class, negative_class)) {
      target_class <- positive_class
    }
  }
  .mlasdo_value_matches_class(
    value = value,
    class_label = target_class,
    positive_class = positive_class,
    negative_class = negative_class,
    positive_median = positive_median,
    negative_median = negative_median,
    threshold = threshold
  )
}

.mlasdo_value_matches_class <- function(value,
                                        class_label,
                                        positive_class,
                                        negative_class,
                                        positive_median,
                                        negative_median,
                                        threshold) {
  if (is.null(class_label) || is.na(class_label) || is.na(value) || is.na(threshold)) {
    return(NA)
  }
  if (identical(class_label, positive_class)) {
    if (is.na(positive_median) || is.na(negative_median)) return(NA)
    if (positive_median >= negative_median) {
      return(value >= threshold)
    } else {
      return(value <= threshold)
    }
  }
  if (identical(class_label, negative_class)) {
    if (is.na(positive_median) || is.na(negative_median)) return(NA)
    if (negative_median >= positive_median) {
      return(value >= threshold)
    } else {
      return(value <= threshold)
    }
  }
  NA
}

.mlasdo_build_individual_interpretability <- function(anomalies,
                                                      clinical_data,
                                                      domain_expert,
                                                      id_column,
                                                      target_column) {
  if (is.null(anomalies) || is.null(anomalies$anomalies) ||
      nrow(anomalies$anomalies) == 0 || is.null(domain_expert) ||
      is.null(domain_expert$table) || nrow(domain_expert$table) == 0 ||
      is.null(clinical_data) || nrow(clinical_data) == 0) {
    return(NULL)
  }

  domain_tbl <- domain_expert$table
  domain_tbl$variable <- as.character(domain_tbl$variable)
  domain_vars <- unique(stats::na.omit(domain_tbl$variable))
  domain_vars <- intersect(domain_vars, names(clinical_data))
  if (length(domain_vars) == 0) {
    return(NULL)
  }

  anomaly_df <- anomalies$anomalies
  matched_idx <- match(anomaly_df$id, clinical_data[[id_column]])
  valid <- !is.na(matched_idx)
  if (!any(valid)) {
    return(NULL)
  }
  anomaly_df <- anomaly_df[valid, , drop = FALSE]
  matched_idx <- matched_idx[valid]
  clinical_subset <- clinical_data[matched_idx, , drop = FALSE]

  build_label <- function(var_name, value, threshold) {
    value_fmt <- if (is.numeric(value)) {
      format(value, digits = 4, trim = TRUE)
    } else {
      as.character(value)
    }
    threshold_fmt <- if (is.numeric(threshold)) {
      format(threshold, digits = 4, trim = TRUE)
    } else {
      as.character(threshold)
    }
    sprintf(
      "<li><strong>%s</strong>: %s (threshold: %s)</li>",
      var_name,
      value_fmt,
      threshold_fmt
    )
  }

  rows <- lapply(seq_len(nrow(anomaly_df)), function(i) {
    sample_id <- anomaly_df$id[i]
    sample_class <- as.character(anomaly_df$diagnosis_true[i])
    sample_pred <- anomaly_df$diagnosis_pred[i]
    sample_row <- clinical_subset[i, , drop = FALSE]
    support_list <- character()
    oppose_list <- character()
    var_available <- FALSE

    for (j in seq_len(nrow(domain_tbl))) {
      var_name <- domain_tbl$variable[j]
      if (!var_name %in% domain_vars) next
      value <- sample_row[[var_name]]
      threshold <- domain_tbl$threshold[j]
      if (is.null(value) || all(is.na(value)) || is.na(threshold)) {
        next
      }
      var_available <- TRUE
      support_flag <- .mlasdo_value_supports_class_change(
        value = value,
        sample_class = sample_class,
        positive_class = domain_expert$positive_class,
        negative_class = domain_expert$negative_class,
        positive_median = domain_tbl$positive_median[j],
        negative_median = domain_tbl$negative_median[j],
        threshold = threshold
      )
      if (is.na(support_flag)) {
        next
      }
      label <- build_label(var_name, value, threshold)
      if (isTRUE(support_flag)) {
        support_list <- c(support_list, label)
      } else {
        oppose_list <- c(oppose_list, label)
      }
    }

    if (!var_available) {
      return(NULL)
    }

    support_count <- length(support_list)
    oppose_count <- length(oppose_list)
    total_count <- support_count + oppose_count
    support_ratio <- if (total_count > 0) support_count / total_count else NA_real_

    data.frame(
      participant_id = sample_id,
      diagnosis_true = as.character(sample_class),
      diagnosis_pred = as.character(sample_pred),
      support_count = support_count,
      oppose_count = oppose_count,
      support_ratio = support_ratio,
      support_details = if (length(support_list) > 0) paste(support_list, collapse = "") else "",
      oppose_details = if (length(oppose_list) > 0) paste(oppose_list, collapse = "") else "",
      stringsAsFactors = FALSE
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(NULL)
  }
  data <- do.call(rbind, rows)
  data <- data[order(data$support_ratio, decreasing = TRUE), , drop = FALSE]
  list(
    data = data,
    positive_class = domain_expert$positive_class,
    negative_class = domain_expert$negative_class,
    table_path = domain_expert$table_path
  )
}

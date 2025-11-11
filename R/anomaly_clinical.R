.mlasdo_anomaly_clinical_analysis <- function(anomalies,
                                              clinical_data,
                                              positive_class,
                                              negative_class,
                                              numeric_vars,
                                              categorical_vars,
                                              sd_label,
                                              output_dir,
                                              id_column = "participant_id",
                                              target_column = "diagnosis_type") {
  if (is.null(anomalies) || nrow(anomalies$anomalies) == 0) {
    return(NULL)
  }

  prepared <- .mlasdo_prepare_anomaly_clinical_data(
    clinical_data = clinical_data,
    anomalies = anomalies$anomalies,
    positive_class = positive_class,
    negative_class = negative_class
  )
  if (is.null(prepared)) {
    return(NULL)
  }

  data_grouped <- prepared$data
  group_levels <- prepared$group_levels
  clinical_dir <- file.path(output_dir, "clinical_tests")
  dir.create(clinical_dir, recursive = TRUE, showWarnings = FALSE)

  inferred <- .mlasdo_infer_clinical_variable_types(
    data = data_grouped,
    numeric_vars = numeric_vars,
    categorical_vars = categorical_vars,
    id_column = id_column,
    target_column = target_column
  )
  numeric_vars <- inferred$numeric_vars
  categorical_vars <- inferred$categorical_vars
  if ((is.null(numeric_vars) || length(numeric_vars) == 0) &&
      (is.null(categorical_vars) || length(categorical_vars) == 0)) {
    return(NULL)
  }

  normalized <- .mlasdo_normalize_clinical_variables(
    data_grouped,
    numeric_vars,
    categorical_vars
  )
  data_grouped <- normalized$data
  valid_numeric <- normalized$numeric_vars
  valid_categorical <- normalized$categorical_vars

  numeric_pairs_neg <- list(
    c(prepared$anomalous_negative, prepared$negative_label),
    c(prepared$anomalous_negative, prepared$positive_label)
  )
  numeric_pairs_pos <- list(
    c(prepared$anomalous_positive, prepared$positive_label),
    c(prepared$anomalous_positive, prepared$negative_label)
  )
  categorical_pairs_neg <- list(
    c(prepared$anomalous_negative, prepared$negative_label)
  )
  categorical_pairs_pos <- list(
    c(prepared$anomalous_positive, prepared$positive_label)
  )

  scenario_neg <- .mlasdo_run_anomaly_comparison(
    data = data_grouped,
    numeric_vars = valid_numeric,
    categorical_vars = valid_categorical,
    focus_label = prepared$anomalous_negative,
    numeric_comparisons = numeric_pairs_neg,
    categorical_comparisons = categorical_pairs_neg,
    group_levels = group_levels,
    sd_label = sd_label,
    plot_prefix = "anomalous_negative",
    output_dir = clinical_dir
  )

  scenario_pos <- .mlasdo_run_anomaly_comparison(
    data = data_grouped,
    numeric_vars = valid_numeric,
    categorical_vars = valid_categorical,
    focus_label = prepared$anomalous_positive,
    numeric_comparisons = numeric_pairs_pos,
    categorical_comparisons = categorical_pairs_pos,
    group_levels = group_levels,
    sd_label = sd_label,
    plot_prefix = "anomalous_positive",
    output_dir = clinical_dir
  )

  list(
    negative = scenario_neg,
    positive = scenario_pos,
    sd_label = sd_label
  )
}

.mlasdo_prepare_anomaly_clinical_data <- function(clinical_data,
                                                  anomalies,
                                                  positive_class,
                                                  negative_class) {
  if (is.null(clinical_data) || nrow(clinical_data) == 0) {
    return(NULL)
  }
  if (!all(c("participant_id", "diagnosis_type") %in% names(clinical_data))) {
    stop("Clinical dataset must contain 'participant_id' and 'diagnosis_type' columns.")
  }

  pos_label <- positive_class %||% "Case"
  neg_label <- negative_class %||% "Control"
  anomalous_positive <- paste("Anomalous", pos_label)
  anomalous_negative <- paste("Anomalous", neg_label)
  group_levels <- c(
    pos_label,
    neg_label,
    anomalous_positive,
    anomalous_negative
  )

  anomaly_ids <- unique(anomalies$id)
  anomaly_positive_ids <- unique(anomalies$id[anomalies$diagnosis_true == pos_label])
  anomaly_negative_ids <- unique(anomalies$id[anomalies$diagnosis_true == neg_label])

  data <- clinical_data
  data$group_label <- pos_label
  data$group_label[data$diagnosis_type == neg_label] <- neg_label
  data$group_label[data$participant_id %in% anomaly_positive_ids &
                     data$diagnosis_type == pos_label] <- anomalous_positive
  data$group_label[data$participant_id %in% anomaly_negative_ids &
                     data$diagnosis_type == neg_label] <- anomalous_negative

  data$group_label <- factor(
    data$group_label,
    levels = group_levels
  )

  valid_rows <- !is.na(data$group_label)
  data <- data[valid_rows, , drop = FALSE]
  if (nrow(data) == 0) {
    return(NULL)
  }

  list(
    data = data,
    group_levels = group_levels,
    positive_label = pos_label,
    negative_label = neg_label,
    anomalous_positive = anomalous_positive,
    anomalous_negative = anomalous_negative
  )
}

.mlasdo_normalize_clinical_variables <- function(data,
                                                 numeric_vars,
                                                 categorical_vars) {
  available_vars <- names(data)
  valid_numeric <- intersect(numeric_vars %||% character(0), available_vars)
  valid_categorical <- intersect(categorical_vars %||% character(0), available_vars)

  if (length(valid_numeric) > 0) {
    data[valid_numeric] <- lapply(data[valid_numeric], function(col) {
      suppressWarnings(as.numeric(col))
    })
  }
  if (length(valid_categorical) > 0) {
    data[valid_categorical] <- lapply(data[valid_categorical], function(col) {
      if (is.factor(col)) {
        col
      } else {
        factor(col)
      }
    })
  }

  list(
    data = data,
    numeric_vars = valid_numeric,
    categorical_vars = valid_categorical
  )
}

.mlasdo_infer_clinical_variable_types <- function(data,
                                                  numeric_vars,
                                                  categorical_vars,
                                                  id_column,
                                                  target_column) {
  available <- names(data)
  blacklist <- unique(stats::na.omit(c(id_column, target_column, "group_label")))
  candidate_vars <- setdiff(available, blacklist)

  inferred_numeric <- numeric_vars
  if (is.null(inferred_numeric) || length(inferred_numeric) == 0) {
    inferred_numeric <- candidate_vars[vapply(
      candidate_vars,
      function(col) is.numeric(data[[col]]),
      logical(1)
    )]
  } else {
    inferred_numeric <- intersect(inferred_numeric, available)
  }

  inferred_categorical <- categorical_vars
  if (is.null(inferred_categorical) || length(inferred_categorical) == 0) {
    inferred_categorical <- candidate_vars[vapply(
      candidate_vars,
      function(col) {
        is.factor(data[[col]]) ||
          (is.character(data[[col]]) && length(unique(stats::na.omit(data[[col]]))) > 1)
      },
      logical(1)
    )]
  } else {
    inferred_categorical <- intersect(inferred_categorical, available)
  }

  list(
    numeric_vars = inferred_numeric,
    categorical_vars = inferred_categorical
  )
}

.mlasdo_run_anomaly_comparison <- function(data,
                                           numeric_vars,
                                           categorical_vars,
                                           focus_label,
                                           numeric_comparisons,
                                           categorical_comparisons,
                                           group_levels,
                                           sd_label,
                                           plot_prefix,
                                           output_dir) {
  subset_data <- data
  subset_data$group_label <- factor(subset_data$group_label, levels = group_levels)
  if (nrow(subset_data) == 0) {
    return(NULL)
  }

  numeric_results <- NULL
  sig_vars <- character(0)
  if (length(numeric_vars) > 0) {
    numeric_results <- .mlasdo_compute_numeric_tests(subset_data, numeric_vars, numeric_comparisons)
    sig_vars <- unique(numeric_results$variable[
      !is.na(numeric_results$p_value) & numeric_results$p_value < 0.05
    ])
    sig_vars <- unique(stats::na.omit(as.character(sig_vars)))
  }

  plot_path <- NULL
  if (length(sig_vars) > 0) {
    plot_path <- file.path(output_dir, paste0(plot_prefix, "_numeric_plot.png"))
    plot_path <- .mlasdo_save_anomaly_numeric_plot(
      subset_data = subset_data,
      variables = sig_vars,
      group_levels = group_levels,
      focus_label = focus_label,
      comparison_pairs = numeric_comparisons,
      sd_label = sd_label,
      output_path = plot_path
    )
  }

  categorical_results <- NULL
  categorical_sig <- NULL
  categorical_plot_path <- NULL
  if (length(categorical_vars) > 0) {
    categorical_results <- .mlasdo_compute_categorical_tests(
      subset_data,
      categorical_vars,
      categorical_comparisons
    )
    if (!is.null(categorical_results)) {
      categorical_sig <- categorical_results[!is.na(categorical_results$p_value) &
                                               categorical_results$p_value < 0.05, , drop = FALSE]
    }
    if (!is.null(categorical_sig) && nrow(categorical_sig) > 0) {
      categorical_plot_path <- file.path(output_dir, paste0(plot_prefix, "_categorical_plot.png"))
      categorical_plot_path <- .mlasdo_save_anomaly_categorical_plot(
        results_df = categorical_sig,
        focus_label = focus_label,
        output_path = categorical_plot_path
      )
    }
  }

  list(
    focus_label = focus_label,
    numeric_comparisons = numeric_comparisons,
    categorical_comparisons = categorical_comparisons,
    numeric_pvalues = numeric_results,
    categorical_pvalues = categorical_results,
    categorical_significant = categorical_sig,
    plot_path = plot_path,
    categorical_plot_path = categorical_plot_path
  )
}

.mlasdo_compute_numeric_tests <- function(data, numeric_vars, comparison_pairs) {
  res_list <- list()
  for (var in numeric_vars) {
    for (comp in comparison_pairs) {
      group_a <- comp[1]
      group_b <- comp[2]
      x <- data[[var]][data$group_label == group_a]
      y <- data[[var]][data$group_label == group_b]
      if (sum(!is.na(x)) > 1 && sum(!is.na(y)) > 1) {
        p_value <- stats::wilcox.test(x, y, exact = FALSE)$p.value
      } else {
        p_value <- NA_real_
      }
      res_list[[length(res_list) + 1L]] <- data.frame(
        variable = var,
        comparison = paste(group_a, "vs", group_b),
        p_value = p_value,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, res_list)
}

.mlasdo_compute_categorical_tests <- function(data, categorical_vars, comparison_pairs) {
  categorical_vars <- setdiff(categorical_vars, "group_label")
  if (length(categorical_vars) == 0) {
    return(NULL)
  }

  all_results <- list()
  for (var in categorical_vars) {
    var_data <- data[[var]]
    if (all(is.na(var_data))) next
    var_factor <- factor(var_data)
    for (comp in comparison_pairs) {
      subset_idx <- data$group_label %in% comp
      subset_factor <- droplevels(var_factor[subset_idx])
      subset_groups <- droplevels(data$group_label[subset_idx])
      if (length(subset_factor) == 0) next
      levels_var <- levels(subset_factor)
      for (lev in levels_var) {
        binary_var <- factor(
          ifelse(!is.na(subset_factor) & subset_factor == lev, lev, "rest"),
          levels = c(lev, "rest")
        )
        tbl <- table(subset_groups, binary_var, useNA = "no")
        if (nrow(tbl) < 2 || ncol(tbl) < 2) next
        chisq_res <- suppressWarnings(stats::chisq.test(tbl, correct = FALSE))
        expected <- chisq_res$expected
        if (any(expected < 5)) {
          fisher_res <- stats::fisher.test(tbl)
          p_value <- fisher_res$p.value
          method <- "Fisher's exact test"
        } else {
          p_value <- chisq_res$p.value
          method <- "Chi-square test"
        }
        group1 <- comp[1]
        group2 <- comp[2]
        row_levels_tbl <- rownames(tbl)
        col_levels_tbl <- colnames(tbl)
        row_idx_g1 <- match(group1, row_levels_tbl)
        row_idx_g2 <- match(group2, row_levels_tbl)
        col_idx_lev <- match(lev, col_levels_tbl)
        col_idx_rest <- match("rest", col_levels_tbl)
        count_lev_g1 <- if (!is.na(row_idx_g1) && !is.na(col_idx_lev)) tbl[row_idx_g1, col_idx_lev] else 0
        count_lev_g2 <- if (!is.na(row_idx_g2) && !is.na(col_idx_lev)) tbl[row_idx_g2, col_idx_lev] else 0
        count_rest_g1 <- if (!is.na(row_idx_g1) && !is.na(col_idx_rest)) tbl[row_idx_g1, col_idx_rest] else 0
        count_rest_g2 <- if (!is.na(row_idx_g2) && !is.na(col_idx_rest)) tbl[row_idx_g2, col_idx_rest] else 0
        all_results[[length(all_results) + 1L]] <- data.frame(
          variable = var,
          level = lev,
          comparison = paste(comp[1], "vs", comp[2]),
          group1 = comp[1],
          group2 = comp[2],
          method = method,
          p_value = p_value,
          level_count_group1 = count_lev_g1,
          level_count_group2 = count_lev_g2,
          rest_count_group1 = count_rest_g1,
          rest_count_group2 = count_rest_g2,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(all_results) == 0) {
    return(NULL)
  }
  do.call(rbind, all_results)
}

.mlasdo_save_anomaly_numeric_plot <- function(subset_data,
                                              variables,
                                              group_levels,
                                              focus_label,
                                              comparison_pairs,
                                              sd_label,
                                              output_path) {
  if (length(variables) == 0) {
    return(NULL)
  }
  color_map <- stats::setNames(
    c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"),
    group_levels
  )

  make_plot <- function(var_name) {
    df <- subset_data[, c("group_label", var_name), drop = FALSE]
    names(df) <- c("group_label", "value")
    df <- df[!is.na(df$value), , drop = FALSE]
    if (nrow(df) == 0) {
      return(NULL)
    }
    valid_comps <- Filter(
      function(comp) {
        sum(df$group_label == comp[1]) > 0 && sum(df$group_label == comp[2]) > 0
      },
      comparison_pairs
    )

    p <- ggpubr::ggboxplot(
      df,
      x = "group_label",
      y = "value",
      fill = "group_label",
      palette = color_map[group_levels],
      add = "jitter",
      notch = FALSE
    ) +
      ggpubr::rotate_x_text(angle = 45) +
      ggplot2::labs(
        title = var_name,
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold", size = 12)
      ) +
      ggplot2::scale_x_discrete(drop = FALSE, limits = group_levels)

    if (length(valid_comps) > 0) {
      p <- p + ggpubr::stat_compare_means(
        comparisons = valid_comps,
        method = "wilcox.test",
        label = "p.signif",
        hide.ns = TRUE
      )
    }
    p
  }

  plot_list <- lapply(variables, make_plot)
  plot_list <- Filter(Negate(is.null), plot_list)
  if (length(plot_list) == 0) {
    return(NULL)
  }

  n_plots <- length(plot_list)
  ncol <- min(3, max(1, ceiling(sqrt(n_plots))))
  nrow <- ceiling(n_plots / ncol)
  width <- max(6, ncol * 4.5)
  height <- max(5, nrow * 5)

  grid_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow)
  header <- paste(
    "Significant numerical covariates -", focus_label,
    "| Threshold:", sd_label
  )
  grid_plot <- ggpubr::annotate_figure(
    grid_plot,
    top = ggpubr::text_grob(header, face = "bold", size = 12)
  )

  ggplot2::ggsave(filename = output_path, plot = grid_plot, width = width, height = height, dpi = 300)
  normalizePath(output_path, winslash = "/", mustWork = FALSE)
}

.mlasdo_save_anomaly_categorical_plot <- function(results_df,
                                                  focus_label,
                                                  output_path) {
  if (is.null(results_df) || nrow(results_df) == 0) {
    return(NULL)
  }

  format_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 1e-4) format(p, scientific = TRUE, digits = 2) else sprintf("%.4f", p)
  }

  plot_entries <- lapply(seq_len(nrow(results_df)), function(i) {
    row <- results_df[i, ]
    facet_label <- sprintf(
      "%s = %s\n%s (p = %s)",
      row$variable,
      row$level,
      row$comparison,
      format_p(row$p_value)
    )
    data.frame(
      facet_label = facet_label,
      group_label = factor(
        c(row$group1, row$group1, row$group2, row$group2),
        levels = c(row$group1, row$group2)
      ),
      category = factor(
        c("Level", "Rest", "Level", "Rest"),
        levels = c("Level", "Rest")
      ),
      count = c(
        row$level_count_group1,
        row$rest_count_group1,
        row$level_count_group2,
        row$rest_count_group2
      ),
      stringsAsFactors = FALSE
    )
  })

  tile_df <- do.call(rbind, plot_entries)
  tile_df$count[is.na(tile_df$count)] <- 0

  n_facets <- length(unique(tile_df$facet_label))
  ncol <- min(3, max(1, ceiling(sqrt(n_facets))))
  nrow <- ceiling(n_facets / ncol)

  p <- ggplot2::ggplot(tile_df, ggplot2::aes(x = category, y = group_label, fill = count)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = count), fontface = "bold", size = 3) +
    ggplot2::scale_fill_gradient(low = "#fef0d9", high = "#b30000") +
    ggplot2::facet_wrap(~ facet_label, scales = "fixed", ncol = 2) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 6),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 9),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", size = 11)
    ) +
    ggplot2::labs(title = paste("Categorical tests -", focus_label), fill = "Count") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"))

  width <- max(5, ncol * 3.5)
  height <- max(4, nrow * 3.5)
  ggplot2::ggsave(filename = output_path, plot = p, width = width, height = height, dpi = 300)
  normalizePath(output_path, winslash = "/", mustWork = FALSE)
}

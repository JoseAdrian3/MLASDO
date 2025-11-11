.mlasdo_detect_anomalies <- function(distances,
                                     method = c("sd", "zscore"),
                                     sd_multiplier = 2,
                                     z_threshold = 1.65,
                                     output_dir,
                                     positive_class = NULL,
                                     negative_class = NULL) {
  method <- match.arg(method)
  if (missing(distances)) {
    stop("Argument 'distances' must be provided.")
  }

  if (is.character(distances) && length(distances) == 1) {
    distances <- readRDS(distances)
  }
  distances <- as.data.frame(distances)

  required_cols <- c("id", "diagnosis_true", "diagnosis_pred", "signed_margin")
  missing_cols <- setdiff(required_cols, names(distances))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Distances object is missing required column(s): %s",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  if (!"fold" %in% names(distances)) {
    distances$fold <- NA_integer_
  }

  distances$diagnosis_true <- droplevels(factor(distances$diagnosis_true))
  distances$diagnosis_pred <- droplevels(factor(distances$diagnosis_pred))
  class_levels <- levels(distances$diagnosis_true)
  if (length(class_levels) < 2) {
    stop("Anomaly detection requires at least two classes in 'diagnosis_true'.")
  }

  stats_list <- lapply(class_levels, function(cls) {
    idx <- distances$diagnosis_true == cls
    vals <- distances$signed_margin[idx]
    data.frame(
      diagnosis_true = cls,
      mean_margin = mean(vals, na.rm = TRUE),
      sd_margin = stats::sd(vals, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  stats_df <- do.call(rbind, stats_list)
  stats_df$sd_margin[is.na(stats_df$sd_margin)] <- 0

  stats_df$role <- "other"
  if (!is.null(positive_class) && positive_class %in% stats_df$diagnosis_true) {
    stats_df$role[stats_df$diagnosis_true == positive_class] <- "positive"
  }
  if (!is.null(negative_class) && negative_class %in% stats_df$diagnosis_true) {
    stats_df$role[stats_df$diagnosis_true == negative_class] <- "negative"
  }
  if (all(stats_df$role == "other") && nrow(stats_df) >= 2) {
    max_mean <- max(stats_df$mean_margin, na.rm = TRUE)
    min_mean <- min(stats_df$mean_margin, na.rm = TRUE)
    stats_df$role <- ifelse(stats_df$mean_margin == max_mean, "positive",
                            ifelse(stats_df$mean_margin == min_mean, "negative", "other"))
  } else if (all(stats_df$role == "other") && nrow(stats_df) == 1) {
    stats_df$role <- "positive"
  }

  resolved_positive <- positive_class
  resolved_negative <- negative_class
  if (is.null(resolved_positive)) {
    pos_candidate <- stats_df$diagnosis_true[stats_df$role == "positive"]
    if (length(pos_candidate) == 0) {
      pos_candidate <- stats_df$diagnosis_true
    }
    resolved_positive <- pos_candidate[1]
  }
  if (is.null(resolved_negative)) {
    neg_candidate <- stats_df$diagnosis_true[stats_df$role == "negative" & stats_df$diagnosis_true != resolved_positive]
    if (length(neg_candidate) == 0) {
      neg_candidate <- setdiff(stats_df$diagnosis_true, resolved_positive)
    }
    if (length(neg_candidate) > 0) {
      resolved_negative <- neg_candidate[1]
    } else {
      resolved_negative <- resolved_positive
    }
  }
  stats_df$role[stats_df$diagnosis_true == resolved_positive] <- "positive"
  stats_df$role[stats_df$diagnosis_true == resolved_negative] <- "negative"

  multiplier <- if (method == "sd") sd_multiplier else z_threshold
  stats_df$upper_threshold <- NA_real_
  stats_df$lower_threshold <- NA_real_
  neg_idx_stats <- stats_df$role == "negative"
  pos_idx_stats <- stats_df$role == "positive"
  stats_df$upper_threshold[neg_idx_stats] <- stats_df$mean_margin[neg_idx_stats] +
    sd_multiplier * stats_df$sd_margin[neg_idx_stats]
  stats_df$lower_threshold[pos_idx_stats] <- stats_df$mean_margin[pos_idx_stats] -
    sd_multiplier * stats_df$sd_margin[pos_idx_stats]

  distances$class_mean <- stats_df$mean_margin[match(distances$diagnosis_true, stats_df$diagnosis_true)]
  distances$class_sd <- stats_df$sd_margin[match(distances$diagnosis_true, stats_df$diagnosis_true)]
  valid_sd <- !is.na(distances$class_sd) & distances$class_sd > 0
  distances$z_score <- NA_real_
  distances$z_score[valid_sd] <- (distances$signed_margin[valid_sd] - distances$class_mean[valid_sd]) /
    distances$class_sd[valid_sd]

  distances$class_role <- stats_df$role[match(distances$diagnosis_true, stats_df$diagnosis_true)]
  distances$class_upper <- stats_df$upper_threshold[match(distances$diagnosis_true, stats_df$diagnosis_true)]
  distances$class_lower <- stats_df$lower_threshold[match(distances$diagnosis_true, stats_df$diagnosis_true)]

  flag <- rep(FALSE, nrow(distances))
  pos_idx <- distances$class_role == "positive"
  neg_idx <- distances$class_role == "negative"
  if (identical(method, "sd")) {
    threshold_pos <- distances$signed_margin < distances$class_lower
    threshold_neg <- distances$signed_margin > distances$class_upper
  } else {
    threshold_pos <- distances$z_score < -z_threshold
    threshold_neg <- distances$z_score > z_threshold
  }
  misclassified_pos <- distances$diagnosis_pred == resolved_negative
  misclassified_neg <- distances$diagnosis_pred == resolved_positive
  flag[pos_idx] <- threshold_pos[pos_idx] & misclassified_pos[pos_idx]
  flag[neg_idx] <- threshold_neg[neg_idx] & misclassified_neg[neg_idx]
  flag[is.na(flag)] <- FALSE

  anomalies <- distances[flag, , drop = FALSE]
  anomaly_counts <- data.frame(
    diagnosis_true = class_levels,
    anomalous = vapply(class_levels, function(cls) sum(anomalies$diagnosis_true == cls), integer(1)),
    stringsAsFactors = FALSE
  )

  positive_label <- if (!is.null(resolved_positive)) resolved_positive else "Positive"
  negative_label <- if (!is.null(resolved_negative)) resolved_negative else "Negative"
  distances$class_status <- ifelse(
    as.character(distances$diagnosis_true) == as.character(distances$diagnosis_pred),
    "Correctly classified",
    "Misclassified"
  )
  distances$group <- with(
    distances,
    ifelse(
      diagnosis_true == resolved_positive & class_status == "Correctly classified",
      sprintf("Correctly classified %s", positive_label),
      ifelse(
        diagnosis_true == resolved_positive & class_status == "Misclassified",
        sprintf("Misclassified %s", positive_label),
        ifelse(
          diagnosis_true == resolved_negative & class_status == "Correctly classified",
          sprintf("Correctly classified %s", negative_label),
          sprintf("Misclassified %s", negative_label)
        )
      )
    )
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  anomalies_path <- file.path(output_dir, "anomalous_samples.rds")
  thresholds_path <- file.path(output_dir, "anomaly_thresholds.rds")
  counts_path <- file.path(output_dir, "anomaly_counts.rds")

  saveRDS(anomalies, anomalies_path)
  saveRDS(stats_df, thresholds_path)
  saveRDS(anomaly_counts, counts_path)

  plot_path <- file.path(output_dir, "anomaly_distribution.png")
  plot_saved <- .mlasdo_save_anomaly_plot(
    distances_df = distances,
    thresholds_df = stats_df,
    method = method,
    multiplier = multiplier,
    plot_path = plot_path,
    positive_class = resolved_positive,
    negative_class = resolved_negative
  )

  list(
    method = method,
    sd_multiplier = sd_multiplier,
    z_threshold = z_threshold,
    thresholds = stats_df,
    counts = anomaly_counts,
    anomalies = anomalies,
    plot_path = if (isTRUE(plot_saved)) normalizePath(plot_path, winslash = "/", mustWork = FALSE) else NULL,
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    anomalies_path = normalizePath(anomalies_path, winslash = "/", mustWork = FALSE),
    thresholds_path = normalizePath(thresholds_path, winslash = "/", mustWork = FALSE),
    counts_path = normalizePath(counts_path, winslash = "/", mustWork = FALSE),
    positive_class = resolved_positive,
    negative_class = resolved_negative
  )
}

.mlasdo_save_anomaly_plot <- function(distances_df,
                                      thresholds_df,
                                      method,
                                      multiplier,
                                      plot_path,
                                      positive_class = NULL,
                                      negative_class = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is required to generate the anomaly plot.")
    return(FALSE)
  }

  binwidth <- diff(range(distances_df$signed_margin, na.rm = TRUE))
  if (!is.finite(binwidth) || binwidth <= 0) {
    binwidth <- 0.1
  } else {
    binwidth <- binwidth / 50
  }

  lines_df <- data.frame()
  if (any(thresholds_df$role == "negative")) {
    neg_rows <- thresholds_df[thresholds_df$role == "negative", , drop = FALSE]
    lines_df <- rbind(lines_df, data.frame(
      diagnosis_true = neg_rows$diagnosis_true,
      threshold_value = neg_rows$upper_threshold,
      label = sprintf("%s upper threshold", neg_rows$diagnosis_true),
      stringsAsFactors = FALSE
    ))
  }
  if (any(thresholds_df$role == "positive")) {
    pos_rows <- thresholds_df[thresholds_df$role == "positive", , drop = FALSE]
    lines_df <- rbind(lines_df, data.frame(
      diagnosis_true = pos_rows$diagnosis_true,
      threshold_value = pos_rows$lower_threshold,
      label = sprintf("%s lower threshold", pos_rows$diagnosis_true),
      stringsAsFactors = FALSE
    ))
  }

  pos_label <- positive_class
  neg_label <- negative_class
  if (is.null(pos_label) && any(thresholds_df$role == "positive")) {
    pos_label <- thresholds_df$diagnosis_true[thresholds_df$role == "positive"][1]
  }
  if (is.null(neg_label) && any(thresholds_df$role == "negative")) {
    neg_label <- thresholds_df$diagnosis_true[thresholds_df$role == "negative"][1]
  }
  if (is.null(pos_label)) {
    pos_label <- setdiff(unique(distances_df$diagnosis_true), neg_label)[1]
  }
  if (is.null(neg_label)) {
    neg_label <- setdiff(unique(distances_df$diagnosis_true), pos_label)[1]
  }
  if (is.null(pos_label)) pos_label <- "Positive"
  if (is.null(neg_label)) neg_label <- "Negative"

  palette <- c(
    stats::setNames("#1B9E77", sprintf("Correctly classified %s", pos_label)),
    stats::setNames("#D95F02", sprintf("Correctly classified %s", neg_label)),
    stats::setNames("#7570B3", sprintf("Misclassified %s", pos_label)),
    stats::setNames("#E6AB02", sprintf("Misclassified %s", neg_label))
  )
  groups <- unique(distances_df$group)
  missing_groups <- setdiff(groups, names(palette))
  if (length(missing_groups) > 0) {
    extra_cols <- grDevices::rainbow(length(missing_groups))
    names(extra_cols) <- missing_groups
    palette <- c(palette, extra_cols)
  }
  ordered_groups <- names(palette)
  distances_df$group <- factor(distances_df$group, levels = ordered_groups)

  if (method == "sd") {
    multiplier_label <- sprintf("mean ± %.2f SD", multiplier)
    subtitle_label <- sprintf("Anomaly method: mean ± %.2f SD", multiplier)
  } else {
    multiplier_label <- sprintf("|z| > %.2f", multiplier)
    subtitle_label <- sprintf("Anomaly method: |z| > %.2f", multiplier)
  }

  title_label <- "Distribution of signed distances to the SVM hyperplane"

  p <- ggplot2::ggplot(distances_df, ggplot2::aes(x = signed_margin, fill = group)) +
    ggplot2::geom_histogram(binwidth = binwidth, color = "black", alpha = 0.4, position = "identity") +
    ggplot2::labs(
      title = title_label,
      subtitle = subtitle_label,
      x = "Signed distance to the hyperplane",
      y = "Frequency",
      fill = "Group",
      color = "Threshold"
    ) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.3),
      legend.direction = "vertical",
      plot.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1, byrow = TRUE))

  if (nrow(lines_df) > 0) {
    p <- p +
      ggplot2::geom_vline(
        data = lines_df,
        ggplot2::aes(xintercept = threshold_value, color = label),
        linewidth = 0.9,
        alpha = 0.8
      ) +
      ggplot2::scale_color_manual(
        values = rep("#2B8CBE", nrow(lines_df)),
        labels = rep(multiplier_label, nrow(lines_df))
      )
  }

  ggplot2::ggsave(plot_path, plot = p, width = 10, height = 6, dpi = 300)
  TRUE
}

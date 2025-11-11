.mlasdo_source_label <- function(x) {
  if (is.character(x) && length(x) == 1) {
    normalizePath(x, winslash = "/", mustWork = FALSE)
  } else {
    "in-memory object"
  }
}

.mlasdo_prepare_directory <- function(path) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.mlasdo_prepare_subdir <- function(parent_dir, subdir) {
  dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(parent_dir, subdir)
  .mlasdo_prepare_directory(path)
}

.mlasdo_save_preprocessing_data <- function(preprocess_dir, omic, clinical, verbose) {
  files <- list(
    omic = file.path(preprocess_dir, "omic_aligned.rds"),
    clinical = file.path(preprocess_dir, "clinical_aligned.rds")
  )
  saveRDS(omic, files$omic)
  saveRDS(clinical, files$clinical)
  if (isTRUE(verbose)) {
    message(sprintf("Aligned datasets saved under %s", preprocess_dir))
  }
  files
}

.mlasdo_save_metadata <- function(output_dir, metadata) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  metadata_path <- file.path(output_dir, "run_metadata.rds")
  saveRDS(metadata, metadata_path)
  normalizePath(metadata_path, winslash = "/", mustWork = TRUE)
}

.mlasdo_stop_if_missing <- function(x, label) {
  if (missing(x) || length(x) == 0 || is.null(x)) {
    stop(sprintf("Argument '%s' must be provided.", label))
  }
}

.mlasdo_select_cores <- function(requested) {
  available <- max(1L, parallel::detectCores() - 1L)
  if (is.null(requested)) {
    return(available)
  }
  as.integer(max(1L, min(requested, available)))
}

.mlasdo_copy_directory <- function(src, dst) {
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(src, full.names = TRUE, all.files = FALSE, recursive = FALSE)
  if (length(files) == 0) {
    return(invisible(dst))
  }
  file.copy(files, dst, recursive = TRUE, overwrite = TRUE)
  invisible(dst)
}

`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs) && length(lhs) > 0) {
    lhs
  } else {
    rhs
  }
}

.mlasdo_log_step <- function(fmt, ...) {
  msg <- if (missing(fmt)) "" else sprintf(fmt, ...)
  cat(sprintf("[MLASDO] %s\n", msg))
}

.mlasdo_configure_target <- function(vec, positive_class = NULL, negative_class = NULL) {
  if (!is.factor(vec)) {
    vec <- factor(vec)
  }
  lvls <- levels(vec)
  if (!is.null(positive_class) && !positive_class %in% lvls) {
    stop(sprintf("Positive class '%s' not found in target levels.", positive_class))
  }
  if (!is.null(negative_class) && !negative_class %in% lvls) {
    stop(sprintf("Negative class '%s' not found in target levels.", negative_class))
  }
  if (!is.null(positive_class) && !is.null(negative_class) &&
      positive_class == negative_class) {
    stop("Positive and negative classes must be different.")
  }
  resolved_positive <- positive_class
  resolved_negative <- negative_class
  if (is.null(resolved_positive) && is.null(resolved_negative)) {
    return(list(values = vec, positive = NULL, negative = NULL))
  }
  if (is.null(resolved_positive)) {
    resolved_positive <- setdiff(lvls, resolved_negative)[1]
  }
  if (is.null(resolved_negative)) {
    resolved_negative <- setdiff(lvls, resolved_positive)[1]
  }
  ordered_levels <- unique(c(resolved_negative, resolved_positive, lvls))
  list(
    values = factor(vec, levels = ordered_levels),
    positive = resolved_positive,
    negative = resolved_negative
  )
}

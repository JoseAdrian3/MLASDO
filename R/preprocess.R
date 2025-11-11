#' Prepare aligned datasets for MLASDO
#'
#' Centralises all the logic that loads the omic/clinical inputs, validates the
#' ID/target columns, and produces two aligned data frames ready for the rest of
#' the pipeline.
#'
#' @return A list with `omic`, `clinical`, and `sources` entries.
prepare_mlasdo_data <- function(omic_data,
                                clinical_data,
                                id_column,
                                target_column,
                                verbose = interactive()) {
  .mlasdo_stop_if_missing(id_column, "id_column")
  .mlasdo_stop_if_missing(target_column, "target_column")

  omic_df <- .mlasdo_as_data_frame(omic_data, "omic_data")
  clinical_df <- .mlasdo_as_data_frame(clinical_data, "clinical_data")

  .mlasdo_require_columns(omic_df, id_column, "omic dataset")
  .mlasdo_require_columns(clinical_df, id_column, "clinical dataset")

  target_in_omic <- target_column %in% names(omic_df)
  target_in_clinical <- target_column %in% names(clinical_df)
  if (!target_in_omic && !target_in_clinical) {
    stop(
      sprintf(
        "Column '%s' not found in either dataset. Provide a valid target_column.",
        target_column
      )
    )
  }

  alignment <- .mlasdo_align_by_id(omic_df, clinical_df, id_column)
  omic_aligned <- alignment$omic
  clinical_aligned <- alignment$clinical

  primary_target <- NULL
  secondary_target <- NULL
  prefer_clinical <- target_in_clinical && target_column %in% names(clinical_aligned)
  if (prefer_clinical) {
    primary_target <- clinical_aligned[[target_column]]
    if (target_in_omic && target_column %in% names(omic_aligned)) {
      secondary_target <- omic_aligned[[target_column]]
    }
  } else if (target_column %in% names(omic_aligned)) {
    primary_target <- omic_aligned[[target_column]]
    if (target_in_clinical && target_column %in% names(clinical_aligned)) {
      secondary_target <- clinical_aligned[[target_column]]
    }
  }

  if (is.null(primary_target)) {
    stop(
      sprintf(
        "Unable to locate column '%s' after aligning the datasets.",
        target_column
      )
    )
  }

  if (!is.null(secondary_target)) {
    mismatch_idx <- which(
      !is.na(primary_target) &
        !is.na(secondary_target) &
        as.character(primary_target) != as.character(secondary_target)
    )
    if (length(mismatch_idx) > 0 && prefer_clinical && isTRUE(verbose)) {
      message(
        sprintf(
          "Detected %d '%s' mismatches; keeping values from the clinical dataset.",
          length(mismatch_idx),
          target_column
        )
      )
    }
    missing_idx <- which(is.na(primary_target) & !is.na(secondary_target))
    if (length(missing_idx) > 0) {
      primary_target[missing_idx] <- secondary_target[missing_idx]
    }
  }

  if (!is.factor(primary_target)) {
    primary_target <- factor(primary_target)
  }
  omic_aligned[[target_column]] <- primary_target
  clinical_aligned[[target_column]] <- primary_target

  shared_columns <- intersect(names(omic_aligned), names(clinical_aligned))
  shared_columns <- setdiff(shared_columns, c(id_column, target_column))
  if (length(shared_columns) > 0) {
    clinical_aligned <- clinical_aligned[
      setdiff(names(clinical_aligned), shared_columns)
    ]
  }

  list(
    omic = omic_aligned,
    clinical = clinical_aligned,
    sources = list(
      omic = .mlasdo_source_label(omic_data),
      clinical = .mlasdo_source_label(clinical_data)
    )
  )
}

.mlasdo_as_data_frame <- function(x, label) {
  if (is.data.frame(x)) {
    return(x)
  }

  if (length(x) != 1 || !is.character(x)) {
    stop(sprintf("Argument '%s' must be a data frame or a single file path.", label))
  }

  path <- normalizePath(x, winslash = "/", mustWork = TRUE)
  ext <- tolower(tools::file_ext(path))

  if (ext == "rds") {
    obj <- readRDS(path)
  } else {
    stop(sprintf("Unsupported file type for '%s': %s", label, ext))
  }

  if (!is.data.frame(obj)) {
    stop(sprintf("Object loaded from '%s' is not a data frame.", path))
  }
  obj
}

.mlasdo_require_columns <- function(data, columns, dataset_label) {
  missing_cols <- setdiff(columns, names(data))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Column(s) %s not found in %s.",
        paste(shQuote(missing_cols), collapse = ", "),
        dataset_label
      )
    )
  }
}

.mlasdo_align_by_id <- function(omic_df, clinical_df, id_column) {
  ids_omic <- omic_df[[id_column]]
  ids_clinic <- clinical_df[[id_column]]

  if (anyDuplicated(ids_omic)) {
    stop(sprintf("Duplicate IDs detected in omic dataset for column '%s'.", id_column))
  }

  clinic_unique <- clinical_df[!duplicated(ids_clinic), , drop = FALSE]
  match_idx <- match(ids_omic, clinic_unique[[id_column]])
  if (anyNA(match_idx)) {
    missing_ids <- unique(ids_omic[is.na(match_idx)])
    stop(
      sprintf(
        "The following IDs are missing in the clinical dataset: %s",
        paste(head(missing_ids, 10), collapse = ", ")
      )
    )
  }

  clinical_ordered <- clinic_unique[match_idx, , drop = FALSE]
  list(omic = omic_df, clinical = clinical_ordered)
}

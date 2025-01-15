#' Perform Statistical Analysis for User Variable Classification
#'
#' This function performs statistical analysis to evaluate the relationships between predictors and user-defined classes.
#' It calculates odds ratios, statistical tests, and Wilcoxon tests for categorical and numeric predictors.
#'
#' @param saving_name Name to save the output files.
#' @param changed_clinic_data Modified clinic dataset for analysis.
#' @param first_group Name of the first group for classification.
#' @param second_group Name of the second group for classification.
#' @param active_predictors List of all predictors used in the analysis.
#' @param categoric_active_predictors List of categorical predictors.
#' @param numeric_active_predictors List of numeric predictors.
#' @param class_variable Name of the variable defining the classification groups.
#' @param default_odds_predictors Predictors for which default odds ratio calculations apply.
#' @param method Statistical method to apply (e.g., Fisher's test, Chi-square).
#' @param dir_path Path to the directory for saving results.
#'
#' @return Saves results including odds ratios and Wilcoxon test results to specified directory.
#' @export
perform_statistic_analysis_user_variable_classification <- function(
    saving_name,
    changed_clinic_data,
    first_group,
    second_group,
    active_predictors,
    categoric_active_predictors,
    numeric_active_predictors,
    class_variable,
    default_odds_predictors,
    method,
    dir_path
) {
  # DATA READING
  clinic <- changed_clinic_data

  # DATA PROCESSING: Replace NA values in categorical predictors with "Empty"
  for (active_predictor in categoric_active_predictors) {
    clinic[[active_predictor]] <- replace(clinic[[active_predictor]], is.na(clinic[[active_predictor]]), "Empty")
  }

  # PATIENT GROUP CREATION
  c_prima <- subset(clinic, clinic[[class_variable]] == "Control")
  c_prima_prima <- subset(clinic, clinic[[class_variable]] == first_group)
  e_prima <- subset(clinic, clinic[[class_variable]] == "Case")
  e_prima_prima <- subset(clinic, clinic[[class_variable]] == second_group)

  c <- rbind(c_prima, e_prima_prima)
  e <- rbind(e_prima, c_prima_prima)

  # Initialize data frames for results
  total_odds_controls <- data.frame(
    OddRatio = numeric(),
    Fisher_chi.square = numeric(),
    Feature = character(),
    FeatureValue = character(),
    FeatureRestValues = character(),
    SumCCFirst = integer(),
    SumCCSecond = integer(),
    SumCMFirst = integer(),
    SumCMSecond = integer()
  )
  total_odds_cases <- data.frame(
    OddRatio = numeric(),
    Fisher_chi.square = numeric(),
    Feature = character(),
    FeatureValue = character(),
    FeatureRestValues = character(),
    SumCCFirst = integer(),
    SumCCSecond = integer(),
    SumCMFirst = integer(),
    SumCMSecond = integer()
  )

  total_wilcox_controls <- data.frame(value = numeric(0), group = numeric(0), predictor = character(0))
  total_wilcox_cases <- data.frame(value = numeric(0), group = numeric(0), predictor = character(0))

  # ANALYSIS OF PREDICTORS
  for (active_predictor in active_predictors) {
    if (active_predictor %in% numeric_active_predictors) {
      # Handle numeric predictors and compute Wilcoxon test data

      num_patients_cc <- nrow(e_prima_prima)
      cc <- if (num_patients_cc == 0) 0 else e_prima_prima[[active_predictor]]

      num_patients_cm <- nrow(c_prima)
      cm <- if (num_patients_cm == 0) 0 else c_prima[[active_predictor]]

      num_patients_ec <- nrow(c_prima_prima)
      ec <- if (num_patients_ec == 0) 0 else c_prima_prima[[active_predictor]]

      num_patients_em <- nrow(e_prima)
      em <- if (num_patients_em == 0) 0 else e_prima[[active_predictor]]

      # Add data to controls and cases
      total_wilcox_controls <- rbind(total_wilcox_controls, data.frame(value = cc, group = "Controls2Cases", predictor = active_predictor))
      total_wilcox_controls <- rbind(total_wilcox_controls, data.frame(value = cm, group = "Controls", predictor = active_predictor))
      total_wilcox_cases <- rbind(total_wilcox_cases, data.frame(value = ec, group = "Cases2Controls", predictor = active_predictor))
      total_wilcox_cases <- rbind(total_wilcox_cases, data.frame(value = em, group = "Cases", predictor = active_predictor))

    } else {
      # Handle categorical predictors and compute odds ratios
      levels_categorical <- levels(clinic[[active_predictor]])

      for (level in levels_categorical) {
        sum_cc_first <- sum(e_prima_prima[[active_predictor]] == level, na.rm = TRUE)
        sum_cc_second <- sum(e_prima_prima[[active_predictor]] != level, na.rm = TRUE)
        sum_cm_first <- sum(c_prima[[active_predictor]] == level, na.rm = TRUE)
        sum_cm_second <- sum(c_prima[[active_predictor]] != level, na.rm = TRUE)

        contingency_table <- matrix(c(sum_cc_first, sum_cc_second, sum_cm_first, sum_cm_second), nrow = 2)

        test_result <- if (any(contingency_table < 10)) fisher.test(contingency_table) else chisq.test(contingency_table)
        odd_ratio <- if (!is.null(test_result$estimate)) test_result$estimate else NA

        total_odds_controls <- rbind(total_odds_controls, data.frame(
          OddRatio = odd_ratio,
          Fisher_chi.square = test_result$p.value,
          Feature = active_predictor,
          FeatureValue = level,
          FeatureRestValues = paste(setdiff(levels_categorical, level), collapse = ","),
          SumCCFirst = sum_cc_first,
          SumCCSecond = sum_cc_second,
          SumCMFirst = sum_cm_first,
          SumCMSecond = sum_cm_second
        ))
      }
    }
  }

  # Saving Results
  dir_path <- paste(saving_name, dir_path, sep = "/")

  dir_path <- paste(dir_path, saving_name, sep = "/")

  gaPathOddsControls <- paste(dir_path, "OddsRatiosControls.rds", sep="_")
  gaPathOddsCases <- paste(dir_path, "OddsRatiosCases.rds", sep="_")

  saveRDS(total_odds_controls, gaPathOddsControls)
  saveRDS(total_odds_cases, gaPathOddsCases)

  gaPathWilcoxControls <- paste(dir_path, "WilcoxControls.rds", sep="_")
  saveRDS(total_wilcox_controls, gaPathWilcoxControls)

  gaPathWilcoxCases <- paste(dir_path, "WilcoxCases.rds", sep="_")
  saveRDS(total_wilcox_cases, gaPathWilcoxCases)
}

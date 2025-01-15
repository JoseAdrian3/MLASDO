#' Statistical Calculations for Omic and Clinical Data
#'
#' This function performs data processing and statistical analysis on omic and clinical datasets.
#'
#' @param method The statistical method to be used.
#' @param omic_data Data frame containing omic data.
#' @param active_predictors Vector of active predictor variable names.
#' @param class_variable Name of the class variable.
#' @param clinic_data Data frame containing clinical data.
#' @param id_column Column name used to identify rows.
#' @param solution Binary vector representing a solution from a genetic algorithm.
#' @param categoric_active_predictors Optional vector of categorical active predictors.
#' @param numeric_active_predictors Optional vector of numeric active predictors.
#' @param default_odds_predictors Default predictors for odds ratio analysis.
#' @param saving_name Name used for saving analysis results.
#' @param dir_path Directory path for saving outputs.
#' @param seed Random seed for reproducibility.
#' @export
statistical_calculations <- function(method,
                                     omic_data,
                                     active_predictors,
                                     class_variable,
                                     clinic_data,
                                     id_column,
                                     solution,
                                     categoric_active_predictors = NULL,
                                     numeric_active_predictors = NULL,
                                     default_odds_predictors,
                                     saving_name,
                                     dir_path,
                                     seed) {

  # Set the random seed for reproducibility
  set.seed(seed)

  #### DATA PROCESSING ####

  # Remove active predictors from the omic data
  omic_data <- omic_data[, !(names(omic_data) %in% active_predictors)]

  # Filter out the class variable from the omic data
  omics_filtered <- omic_data[, !(names(omic_data) %in% class_variable)]

  # Merge clinical data with the filtered omic data using the ID column
  valid_clinic_data <- merge(clinic_data, omics_filtered, by = id_column)

  # Retain only the active predictors, ID column, and class variable in the clinical data
  valid_clinic_data <- valid_clinic_data[, (names(valid_clinic_data) %in% c(active_predictors, id_column, class_variable))]

  # Copy the original class variable for further modifications
  changed_diagnoses <- omic_data[[class_variable]]

  # Create strings representing diagnostic changes
  first_group <- paste("Case", "Control", sep = "2")
  second_group <- paste("Control", "Case", sep = "2")

  # Convert original diagnoses to numeric values (1 for "Case", 0 for "Control")
  changed_diagnoses <- ifelse(changed_diagnoses == "Case", 1, 0)

  # Apply the genetic algorithm solution using a bitwise XOR operation
  changed_diagnoses <- bitwXor(changed_diagnoses, solution)

  # Get indices where changes are indicated by the solution
  change_indices <- which(solution == 1)

  # Iterate over the change indices to modify the diagnoses
  for (i in change_indices) {
    if (changed_diagnoses[[i]] == 1) {
      changed_diagnoses[[i]] <- second_group
    } else {
      changed_diagnoses[[i]] <- first_group
    }
  }

  # Revert numeric diagnoses back to textual representation
  changed_diagnoses <- ifelse(changed_diagnoses == 0, "Control", changed_diagnoses)
  changed_diagnoses <- ifelse(changed_diagnoses == 1, "Case", changed_diagnoses)

  # Update the clinical data with the changed diagnoses
  changed_clinic_data <- valid_clinic_data
  changed_clinic_data[[class_variable]] <- changed_diagnoses

  #### STATISTICAL ANALYSIS ####

  # Perform statistical analysis using the MLASDO package
  MLASDO::perform_statistic_analysis_user_variable_classification(
    saving_name = saving_name,
    changed_clinic_data = changed_clinic_data,
    first_group = first_group,
    second_group = second_group,
    active_predictors = active_predictors,
    categoric_active_predictors = categoric_active_predictors,
    numeric_active_predictors = numeric_active_predictors,
    class_variable = class_variable,
    default_odds_predictors = default_odds_predictors,
    method = method,
    dir_path = dir_path
  )
}

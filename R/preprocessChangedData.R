#' Preprocess Changed Data
#'
#' This function preprocesses omic and clinical data based on the specified solution.
#' It also modifies the class variable and creates an additional dataset
#' `all_changed_omic_data` where class changes apply only to the indices
#' provided in `svm_data$selected_samples$index`.
#'
#' @param omic Data frame containing omic data.
#' @param clinic_data Data frame containing clinical data.
#' @param solution Integer indicating the solution type (e.g., 1 for swapping classes).
#' @param active_predictors Vector of predictor variable names to exclude from the omic data.
#' @param class_variable String of the class variable.
#' @param id_column String of the IDs variable.
#' @param seed Seed for reproducibility.
#' @param svm_data List with SVM-related information, including `selected_samples$index`.
#'
#' @return A list containing four data frames:
#'   - `changed_omic_data`: Omic dataset modified according to `solution`.
#'   - `changed_clinic_data`: Clinical dataset modified according to `solution`.
#'   - `transition_clinic_data`: Clinical dataset with transition labels (e.g., "Case2Control").
#'   - `all_changed_omic_data`: Omic dataset where classes are changed only for the
#'     indices specified in `svm_data$selected_samples$index`.
#'
#' @export
preprocess_changed_data <- function(omic,
                                    clinic_data,
                                    solution,
                                    active_predictors,
                                    class_variable,
                                    id_column,
                                    seed,
                                    svm_data) {

  # Set the random seed for reproducibility
  set.seed(seed)

  # Store the original class labels as characters
  original_diagnoses <- omic[[class_variable]]
  original_diagnoses_chr <- as.character(original_diagnoses)

  # 1) Remove active predictors from the omic dataset
  omic_id <- omic[[id_column]]
  omic <- omic[, !(names(omic) %in% active_predictors)]

  # 2) Merge to ensure "clinic_data" and "omic" are aligned
  # (only keeping id_column, class_variable, and active_predictors)
  valid_clinic_data <- merge(clinic_data, omic, by = id_column)
  vars_to_keep <- c(active_predictors, id_column, class_variable)
  valid_clinic_data <- valid_clinic_data[, names(valid_clinic_data) %in% vars_to_keep]

  # 3) Modify classes based on "solution"
  # (For instance, solution == 1 => swap "Case" <-> "Control")
  changed_diagnoses <- ifelse(
    solution == 1 & original_diagnoses_chr == "Case", "Control",
    ifelse(solution == 1 & original_diagnoses_chr == "Control", "Case", original_diagnoses_chr)
  )

  # Create omic and clinic data frames with modified classes
  changed_omic_data <- omic
  changed_omic_data[[class_variable]] <- changed_diagnoses

  changed_clinic_data <- valid_clinic_data
  changed_clinic_data[[class_variable]] <- changed_diagnoses

  # Create transition labels
  transition_labels <- ifelse(
    solution == 1 & original_diagnoses_chr == "Case", "Case2Control",
    ifelse(solution == 1 & original_diagnoses_chr == "Control", "Control2Case", original_diagnoses_chr)
  )

  transition_clinic_data <- valid_clinic_data
  transition_clinic_data[[class_variable]] <- transition_labels

  # 4) Return all modified datasets in a list
  list(
    changed_omic_data = changed_omic_data,
    changed_clinic_data = changed_clinic_data,
    transition_clinic_data = transition_clinic_data
  )
}

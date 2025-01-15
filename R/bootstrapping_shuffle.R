#' Bootstrapping Shuffle Function
#'
#' Perform a bootstrapping-based shuffle process.
#' The function performs the following steps:
#' - Normalizes the distances
#' - Repeats the process 100,000 times, assigning random labels between Control and Case
#'   and calculates the threshold of that distribution
#' - Computes the average of all thresholds
#' - Selects the anomalous samples based on the average of the thresholds
#'
#' @param clinic_data Clinical dataset.
#' @param svm_data Transcriptomic dataset containing misclassified samples and their distances.
#' @param seed Seed for reproducibility.
#' @param saving_name String representing the base name for saving output files.
#'
#' @return The function saves bootstrapping results, thresholds,
#' and selected samples in the specified output directory.
#'
#' @export
bootstrapping_shuffle <- function(clinic_data, svm_data, seed, saving_name) {

  # --- Step 1: Extract and integrate data ---------------------------------------
  # Extract misclassified samples and combine them with clinical data.
  # Add an index column to the clinical dataset for later reference.
  missclassified_samples <- svm_data$selected_samples
  clinic_data$index <- seq_len(nrow(clinic_data))
  merged_data <- merge(missclassified_samples, clinic_data, by = "index")

  # Ensure the distance column is correctly integrated from the misclassified samples.
  merged_data$distance <- missclassified_samples$distance

  # --- Step 2: Find the "elbow" for cases ---------------------------------------
  # Filter distances for cases and sort them in descending order to find the threshold.
  distances_cases <- merged_data$distance[merged_data$diagnosis_type == "Case"]
  sorted_importance_cases <- sort(distances_cases, decreasing = TRUE)
  x_values_cases <- seq_along(sorted_importance_cases)
  y_values_cases <- sorted_importance_cases

  # Create a data frame for case distances to calculate the "elbow point".
  importance_elbow_cases <- data.frame(
    x = x_values_cases,
    y = y_values_cases
  )

  # Find the elbow point for cases using the helper function `find_curve_elbow`.
  elbow_cases <- find_curve_elbow(importance_elbow_cases)

  # Define the threshold for cases using the distance at the elbow point.
  threshold_case <- sorted_importance_cases[elbow_cases]

  # Select case samples with distances greater than or equal to the threshold.
  selected_case_samples <- merged_data[
    merged_data$diagnosis_type == "Case" &
      merged_data$distance >= threshold_case,
  ]

  # --- Step 3: Find the "elbow" for controls -----------------------------------
  # Filter distances for controls and sort them in descending order to find the threshold.
  distances_controls <- merged_data$distance[merged_data$diagnosis_type == "Control"]
  sorted_importance_controls <- sort(distances_controls, decreasing = TRUE)
  x_values_controls <- seq_along(sorted_importance_controls)
  y_values_controls <- sorted_importance_controls

  # Create a data frame for control distances to calculate the "elbow point".
  importance_elbow_controls <- data.frame(
    x = x_values_controls,
    y = y_values_controls
  )

  # Find the elbow point for controls using the helper function `find_curve_elbow`.
  elbow_controls <- find_curve_elbow(importance_elbow_controls)

  # Define the threshold for controls using the distance at the elbow point.
  threshold_control <- sorted_importance_controls[elbow_controls]

  # Select control samples with distances greater than or equal to the threshold.
  selected_control_samples <- merged_data[
    merged_data$diagnosis_type == "Control" &
      merged_data$distance >= threshold_control,
  ]

  # --- Step 4: Create a solution vector ----------------------------------------
  # Initialize a vector of zeros with the same length as the clinical data.
  # Mark selected samples (cases and controls) with "1" using their indices.
  solution_bootstrapping <- rep(0, nrow(clinic_data))
  solution_bootstrapping[selected_case_samples$index] <- 1
  solution_bootstrapping[selected_control_samples$index] <- 1

  # --- Step 5: Save results to files -------------------------------------------
  # Create directories for saving results if they do not exist.
  dir_path <- paste(saving_name, "bootstrapping_shuffle", sep = "/")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }

  # Save case samples to a file.
  case_data_path <- paste(saving_name, "case_data.rds", sep = "_")
  saveRDS(selected_case_samples, file = file.path(dir_path, case_data_path))

  # Save control samples to a file.
  control_data_path <- paste(saving_name, "control_data.rds", sep = "_")
  saveRDS(selected_control_samples, file = file.path(dir_path, control_data_path))

  # Save the solution vector to a file.
  solution_path <- paste(saving_name, "solution_bootstrapping.rds", sep = "_")
  saveRDS(solution_bootstrapping, file = file.path(dir_path, solution_path))

  # Save thresholds for controls and cases to separate files.
  threshold_control_path <- paste(saving_name, "threshold_control_bootstrapping.rds", sep = "_")
  saveRDS(threshold_control, file = file.path(dir_path, threshold_control_path))

  threshold_case_path <- paste(saving_name, "threshold_case_bootstrapping.rds", sep = "_")
  saveRDS(threshold_case, file = file.path(dir_path, threshold_case_path))

  # Save the merged data for reference.
  merge_data_path <- paste(saving_name, "merge_data_bootstrapping.rds", sep = "_")
  saveRDS(merged_data, file = file.path(dir_path, merge_data_path))
}

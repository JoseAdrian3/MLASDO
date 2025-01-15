#' Detect Anomalies
#'
#' Methodology based on the combination of a genetic algorithm and a Machine Learning technique to mutate the final
#' diagnosis of patients, detecting anomalies in them.
#'
#' These discrepancies may have their origin in the evolution of the subject itself, leading it from one group to
#' another in biological terms, or they may derive from human error in the labeling of the samples, among other causes.
#'
#' The identification and correction of these erroneous situations becomes essential to preserve the integrity
#' and accuracy of the data.
#'
#' @param num_trees Integer | Number of trees of the Random Forest model. Default value: 100.
#' @param mtry Integer | Number of predictors evaluated at each node of each tree. Default value: 225.
#' @param lambdas Numeric vector | Sequence of lambda values for Lasso regression. Default value: 10^seq(3, -3, length = 100).
#' @param omic_data_path String | Path to omic data. If not specified, sample data will be used.
#' @param clinic_data_path String | Path to clinical data. If not specified, sample data will be used.
#' @param omic_predictors_to_ignore Character vector | Variables to remove from the omic dataset.
#' @param clinic_predictors_to_ignore Character vector | Variables to remove from the clinical dataset.
#' @param id_column String | Identifier variable in both datasets. Default: "Trial" (for sample data).
#' @param active_predictors Character vector | Predictors for ratio study after genetic algorithm execution.
#' @param numeric_active_predictors Character vector | Numerical predictors among active predictors.
#' @param categoric_active_predictors Character vector | Categorical predictors among active predictors.
#' @param class_variable String | Target variable, binary with two possible values. Default: "Ca.Co.Last" (for sample data).
#' @param saving_name String | Name for saving the model and solution. Default: Current date as string.
#' @param n_cores Integer | Number of cores for parallelization. Default value: 6.
#' @param partition_percentage Numeric | Fraction of data split into training and test sets. Default value: 0.9 (90%).
#' @param diagnostic_change_probability Numeric | Probability of each gene in solutions to change. Default value: 0.1 (10%).
#' @param seed Integer | Random seed for reproducibility. Default value: 1234.
#' @param odds_default_values Character vector | Default values for categorical features, formatted as "feature:value".
#' @param weights Numeric vector | Weights assigned to transcriptomics predictors.
#' @param c Numeric vector | Cost parameter for SVM. Default value: c(1).
#' @param gamma Numeric vector | Gamma parameter for SVM. Default value: c(1).
#' @param kernel String | Kernel type for SVM. Default value: "lineal".
#' @param just_analysis Logical | If TRUE, only analysis is performed without modifications.
#' @param enrichment String | Type of enrichment analysis. Default value: "transcriptomic".
#' @param absolute_path String | Absolute path for file outputs.
#' @param omic String | Omic type to be analyzed.
#' @param domain_table Data frame | Domain table used for enrichment analysis.
#'
#' @export
#'
#' @examples
#' detect_anomalies(saving_name = "DefaultExecutionLasso")
#' detect_anomalies(saving_name = "DefaultExecutionRF")
#' detect_anomalies(saving_name = "QuickExecution", active_predictors = c("sex", "age", "mutation", "ethnicity"))
#' detect_anomalies(saving_name = "ExecutionWithOwnData", omic_data_path = "./myOmicData.tsv", clinic_data_path = "./myClinicData.tsv", id_column = "patient_id", class_variable = "diagnosis", active_predictors = c("sex", "age", "ethnicity"), numeric_active_predictors = c("age"), categoric_active_predictors = c("sex", "ethnicity"))

detect_anomalies <- function(
    num_trees = 500,
    mtry = 225,
    lambdas = 10^seq(3, -3, length = 100),
    omic_data_path = "",
    clinic_data_path = "",
    omic_predictors_to_ignore = NULL,
    clinic_predictors_to_ignore = NULL,
    id_column = NULL,
    active_predictors = NULL,
    numeric_active_predictors = NULL,
    categoric_active_predictors = NULL,
    class_variable = NULL,
    saving_name = "",
    n_cores = 6,
    partition_percentage = 0.9,
    diagnostic_change_probability = 0.1,
    seed = 1234,
    odds_default_values = character(0),
    weights = numeric(0),
    c = c(1),
    gamma = c(1),
    kernel = "lineal",
    just_analysis = FALSE,
    enrichment = "transcriptomic",
    absolute_path = "",
    omic = "",
    domain_table = NULL
){

  #### Seed for reproducibility ####
  set.seed(seed)

  #### Check and load required libraries ####
  required_packages <- c(
    "ranger", "dplyr", "plotly", "DT", "jsmodule", "ggpubr", "caret",
    "matrixStats", "e1071", "doParallel", "ggtext", "epitools",
    "clusterProfiler", "org.Hs.eg.db", "stringr", "kernlab",
    "pathviewr", "MetaboAnalystR", "umap", "dynamicTreeCut",
    "colorspace", "dendextend", "reactable", "UpSetR", "VennDiagram"
  )

  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("The package", pkg, "is not installed."))
    }
  }

  #' #### Data Reading ####

  if (omic_data_path == "") {
    omic_route <- system.file("data", "omicData.tsv", package = "MLASDO")
    class_variable <- "Ca.Co.Last"
    id_column <- "Trial"
  } else {
    omic_route <- omic_data_path
  }

  if (clinic_data_path == "") {
    clinic_route <- system.file("data", "clinicData.tsv", package = "MLASDO")
    class_variable <- "Ca.Co.Last"
    id_column <- "Trial"
  } else {
    clinic_route <- clinic_data_path
  }

  if (grepl("\\.tsv", omic_route)) {
    omic_data <- read.table(omic_route, header = TRUE, sep = "\t", row.names = 1)
  } else if (grepl("\\.csv", omic_route)) {
    omic_data <- read.csv(omic_route)
  } else if (grepl("\\.rds", omic_route)) {
    omic_data <- readRDS(omic_route)
  } else {
    return("The only valid formats for omic data are: .tsv and .csv")
  }

  if (grepl("\\.tsv", clinic_route)) {
    clinic_data <- read.table(clinic_route, header = TRUE, sep = "\t", row.names = 1)
  } else if (grepl("\\.csv", clinic_route)) {
    clinic_data <- read.csv(clinic_route)
  } else if (grepl("\\.rds", clinic_route)) {
    clinic_data <- readRDS(clinic_route)
  } else {
    return("The only valid formats for clinic data are: .tsv and .csv")
  }

  #' #### Relevel Clinic Data ####

  # Variable to store the features with a default value
  default_odds_predictors <- character(0)

  # Function to convert a dataframe's (df) feature (predictor) to a factor with the first value as the default value (default_value)
  convert_to_factor <- function(df, predictor, default_value) {
    df[[predictor]] <- factor(df[[predictor]], levels = c(default_value, setdiff(unique(df[[predictor]]), default_value)))
    return(df)
  }

  # From each string, separate by ':', with the left side being the predictor and the right side being the desired default value (future first level of the feature)
  for (entry in odds_default_values) {
    parts <- unlist(strsplit(entry, ":"))
    predictor <- parts[1]
    default_value <- parts[2]

    # Checking feature is not a numeric feature
    if (!predictor %in% numeric_active_predictors) {

      # Updating variable of features with default predictors
      default_odds_predictors <- c(default_odds_predictors, predictor)

      clinic_data <- convert_to_factor(clinic_data, predictor, default_value)

    } else {
      print("You have attempted to convert a numerical predictor to have a default value")
    }
  }

  #### Checking Parameters ####

  if (!is.null(omic_predictors_to_ignore)) {
    valid_columns_to_remove <- names(omic_data) %in% omic_predictors_to_ignore
    if (any(valid_columns_to_remove)) {
      omic_data <- omic_data[, !valid_columns_to_remove]
    }
  }

  if (!is.null(clinic_predictors_to_ignore)) {
    valid_columns_to_remove <- names(clinic_data) %in% clinic_predictors_to_ignore
    if (any(valid_columns_to_remove)) {
      clinic_data <- clinic_data[, !valid_columns_to_remove]
    }
  }

  if (is.null(active_predictors) & !is.null(numeric_active_predictors) & !is.null(categoric_active_predictors)) {
    active_predictors <- c(numeric_active_predictors, categoric_active_predictors)
  } else if (is.null(active_predictors)) {
    active_predictors <- names(clinic_data)
    active_predictors <- active_predictors[active_predictors != class_variable]
    active_predictors <- active_predictors[active_predictors != id_column]
  }

  # Checking if the patients in the three datasets are the same
  omic_ids <- sort(omic_data[[id_column]])
  clinic_ids <- sort(clinic_data[[id_column]])

  if (!identical(omic_ids, clinic_ids)) {
    return("The IDs of the patients don't match between omic and clinic datasets")
  }

  omic_class <- omic_data[c(id_column, class_variable)]
  clinic_class <- clinic_data[c(id_column, class_variable)]

  omic_class <- omic_class[order(omic_class[[id_column]]), ]
  clinic_class <- clinic_class[order(clinic_class[[id_column]]), ]

  rownames(omic_class) <- NULL
  rownames(clinic_class) <- NULL

  if (!identical(omic_class, clinic_class)) {
    return("The diagnostic of the patients don't match between omic and clinic datasets")
  }

  # Check if the selected variable as a predictor is included in the variables to be studied
  if (class_variable %in% active_predictors) {
    return("The selected variable for prediction cannot be included in the variables chosen for analysis")
  }

  if (id_column %in% active_predictors) {
    return("The variable that identifies the patients cannot be included in the variables chosen for analysis")
  }

  if (id_column %in% omic_predictors_to_ignore) {
    return("The variable that identifies the patients cannot be included in the variables chosen to ignore")
  }

  if (id_column %in% clinic_predictors_to_ignore) {
    return("The variable that identifies the patients cannot be included in the variables chosen to ignore")
  }

  # Check if the selected variable as the class variable for prediction is valid
  if (!class_variable %in% names(omic_data)) {
    return("The selected variable for prediction does not exist in the omic dataframe")
  }

  if (!class_variable %in% names(clinic_data)) {
    return("The selected variable for prediction does not exist in the clinic dataframe")
  }

  if (!id_column %in% names(clinic_data)) {
    return("The variable that identifies the patients does not exist in the clinic dataframe")
  }

  if (!id_column %in% names(omic_data)) {
    return("The variable that identifies the patients does not exist in the omic dataframe")
  }

  # Check if the selected variable as the class variable for prediction is binary
  if (length(unique(omic_data[[class_variable]])) > 2) {
    return("The selected variable for prediction is not binary in the omic dataset")
  }

  if (length(unique(clinic_data[[class_variable]])) > 2) {
    return("The selected variable for prediction is not binary")
  }

  # Checking if the class variable is well coded
  if (all(!(omic_data[[class_variable]] %in% c(0, 1))) & all(!(omic_data[[class_variable]] %in% c("Case", "Control")))) {
    return("The class variable is not well coded in the omic dataset")
  }

  if (all(!(clinic_data[[class_variable]] %in% c(0, 1))) & all(!(clinic_data[[class_variable]] %in% c("Case", "Control")))) {
    return("The class variable is not well coded")
  }

  # Checking if the coding of the class variable in the omic dataset is correct
  if (0 %in% unique(omic_data[[class_variable]])) {
    cat("The class variable is coded with 0, those values are changed to Control.\n")
    omic_data[[class_variable]][omic_data[[class_variable]] == 0] <- "Control"
    clinic_data[[class_variable]][clinic_data[[class_variable]] == 0] <- "Control"
  }

  if (1 %in% unique(omic_data[[class_variable]])) {
    cat("The class variable is coded with 1, those values are changed to Case.\n")
    omic_data[[class_variable]][omic_data[[class_variable]] == 1] <- "Case"
    clinic_data[[class_variable]][clinic_data[[class_variable]] == 1] <- "Case"
  }

  if (0 %in% unique(clinic_data[[class_variable]])) {
    cat("The class variable is coded with 0, those values are changed to Control.\n")
    omic_data[[class_variable]][omic_data[[class_variable]] == 0] <- "Control"
    clinic_data[[class_variable]][clinic_data[[class_variable]] == 0] <- "Control"
  }

  if (1 %in% unique(clinic_data[[class_variable]])) {
    cat("The class variable is coded with 1, those values are changed to Case.\n")
    omic_data[[class_variable]][omic_data[[class_variable]] == 1] <- "Case"
    clinic_data[[class_variable]][clinic_data[[class_variable]] == 1] <- "Case"
  }

  # Finally, check if the variables selected for statistical analysis are valid
  invalid_predictors <- active_predictors[!(active_predictors %in% names(clinic_data))]

  if (length(invalid_predictors) > 0) {
    cat("The following variables selected for analysis do not exist in the dataframe:\n")
    cat(paste(invalid_predictors, collapse = ", "))
    cat("\n")
    active_predictors <- active_predictors[!(active_predictors %in% invalid_predictors)]
  }

  if (is.null(categoric_active_predictors) & !is.null(numeric_active_predictors)) {
    cat("You have only indicated which active predictors are numeric, all other active predictors will be treated as categorical.\n")
    rest_of_predictors <- active_predictors[!(active_predictors %in% numeric_active_predictors)]
    categoric_active_predictors <- rest_of_predictors
    invalid_predictors <- numeric_active_predictors[!(numeric_active_predictors %in% names(clinic_data))]
    numeric_active_predictors <- numeric_active_predictors[!(numeric_active_predictors %in% invalid_predictors)]
    cat("Active numerical predictors: ")
    cat(paste(numeric_active_predictors, collapse = ", "))
    cat("\n")
    cat("Active categoric predictors: ")
    cat(paste(categoric_active_predictors, collapse = ", "))
    cat("\n")
  }

  if (!is.null(categoric_active_predictors) & is.null(numeric_active_predictors)) {
    cat("You have only indicated which active predictors are categorical, all other active predictors will be treated as numeric.\n")
    rest_of_predictors <- active_predictors[!(active_predictors %in% categoric_active_predictors)]
    numeric_active_predictors <- rest_of_predictors
    invalid_predictors <- categoric_active_predictors[!(categoric_active_predictors %in% names(clinic_data))]
    categoric_active_predictors <- categoric_active_predictors[!(categoric_active_predictors %in% invalid_predictors)]
    cat("Active numerical predictors: ")
    cat(paste(numeric_active_predictors, collapse = ", "))
    cat("\n")
    cat("Active categoric predictors: ")
    cat(paste(categoric_active_predictors, collapse = ", "))
    cat("\n")
  }

  if (!is.null(categoric_active_predictors) & !is.null(numeric_active_predictors)) {
    repeated_predictors <- intersect(numeric_active_predictors, categoric_active_predictors)
    if (length(repeated_predictors) > 0) {
      cat("The following variables selected for analysis are listed as categorical and numerical variables:\n")
      cat(paste(repeated_predictors, collapse = ", "))
      return("Please remove repeated variables.")
    }
  }

  # If the name is empty, we need to establish a default one
  if (saving_name == "") {
    saving_name <- format(Sys.time(), "%d-%m-%H-%M")
  }

  #### SVM METHOD ####
  # Check if analysis mode is enabled (if not, execute the SVM algorithm and prepare data)
  if (!just_analysis) {

    # Create a directory to save SVM-related data
    dir.create(paste(saving_name, "svm_data", sep = "/"))

    cat("\nExecuting the SVM algorithm.\n")

    # Execute the SVM algorithm using the provided parameters and save results
    MLASDO::execute_svm(
      saving_name = saving_name,
      omic_data = omic_data,
      class_variable = class_variable,
      id_column = id_column,
      n_cores = n_cores,
      seed = seed,
      c = c_value,       # `C` renamed to `c_value` for clarity
      kernel = kernel,
      gamma = gamma
    )
  }

  # Define the path for the SVM data directory
  svm_data_dir_path <- paste(saving_name, "svm_data", saving_name, sep = "/")

  # Reading the primary SVM data
  svm_data_path <- paste(svm_data_dir_path, "svm_data.rds", sep = "_")
  svm_data <- readRDS(svm_data_path)

  # Reading the complete SVM data (all results)
  svm_data_all_path <- paste(svm_data_dir_path, "svm_data_all.rds", sep = "_")
  svm_data_all <- readRDS(svm_data_all_path)

  ### OUTLIERS SAMPLES ###
  # Create a directory for bootstrapping and shuffle results
  dir.create(paste(saving_name, "bootstrapping_shuffle", sep = "/"))

  # Perform bootstrapping and shuffling to identify outlier samples
  MLASDO::bootstrapping_shuffle(
    clinic_data = clinic_data, # `clinicData` renamed to `clinic_data`
    svm_data = svm_data,
    seed = seed,
    saving_name = saving_name
  )

  # Define the path for the bootstrapping results directory
  bootstrapping_dir_path <- paste(saving_name, "bootstrapping_shuffle", saving_name, sep = "/")

  # Load the solution generated by bootstrapping
  bootstrapping_solution_path <- paste(bootstrapping_dir_path, "solution_bootstrapping.rds", sep = "_")
  bootstrapping_solution <- readRDS(bootstrapping_solution_path)

  #### STATISTICAL TESTS ####

  # Create a directory to save statistical analysis results for the SVM method
  dir.create(paste(saving_name, "statistics_svm", sep = "/"))

  cat("Performing statistical analysis.\n")

  # Perform statistical calculations for the SVM method
  MLASDO::statistical_calculations(
    method = "svm",
    omic = omic_data,
    active_predictors = active_predictors,
    class_variable = class_variable,
    clinic_data = clinic_data,
    id_column = id_column,
    solution = bootstrapping_solution,
    categoric_active_predictors = categoric_active_predictors,
    numeric_active_predictors = numeric_active_predictors,
    default_odds_predictors = default_odds_predictors,
    saving_name = saving_name,
    dir_path = "statistics_svm",
    seed = seed
  )

  #### AFTER AND BEFORE MODELS ####

  # Directory for saving SVM model results
  name_dir <- "models_svm"
  dir.create(paste(saving_name, name_dir, sep = "/"))

  # BEFORE processing: Preprocessing data before making changes
  before_data <- MLASDO::preprocess_data(omic_data, class_variable, partition_percentage, id_column, seed)

  # Path for saving model results
  dir_path <- paste(saving_name, name_dir, saving_name, sep = "/")

  # Indices for training data
  train_indices <- before_data$indices_train

  # Processing data after applying changes from bootstrapping
  changed_data_svm <- MLASDO::preprocess_changed_data(
    omic = omic_data,
    clinic_data = clinic_data,
    solution = bootstrapping_solution,
    active_predictors = active_predictors,
    class_variable = class_variable,
    id_column = id_column,
    svm_data = svm_data,
    seed = seed
  )

  # Extracting changed omic data and splitting into training and test sets
  changed_omic_data <- changed_data_svm$changed_omic_data
  train_data <- changed_omic_data[train_indices, ]
  test_data <- changed_omic_data[-train_indices, ]

  # Removing class and ID columns from training and test sets
  omic_train <- train_data[, !(colnames(train_data) %in% c(class_variable, id_column))]
  omic_test <- test_data[, !(colnames(test_data) %in% c(class_variable, id_column))]

  # Creating a list to store processed data after changes
  after_data <- list(
    omic_train = omic_train,
    omic_test = omic_test,
    omic_train_diagnosis = changed_omic_data[train_indices, class_variable],
    omic_test_diagnosis = changed_omic_data[-train_indices, class_variable]
  )

  # Check if analysis is enabled (skip this if `just_analysis` is TRUE)
  if (!just_analysis) {

    #### BEFORE MODELS ####
    cat("RF BEFORE\n\n")
    svm_model_rf <- MLASDO::train_rf(
      before_data$omic_train,
      before_data$omic_train_diagnosis,
      before_data$omic_test,
      before_data$omic_test_diagnosis,
      num_trees,
      mtry,
      seed,
      dir_path,
      "before",
      n_cores
    )

    cat("LASSO BEFORE\n\n")
    svm_model_lasso <- MLASDO::train_lasso(
      before_data$omic_train,
      before_data$omic_train_diagnosis,
      before_data$omic_test,
      before_data$omic_test_diagnosis,
      lambdas,
      seed,
      dir_path,
      "before",
      n_cores
    )

    cat("SVM LINEAR BEFORE\n\n")
    svm_model_svm_linear <- MLASDO::train_svm_linear(
      before_data$omic_train,
      before_data$omic_train_diagnosis,
      before_data$omic_test,
      before_data$omic_test_diagnosis,
      c_value,
      seed,
      dir_path,
      "before",
      n_cores
    )

    cat("SVM RADIAL BEFORE\n\n")
    svm_model_svm_radial <- MLASDO::train_svm_radial(
      before_data$omic_train,
      before_data$omic_train_diagnosis,
      before_data$omic_test,
      before_data$omic_test_diagnosis,
      c_value,
      seed,
      dir_path,
      "before",
      n_cores
    )

    #### AFTER MODELS ####
    cat("RF AFTER\n\n")
    svm_model_rf <- MLASDO::train_rf(
      after_data$omic_train,
      after_data$omic_train_diagnosis,
      after_data$omic_test,
      after_data$omic_test_diagnosis,
      num_trees,
      mtry,
      seed,
      dir_path,
      "after",
      n_cores
    )

    cat("LASSO AFTER\n\n")
    svm_model_lasso <- MLASDO::train_lasso(
      after_data$omic_train,
      after_data$omic_train_diagnosis,
      after_data$omic_test,
      after_data$omic_test_diagnosis,
      lambdas,
      seed,
      dir_path,
      "after",
      n_cores
    )

    cat("SVM LINEAR AFTER\n\n")
    svm_model_svm_linear <- MLASDO::train_svm_linear(
      after_data$omic_train,
      after_data$omic_train_diagnosis,
      after_data$omic_test,
      after_data$omic_test_diagnosis,
      c_value,
      seed,
      dir_path,
      "after",
      n_cores
    )

    cat("SVM RADIAL AFTER\n\n")
    svm_model_svm_radial <- MLASDO::train_svm_radial(
      after_data$omic_train,
      after_data$omic_train_diagnosis,
      after_data$omic_test,
      after_data$omic_test_diagnosis,
      c_value,
      seed,
      dir_path,
      "after",
      n_cores
    )
  }


  cat("Compiling Markdown file.\n")
  MLASDO::compileMarkdown(
    classVariable = class_variable,
    savingName = saving_name,
    originalDiagnosis = omic_data[[class_variable]],
    clinicDataSVM = changed_data_svm$transition_clinic_data,
    categoricActivePredictors = categoric_active_predictors,
    numericActivePredictors = numeric_active_predictors,
    idColumn = id_column,
    numTrees = num_trees,
    mtrysvm = mtry,
    diagnosticChangeProbability = diagnostic_change_probability,
    nCores = n_cores,
    seed = seed,
    lambdas = lambdas,
    C = c,
    gamma = gamma,
    kernel = kernel,
    partitionPercentage = partition_percentage,
    enrichment = enrichment,
    absolutePath = absolute_path,
    omicData = omic_data,
    omic = omic,
    domainTable = domain_table
  )
}

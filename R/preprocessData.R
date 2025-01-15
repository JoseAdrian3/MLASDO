#' Preprocess Data for Training and Testing
#'
#' This function partitions an omics dataset into training and testing subsets,
#' removing the class variable and a specified ID column from the feature set.
#'
#' @param omic A data frame containing the omics dataset.
#' @param class_variable A string representing the name of the class variable column.
#' @param partition_percentage A numeric value specifying the proportion of the dataset to allocate to training (between 0 and 1).
#' @param id_column A string specifying the name of the ID column to be removed from the dataset.
#' @param seed An integer used to set the random seed for reproducibility.
#'
#' @return A list containing the following elements:
#'   \item{omic_train}{A data frame with training features.}
#'   \item{omic_train_diagnosis}{A vector of class labels for the training set.}
#'   \item{omic_test}{A data frame with testing features.}
#'   \item{omic_test_diagnosis}{A vector of class labels for the testing set.}
#'   \item{indices_train}{Indices of the rows used for training.}
#'
#' @export
preprocess_data <- function(omic, class_variable, partition_percentage, id_column, seed) {

  # Set random seed for reproducibility
  set.seed(seed)

  # Generate indices for partitioning the dataset into training and testing subsets
  subset_train <- createDataPartition(y = omic[[class_variable]],
                                      p = partition_percentage,
                                      list = FALSE)[, 1]

  # Store the original ID column (if needed later)
  omic_id <- omic$id

  # Remove the specified ID column from the dataset
  omic[[id_column]] <- NULL

  # Create the training dataset
  omic_train <- omic[subset_train, ]

  # Extract and remove the class variable from the training dataset
  omic_train_diagnosis <- omic_train[[class_variable]]
  omic_train[[class_variable]] <- NULL

  # Create the testing dataset
  omic_test <- omic[-subset_train, ]

  # Extract and remove the class variable from the testing dataset
  omic_test_diagnosis <- omic_test[[class_variable]]
  omic_test[[class_variable]] <- NULL

  # Return a list with processed data and indices
  list(
    omic_train = omic_train,
    omic_train_diagnosis = omic_train_diagnosis,
    omic_test = omic_test,
    omic_test_diagnosis = omic_test_diagnosis,
    indices_train = subset_train
  )
}

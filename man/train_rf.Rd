% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trainRF.R
\name{train_rf}
\alias{train_rf}
\title{Train a Random Forest Model}
\usage{
train_rf(
  omic_train,
  omic_train_diagnosis,
  omic_test,
  omic_test_diagnosis,
  num_trees,
  mtry,
  seed,
  dir_path,
  file_suffix,
  n_cores
)
}
\arguments{
\item{omic_train}{Transcriptomic train dataframe with features (independent variables, x).}

\item{omic_train_diagnosis}{Vector of training data labels (dependent variables, y).}

\item{omic_test}{Transcriptomic test dataframe with features (independent variables, x).}

\item{omic_test_diagnosis}{Vector of test data labels (diagnoses).}

\item{num_trees}{Integer indicating the number of trees to grow in the random forest.}

\item{mtry}{Integer specifying the number of variables randomly sampled as candidates at each split.}

\item{seed}{Seed for reproducibility.}

\item{dir_path}{Directory path where output files will be saved.}

\item{file_suffix}{Suffix to append to saved files.}

\item{n_cores}{Number of CPU cores to use for parallel processing.}
}
\value{
A trained Random Forest model.
}
\description{
This function trains a Random Forest (RF) model using cross-validation and calculates performance metrics.
The trained model, confusion matrix, and confidence intervals are saved to specified files.
}

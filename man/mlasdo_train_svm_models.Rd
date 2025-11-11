\name{mlasdo_train_svm_models}
\alias{mlasdo_train_svm_models}
\title{Train SVM models on the feature-selected MLASDO dataset}
\usage{
mlasdo_train_svm_models(
  omic_data,
  id_column = "participant_id",
  target_column = "diagnosis_type",
  output_dir,
  kernels = c("linear", "radial", "polynomial"),
  linear_grid = NULL,
  radial_grid = NULL,
  polynomial_grid = NULL,
  linear_cost_values = 10^seq(-6, 6, length.out = 13),
  radial_cost_values = 10^seq(-4, 4, length.out = 9),
  radial_sigma_values = 10^seq(-6, 1, length.out = 8),
  polynomial_cost_values = 10^seq(-4, 4, length.out = 9),
  polynomial_degree_values = 2:4,
  polynomial_scale_values = 10^seq(-6, 0, length.out = 7),
  num_folds = 10,
  folds_seed = 2025,
  n_cores = max(1L, parallel::detectCores() - 1L),
  sampling = "down",
  positive_class = NULL,
  negative_class = NULL
)
}
\arguments{
\item{omic_data}{Data frame or path to an \code{.rds} file with the trimmed
omic dataset. Must contain \code{id_column} and \code{target_column}.}

\item{id_column}{Identifier column (defaults to \code{"participant_id"}).}

\item{target_column}{Target variable column (defaults to
\code{"diagnosis_type"}).}

\item{output_dir}{Directory where all kernel-specific folders and the summary
will be written.}

\item{kernels}{Character vector with any subset of \code{"linear"},
\code{"radial"}, and \code{"polynomial"}. Defaults to all three.}

\item{linear_grid, radial_grid, polynomial_grid}{Optional \code{expand.grid()}
data frames with the hyperparameter combinations evaluated for each kernel.
\code{NULL} triggers the default grids.}

\item{linear_cost_values}{Candidate \code{C} values for the linear kernel when
\code{linear_grid} is \code{NULL}.}

\item{radial_cost_values}{Candidate \code{C} values for the radial kernel when
\code{radial_grid} is \code{NULL}.}

\item{radial_sigma_values}{Candidate \code{sigma} (gamma) values for the radial
kernel when \code{radial_grid} is \code{NULL}.}

\item{polynomial_cost_values}{Candidate \code{C} values for the polynomial kernel
when \code{polynomial_grid} is \code{NULL}.}

\item{polynomial_degree_values}{Candidate polynomial degrees when
\code{polynomial_grid} is \code{NULL}.}

\item{polynomial_scale_values}{Candidate \code{scale} (\code{gamma}) values for
the polynomial kernel when \code{polynomial_grid} is \code{NULL}.}

\item{num_folds}{Number of folds used for cross-validation (defaults to 10).}

\item{folds_seed}{Seed applied before generating the shared folds.}

\item{n_cores}{Number of worker processes used by \code{caret::train}.}

\item{sampling}{Strategy passed to \code{caret::trainControl(sampling = ...)}.}

\item{positive_class}{Label treated as the positive class (used to order the
factor levels and for confusion matrices).}

\item{negative_class}{Label treated as the negative class.}
}
\value{
A list with per-kernel results, the summary table, and the best
kernel metadata. Side effects include all the saved RDS/CSV files.
}
\description{
Executes the post-feature-selection modelling stage: the same 10-fold split
is reused across linear, radial and polynomial kernels, the best
hyperparameters per kernel are learned via cross-validation on the full
dataset, and each fold is then retrained with those hyperparameters to obtain
predictions, distances to the separating hyperplane and Balanced Accuracy
metrics. All artifacts are stored under \code{output_dir}, including a
summary with the kernel that achieved the highest average Balanced Accuracy.
}
\examples{
\dontrun{
mlasdo_train_svm_models(
  omic_data = "data/mlasdo_runs/feature_selection_only/feature_selection/omic_selected_features.rds",
  output_dir = "data/svm/feature_importance_shared_cv"
)
}
}

# Machine Learning based Anomalous Sample Detection on Omics (MLASDO)

MLASDO (Machine Learning based Anomalous Sample Detection on Omics) is an R tool that uncovers patients whose omic profiles disagree with their clinical labels and explains why. The package aligns omic and clinical datasets, prioritises the most informative molecular features, trains Support Vector Machines (SVM) to detect outliers, and summarises the findings in reproducible reports that can be shared with clinical experts.

## Workflow Overview

1. **Data preparation** – harmonise `omic_data` and `clinical_data` for the same set of IDs, validate targets, and preserve required covariates.
2. **Feature selection (optional)** – fit a tuned random forest, measure feature importance, and keep only the molecular signals that best separate the clinical groups.
3. **SVM-based anomaly detection** – train multiple kernels (linear, radial, polynomial) with user-defined grids; detect anomalous subjects based on their distances to each SVM hyperplane.
4. **Anomaly validation** – relabel candidates and test whether discrimination improves using both linear and non-linear learners.
5. **Reporting** – generate HTML/Markdown summaries, feature importance tables, anomaly evidence, and domain expert tables for each patient.

## Installation

```r
install.packages("devtools")

library(devtools)
install_github("JoseAdrian3/MLASDO")

library(MLASDO)
```

MLASDO depends on standard tidyverse tooling plus `caret`, `ranger`, `plotly`, `DT`, `e1071`, `matrixStats`, `doParallel`, `pathviewr`, and related utilities. These packages are pulled automatically through DESCRIPTION.

## Input Requirements

- `omic_data`: either an `.rds` file or an in-memory data frame whose rows represent subjects and columns are omic measurements. You can assume an object named `omic_data` exists when following the examples.
- `clinical_data`: path or data frame with the clinical covariates (e.g., demographics, scores). The examples refer to a `clinical_data` object with the same subjects as `omic_data`.
- `id_column`: column shared by both tables that uniquely identifies each subject.
- `target_column`: factor or character column describing the clinical class (e.g., `"AHC"` vs `"HC"`).

## Quick Start

```r
library(MLASDO)

run_mlasdo(
  omic_data = omic_data,
  clinical_data = clinical_data,
  id_column = "patient_id",
  target_column = "diagnosis",
  output_dir = "mlasdo_runs/basic_demo"
)
```

This single call aligns the provided datasets, runs the default feature selection and SVM configuration, detects anomalies, and stores every artifact under `mlasdo_runs/basic_demo`.

## Example Workflows

### 1. Minimal execution with file paths

```r
run_mlasdo(
  omic_data = "data/omic_data.rds",
  clinical_data = "data/clinical_data.rds",
  id_column = "patient_id",
  target_column = "diagnosis",
  output_dir = "mlasdo_runs/from_files"
)
```

### 2. Reusing an existing feature-selection run (no new feature selection)

If you already executed feature selection and only want to re-run SVMs/anomaly validation, point MLASDO to the saved directory:

```r
run_mlasdo(
  omic_data = omic_data,
  clinical_data = clinical_data,
  id_column = "patient_id",
  target_column = "diagnosis",
  output_dir = "mlasdo_runs/reuse_fs",
  feature_selection_exec = FALSE,
  feature_selection_existing_dir = "mlasdo_runs/baseline_fs/feature_selection",
  svm_exec = TRUE,
  anomalies_exec = TRUE
)
```

### 3. Full run

This snippet mirrors the full testing pipeline stored under `scripts/tests/tests.R` and shows explicit diagnostics before calling MLASDO:

```r
run <- run_mlasdo(
  omic_data                     = omic_data_df,
  clinical_data                 = clinical_data_df,
  id_column                     = "participant_id",
  target_column                 = "diagnosis_type",
  target_positive_class         = "Case",
  target_negative_class         = "Control",
  output_dir                    = "mlasdo_runs/run_ppmi_sc_pc_nm",
  always_split_variables        = c("age_at_baseline", "sex"),
  rf_mtry                       = unique(round(seq(2, floor(sqrt(ncol(omic_data_df))), length.out = 3))),
  rf_splitrule                  = c("gini", "extratrees", "hellinger"),
  rf_min_node_size              = c(1, 5, 10),
  rf_train_method               = "repeatedcv",
  rf_train_number               = 10,
  rf_train_repeats              = 5,
  rf_num_trees                  = 500,
  rf_importance                 = "permutation",
  seed                          = 123,
  svm_kernels                   = c("linear", "radial", "polynomial"),
  svm_linear_cost_values        = 10^seq(-4, 3, length = 8),
  svm_radial_cost_values        = 10^seq(-4, 3, length = 8),
  svm_radial_sigma_values       = 10^seq(-7, 0, length = 8),
  svm_polynomial_cost_values    = 10^seq(-4, 3, length = 8),
  svm_polynomial_scale_values   = 10^seq(-7, 0, length = 8),
  svm_polynomial_degree_values  = 2:4,
  svm_num_folds                 = 10,
  svm_sampling                  = "down",
  anomalies_method              = "sd",
  anomalies_sd_multiplier       = 2
)
```

### 4. Full run regenerating only anomaly outputs

When feature selection and SVM models are already available, you can replay only the anomaly detection (and downstream clinical tests) with a minimal call:

```r
run_mlasdo(
  omic_data                       = omic_data_df,
  clinical_data                   = clinical_data_df,
  id_column                       = "participant_id",
  target_column                   = "diagnosis_type",
  target_positive_class           = "Case",
  target_negative_class           = "Control",
  output_dir                      = "mlasdo_runs/run_ppmi_sc_pc_nm",
  feature_selection_exec          = FALSE,
  feature_selection_existing_dir  = "mlasdo_runs/run_ppmi_sc_pc_nm/feature_selection",
  svm_exec                        = FALSE,
  svm_existing_dir                = "mlasdo_runs/run_ppmi_sc_pc_nm/svm_models",
  anomalies_exec                  = TRUE,
  anomalies_method                = "sd",
  anomalies_sd_multiplier         = 2
)
```

## Outputs

Each run creates timestamped folders under `output_dir` including:

- Preprocessed omic/clinical datasets with aligned IDs.
- Feature-selection artifacts (trained RF model, sorted importance table, trimmed omic matrix).
- SVM grids, distance matrices, and best-kernel summaries.
- Anomaly tables, relabelling evidence, and clinical interpretation statistics.
- HTML/Markdown reports combining all plots and tables for group and individual review.

## Credits

- **Project Leaders**: Juan Antonio Botia Blaya · Alicia Gómez Pascual  
- **Main Developer**: José Adrián Pardo Pérez
- **Preprint**: arxiv.org/abs/2507.03656 
## Contact

`joseadrian.pardop@um.es`

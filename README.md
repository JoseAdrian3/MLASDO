# MLASDO (a software tool to detect and explain omic-based inconsistencies within patient groups in longitudinal cohorts)

## MOTIVATION

In most omics analysis datasets (data with different biological aspects), a relevant problem arises: the presence of samples that, in reality, should be classified as belonging to the alternative group to which they have been initially assigned.

These discrepancies may have their origin in the evolution of the subject itself, leading it from one group to another in biological terms, or they may derive from human error in the labelling of the samples, among other causes.

The identification and correction of these erroneous situations becomes essential to preserve the integrity and accuracy of the data. This problem makes it crucial to develop efficient tools capable of detecting and correcting these classification errors in omics datasets, thus ensuring the reliability of the conclusions derived from biological studies.

## What approach do we follow to solve it

We follow five steps:

### Step 1. Data preparation

Preparation of two datasets with same patients: 

- Dataset with any omics content (transcriptomics, metabolomics, proteomics, methylomics...) of each patient.

- Clinical data dataset, this will contain different clinical data for each patient.

### Step 2. Detection of anomalous samples

The detection of anomalous samples is based on Support Vector Machines (SVM). To optimize the model, MLASDO receive parameters such as ```c```, ```gamma``` and ```kernel```. The detection of anomalous samples is based on their distance to the SVM classifier's hyperplane.

### Step 4. Anomalous samples validation

MLASDO verifies anomalous samples by assessing whether relabeling them improves class separation. It does so by training linear (SVM Linear, Lasso) and non-linear (SVM Radial, Random Forest) models, assigning new labels as determined by MLASDO (e.g., AHC to PD, APD to HC).s

### Step 5. Generating reports to evidence the detection of abnormal cases

Finally, MLASDO generates reports at both group and individual levels for the anomalous samples.

### R packages required
- devtools
- ranger
- dplyr
- plotly
- DT
- jsmodule
- caret
- matrixStats
- e1071
- doParallel
- pathviewr

### Installation
```
    install.packages("devtools")
 
    library(devtools)

    install_github("JoseAdrian3/MLASDO")

    library(MLASDO)
```

### Examples of use of such software

Basic MLASDO execution with all default values and dataset.
```
MLASDO::detectAnomalies(
        saving_name = "example_folder_name", 
        class_variable = "diagnosis", 
        id_column = "id")
```

Execution specifying a dataset.
```
MLASDO::detectAnomalies(
  saving_name = "example_folder_name", 
  class_variable = "diagnosis", 
  id_column = "id")
  omic_data_path = "./omic_data.rds",
  clinic_data_path = "./clinic_data.rds",
)
```

MLASDO executing only the analysis after the SVM.
```
MLASDO::detectAnomalies(
  saving_name = "example_folder_name", 
  class_variable = "diagnosis", 
  id_column = "id")
  omic_data_path = "./omic_data.rds",
  clinic_data_path = "./clinic_data.rds",
  just_analysis = TRUE
)
```

MLASDO indicating the numerical covariates of the clinic dataset.
```
example_numeric_activate_predictors <- clinic_df[, sapply(clinic_df, is.numeric)]

MLASDO::detectAnomalies(
  saving_name = "example_folder_name", 
  class_variable = "diagnosis", 
  id_column = "id")
  omic_data_path = "./omic_data.rds",
  clinic_data_path = clinic_df,
  just_analysis = TRUE,
  numeric_active_predictors = example_numeric_activate_predictors,
)
```

Complete MLASDO execution.
```
example_numeric_activate_predictors <- clinic_df[, sapply(clinic_df, is.numeric)]

MLASDO::detectAnomalies(
  saving_name = "example_folder_name", 
  class_variable = "diagnosis", 
  id_column = "id")
  omic_data_path = "./omic_data.rds",
  clinic_data_path = clinic_df,
  just_analysis = TRUE,
  numeric_active_predictors = example_numeric_activate_predictors,
  partition_percentage = 0.9,
  c = c(1000, 100, 10, 1, 0.1, 0.01, 0.001),
  gamma = c(0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001),
  kernel = c("linear", "radial"),
  n_cores = detectCores() - 2,
  seed = 1234,
  domain_table = "./domain_table.rds"
)
```

## Credits
**Project Leaders:** Juan Antonio Botia Blaya https://github.com/juanbot and Alicia Gómez Pascual https://github.com/aliciagp

**Main developer:** José Adrián Pardo Pérez https://github.com/JoseAdrian3

## Contact
**email address:** joseadrian.pardop@um.es

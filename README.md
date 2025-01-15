# MLASDO (a software tool to detect and explain omic-based inconsistencies within patient groups in longitudinal cohorts)

## What problem are we seeking to solve

In most omics analysis datasets (data with different biological aspects), a relevant problem arises: the presence of samples that, in reality, should be classified as belonging to the alternative group to which they have been initially assigned.

These discrepancies may have their origin in the evolution of the subject itself, leading it from one group to another in biological terms, or they may derive from human error in the labelling of the samples, among other causes.

The identification and correction of these erroneous situations becomes essential to preserve the integrity and accuracy of the data. This problem makes it crucial to develop efficient tools capable of detecting and correcting these classification errors in omics datasets, thus ensuring the reliability of the conclusions derived from biological studies.

## What approach do we follow to solve it

We follow three steps:

### Step 1. Data preparation

This software uses two datasets: 

- Dataset with any omics content (transcriptomics, metabolomics, proteomics, methylomics...) of each patient. This dataset will be used to train the machine learning algorithm in charge of classifying the patients.

- Clinical data dataset, this will contain different clinical data for each patient. This dataset will be used to create a report based on ratios and statistical tests, with this report providing evidence to support that the diagnoses identified as abnormal are indeed abnormal.

Therefore, the first step will be to check that both datasets are composed of the same patients, ensuring consistency between the detection of the anomalous cases and their subsequent evidencing.

After this, it is necessary to make sure that clinical data are not present in the omics dataset, thus avoiding the influence of such data in the detection of the anomalous cases. 

### Step 2. Detection of anomalous PD and HC

The detection of anomalous samples is based on Support Vector Machines (SVM) (Cortes and Vapnik, 1995).

The first step performed by MLASDO is to optimize an SVM to distinguish between the groups as indicted by the target feature, (e.g., healthy controls and cases). To optimize the model, MLASDO receive parameters such as ```c```, ```gamma``` and ```kernel```. The model uses a 10-fold cross-validation approach. MLASDO allows the use of a downsampling approach to address potential class imbalance.

The detection of anomalous samples is based on their distance to the SVM classifier's hyperplane. For the linear kernel, distances to the decision boundary are calculated based on:

\[
D(x) = \sum_{k=1}^p \alpha_k K(x_k, x) + b
\]

\[
|d(x)| = \frac{|D(x)|}{\|w\|}
\]

For the radial kernel, distances are obtained with the expression:

\[
K(x, x') = \exp\left(-\gamma \sum_{i=1}^p (x_{ij} - x'_{ij})^2\right)
\]

Anomalous samples are defined as those misclassified whose distance to the hyperplane is over a threshold. To determine an appropriate threshold, the MLASDO approach employs an elbow rule.

### Step 4. Anomalous samples validation

In order to verify the presence of the anomalous samples, MLASDO verifies whether the changes the labels of the anomalous samples (e.g. AHC changed by case) would improve the separation of the classes. 

To this end, MLASDO trained different linear (SVM Linear and Lasso) and non-linear (SVM Radial and Random Forest) models but this time changing the label of the anomalous samples to the label assigned by MLASDO (i.e. AHC are changed to PD and APD are changed to HC).

### Step 3. Generating reports to evidence the detection of abnormal cases



## What software we have developed for this purpose

To solve this problem, the MLASDO software (a software tool to detect and explain omics-based inconsistencies within patient groups in longitudinal cohorts), which is an R package, has been implemented.

In this software it is possible to configure
- The different parameters of the SVM.
- The different parameters of the models to validate the anomalous samples (Random Forest, Lasso, SVM parameters...).
- The variable on which you want to detect anomalous cases.
- The clinical data on which to perform the analysis.
- The omic data on which to perform the analysis. 
- The percentage split between training and test set.
- If the user wants to perform the entire MLASDO exec or just the analysis

In addition, with this software, it is possible to perform the detection and subsequent analysis on your own datasets.

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
**Project Leader:** Juan Antonio Botia Blaya https://github.com/juanbot

**Main developer:** José Adrián Pardo Pérez https://github.com/JoseAdrian3

## Contact
**email address:** joseadrian.pardop@um.es

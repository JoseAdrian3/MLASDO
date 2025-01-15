# MLASDO (Machine learning based Anomalous Sample Detection on Omics)

## What problem are we seeking to solve

In most omics analysis datasets (data with different biological aspects), a relevant problem arises: the presence of samples that, in reality, should be classified as belonging to the alternative group to which they have been initially assigned.

These discrepancies may have their origin in the evolution of the subject itself, leading it from one group to another in biological terms, or they may derive from human error in the labelling of the samples, among other causes.

The identification and correction of these erroneous situations becomes essential to preserve the integrity and accuracy of the data. This problem makes it crucial to develop efficient tools capable of detecting and correcting these classification errors in omics datasets, thus ensuring the reliability of the conclusions derived from biological studies.

## What approach do we follow to solve it

We follow three steps:

### Step 1. Data preparation

This software uses two datasets: 

- Dataset with any omics content (transcriptomics, metabolomics, proteomics, methylomics...) of each patient. This dataset will be used to train the machine learning algorithm in charge of classifying the patients.

- Clinical data dataset, this will contain different clinical data for each patient. This dataset will be used to create a report based on ratios, with this report providing evidence to support that the diagnoses identified as abnormal are indeed abnormal.

Therefore, the first step will be to check that both datasets are composed of the same patients, ensuring consistency between the detection of the anomalous cases and their subsequent evidencing.

After this, it is necessary to make sure that clinical data are not present in the omics dataset, thus avoiding the influence of such data in the detection of the anomalous cases.

Finally, because a genetic algorithm is used, the predictor referring to the classification of patients must be binary, that is, it must have two possible values.
The following encodings can be used:
- Text strings, "Control" (negative case) and "Case" (positive case).
- 0's and 1's, 0 (negative case) and 1 (positive case).

### Step 2. Detection of anomalous cases

This detection is done by combining a genetic algorithm and a machine learning algorithm, it is possible to use a Lasso model or a Random Forest model.

A genetic algorithm is a technique inspired by biological evolution that simulates the process of natural selection and evolution of species, the execution of a genetic algorithm starts with a randomly generated initial solution set.

- In this case, it is assumed that certain patients are misdiagnosed, so the initial solution set will be created by randomly changing 10% of the patient diagnoses in each of these solutions.

In each generation (iteration) of the genetic algorithm, each of the solutions from the set of possible solutions is used as a training dataset for the ML model, which will classify the patients and whose performance will determine the fitness value (score) of the solution.

The way in which the score for each solution is calculated depends on the model used:

- For the Lasso model, the score of each solution is defined by the average of the balanced accuracies obtained after training 5 Lasso models with that solution as training dataset.

- For RF model, the score of each solution is defined by the balanced accuracy obtained after training the RF model with that solution as training dataset.

After this evaluation, a ranking of the solutions is created: 

- The two best solutions are combined to generate a new solution that is added to the set of possible solutions to the problem. 

- Simultaneously, a certain number of the worst solutions are eliminated from the generation, i.e. they are no longer part of the set of possible solutions to the problem.

This process is repeated for the specified number of generations, until a certain predefined threshold is reached in the score of the solutions or until the value of the best score of the set of solutions is repeated for a certain specified number of generations. 

This procedure is used to determine which samples were anomalous. The genetic algorithm locates these anomalous samples and the machine learning algorithm evidences this detection, validating that these samples belong to the alternative group.

### Step 3. Generating reports to evidence the detection of abnormal cases

After the execution, the software indicates which samples are anomalous cases, the last step is to analyse if these cases detected as anomalous really are, providing evidence based on the clinical data of the patients, for which a report is created.

This report is made up of three sections:

- PCA analysis, which graphically shows the distribution of diagnoses once the abnormal cases are reclassified.

- Table of patients' clinical data, in which the clinical data dataset is displayed.

- Ratio-based analysis, in which 3 graphs are shown:
    - The first will show the ratios of the different clinical data between CONTROL and CONTROL patients detected as abnormal (actually they would be CASE).

    - The second will show the ratios of the different clinical data between CASE and CASE patients detected as abnormal (actually they would be CONTROL).

    - The third will show the ratios of the different clinical data between CONTROL and CASE patients once the abnormal cases are reclassified.

## What software we have developed for this purpose

To solve this problem, the MLASDO software (Machine learning based Anomalous Sample Detection on Omics), which is an R package, has been implemented.

In this software it is possible to configure
- The different parameters on the execution of the genetic algorithm.
- The Machine Learning algorithm to use.
- The different parameters on the execution of the Lasso algorithm.
- The different parameters on the execution of the Random Forest algorithm.
- The variable on which you want to detect anomalous cases.
- The clinical data on which to perform the analysis.

In addition, with this software, it is possible to perform the detection and subsequent analysis on your own datasets, although the software includes a couple of sample datasets.

### R packages required
- devtools
- dplyr
- GA
- glmnet
- caret
- doParallel

### Installation
```
    install.packages("devtools")
 
    library(devtools)

    install_github("tomasBernal/MLASDO")

    library(MLASDO)
```

### Examples of use of such software


Default Execution with Lasso model, the default parameters can be observed in the documentation:
```
    MLASDO::detectAnomalies(
        savingName = "DefaultExecutionLasso", 
        mlAlgorithm = "Lasso", 
        classVariable = "Ca.Co.Last", 
        idColumn = "Trial")
```

Default Execution with Random Forest model, the default parameters can be observed in the documentation:
```
    MLASDO::detectAnomalies(
        savingName = "DefaultExecutionLasso", 
        mlAlgorithm = "RF", 
        classVariable = "Ca.Co.Last", 
        idColumn = "Trial")
```

Test Execution, with this execution, we will be able to verify if we have set the parameters related to the data, active predictors, and the class variable correctly:
```
    MLASDO::detectAnomalies(
        savingName = "QuickExecution",
        mlAlgorithm = "Lasso", 
        nIterations = 3, 
        nStopIter = 1, 
        populationSize = 20, 
        classVariable = "Ca.Co.Last", 
        idColumn = "Trial", 
        activePredictors = c("sex", "age_at_baseline", "Mutation", "Ethnicity")
    )
```


Execution with your own data:
```
    MLASDO::detectAnomalies(
        savingName = "ExecutionWithOwnData", 
        mlAlgorithm = "RF", 
        omicDataPath = "./myOmicData.tsv", 
        clinicDataPath = "./myClinicData.tsv", 
        idColumn = "Patient.Id", 
        nIterations = 3, 
        populationSize = 10, 
        classVariable = "Diagnosis", 
        activePredictors = c("sex", "age", "Ethnicity")
    )
```

## Credits
**Project Leader:** Juan Antonio Botia Blaya https://github.com/juanbot

**Main developer:** Tomás Bernal Beltrán https://github.com/tomasBernal

## Contact
**email address:** tomas.bernalb@um.es
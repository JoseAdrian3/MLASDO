#' executeAndPerformAnalysis
#'
#' Methodology based on the combination of a genetic algorithm and a Machine Learning technique to mutate the final
#' diagnosis of patients, detecting anomalies in them.
#'
#' These discrepancies may have their origin in the evolution of the subject itself, leading it from one group to
#' another in biological terms, or they may derive from human error in the labelling of the samples, among other causes.
#'
#' The identification and correction of these erroneous situations becomes essential to preserve the integrity
#' and accuracy of the data.
#'
#' @param justAnalysis Bool | Indicates whether to perform the analysis directly (TRUE) or to run the genetic algorithm (FALSE). Default value: FALSE.
#' @param geneticPath String | Path to genetic algorithm object.
#' @param solutionPath String | Path to genetic algorithm solution.
#' @param bestModelAfterDetectionPath String | Path to the best model obtained after the detection.
#' @param lassoPredictorsPath String | Path to the mean number of predictors selected by Lasso in each generation.
#'
#' @param mlAlgorithm String | Machine Learning algorithm to be applied, the options are: Lasso or RF (Random Forest).
#'
#' @param predictorsToSelect Integer | Number of predictors to be selected from the most important predictors ranked by the RF model. This parameter is a integer number between 1 and the total number of predictors in the data. Default value: 15.
#' @param numModelExecutions Integer | Number of times the Lasso algorithm is executed. Default value: 5.
#' @param numTrees Integer | Number of trees of the Random Forest model. Default value: 100.
#' @param mtry Integer | Number of predictors that are evaluated at each partition (node) of each tree. Default value: 225.
#' @param splitRule String | This is the rule used by the algorithm to select the predictor and the optimal value to separate a node into two branches during the tree construction. Default value: gini.
#' @param sampleFraction Decimal | Fraction of the training data that will be used to create each of the trees in the forest. Default value: 1.
#' @param maxDepth Integer | Maximum height of each tree in the forest. Default value: 4.
#' @param minNodeSize Integer | Minimum number of observations required in a node to be able to split it. Default value: 30.
#'
#' @param bestLambda Bool | It indicates when to perform cv to find the best lambda for the Lasso model (TRUE) and when not (FALSE). If the value is FALSE, the best lambda found for the Baseline Lasso model will be used for all the Lasso models.
#'
#' @param omicDataPath String | Path to omic data. If the user does not specify a path, the sample data will be used.
#' @param clinicDataPath String | Path to clinic data. If the user does not specify a path, the sample data will be used.
#' @param omicPredictorsToIgnore Array of Strings | Variables to be removed from the omic dataset. These will not be taken into account in the execution.
#' @param clinicPredictorsToIgnore Array of Strings | Variables to be removed from the clinic dataset. These will not be taken into account in the execution.
#'
#' @param idColumn String | Variable that indicates the identifier of each patient in both datasets. If the user does not specify a path to his own data, the value for the sample data, Trial, will be used.
#' @param activePredictors Array of Strings | Predictors on which the study of the ratios will be conducted after the genetic algorithm has been performed. Default value: All the predictors in clinic data, except classVariable and idColumn.
#' @param numericActivePredictors Array of Strings | Numerical predictors of the active predictor list. Default value: All predictors in the active predictor list that return TRUE in the is.numeric function or return TRUE in the is.integer function and have 7 or more distinct values.
#' @param categoricActivePredictors Array of Strings | Categoric predictors of the active predictor list. Default value: All predictors in the active predictor list that return FALSE in the is.numeric function or return TRUE in the is.integer function and have less than 7 distinct values.
#' @param classVariable String | Target variable, which must be binary, meaning it has two possible values. If the user does not specify a path to his own data, the value for the sample data, Ca.Co.Last, will be used.
#' @param savingName String | Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#'
#' @param maximumChangePercentage Decimal | Percentage (expressed as a fraction) indicating the maximum percentage of samples that can be changed. Default value: 0.2 (20%).
#'
#' @param nCores Integer | Number of cores to be used in parallelization. Default value: 6.
#' @param partitionPercentage Decimal | Percentage (expressed as a fraction) with which the data will be split into a training and test set. Default value: 0.9 (90%).
#' @param nIterations Integer | Number of iterations (generations) the genetic algorithm will perform. Default value: 200.
#' @param nStopIter Integer | Number of iterations after which the algorithm will stop if all of them have the same fitness value. Default value: 25.
#' @param populationSize Integer | Number of solutions that will be part of the initial population. Default value: 150.
#' @param diagnosticChangeProbability Decimal | Percentage (expressed as a fraction) indicating the probability of each gene in the solutions to be changed. Default value: 0.1 (10%).
#' @param crossoverOperator String | Name of the crossover operator used in the genetic algorithm. Default value: Single Point Crossover.
#' @param crossoverProbability Decimal | Percentage (expressed as a fraction) indicating the probability of crossover occurrence. Default value: 0.8 (80%).
#' @param selectionOperator String | Name of the selection operator used in the genetic algorithm. Default value: Tournament Selection.
#' @param mutationOperator String | Name of the mutation operator used in the genetic algorithm. Default value: Random Mutation.
#' @param mutationProbability Decimal | Percentage (expressed as a fraction) indicating the probability of mutation occurrence. Default value: 0.1 (10%).
#' @param seed Integer | Seed used for the creation of training and test sets. Default value: 1234.
#'
#' @param pcaAlpha Decimal | Alpha used for the points that don't change in the pca plot. Default value: 0.2.
#' @param pcaSize Decimal | Size used for the points that change in the pca plot. Default value: 1.4.
#'
#' @param oddsDefaultValues String | Parameter indicating for each categorical feature its default value, both feature and value separated by ':'
#'
#' @export
#'
#' @examples
#'
#' MLASDO::detectAnomalies(savingName = "DefaultExecutionLasso", mlAlgorithm = "Lasso",)
#'
#' MLASDO::detectAnomalies(savingName = "DefaultExecutionRF", mlAlgorithm = "RF")
#'
#' MLASDO::detectAnomalies(savingName = "QuickExecution", mlAlgorithm = "Lasso", nIterations = 3, nStopIter = 2, populationSize = 20, activePredictors = c("sex", "age", "Mutation", "Ethnicity"))
#'
#' MLASDO::detectAnomalies(savingName = "ExecutionWithOwnData", mlAlgorithm = "RF", predictorsToSelect = 0.3, omicDataPath = "./myOmicData.tsv", clinicDataPath = "./myClinicData.tsv", idColumn = "Patient.Id", nIterations = 3, populationSize = 10, classVariable = "Diagnosis", activePredictors = c("sex", "age", "Ethnicity"), numericActivePredictors = c("age"), categoricActivePredictors = c("sex", "Ethnicity"))
#'
#' MLASDO::detectAnomalies(justAnalysis = TRUE, mlAlgorithm = "Lasso", geneticPath = "GA.rds", solutionPath = "GA_solution.rds", bestModelAfterDetectionPath = "GA_Best_Model.rds", lassoPredictorsPath = "GA_Lasso_Predictors.rds", savingName = "ExecutionWithOwnData", omicDataPath = "./myOmicData.tsv", clinicDataPath = "./myClinicData.tsv",idColumn = "Patient.Id",classVariable = "Diagnosis", activePredictors = c("sex", "age", "Ethnicity"), oddsDefaultValues <- c("mutation:Healthy Control"))

detectAnomalies <- function(
    justAnalysis = FALSE,
    solutionPath = "",
    geneticPath = "",
    bestModelAfterDetectionPath = "",
    lassoPredictorsPath = "",
    mlAlgorithm,
    predictorsToSelect = 15,
    numModelExecutions = 5,
    numTrees = 100,
    mtry = 225,
    splitRule = "gini",
    sampleFraction = 1,
    maxDepth = 4,
    minNodeSize = 30,
    bestLambda = FALSE,
    omicDataPath = "",
    clinicDataPath = "",
    omicPredictorsToIgnore = NULL,
    clinicPredictorsToIgnore = NULL,
    idColumn,
    activePredictors = NULL,
    numericActivePredictors = NULL,
    categoricActivePredictors = NULL,
    classVariable,
    savingName = "",
    maximumChangePercentage = 0.2,
    nCores = 6,
    partitionPercentage = 0.9,
    nIterations = 200,
    nStopIter = 25,
    populationSize = 150,
    diagnosticChangeProbability = 0.1,
    crossoverOperator = "gabin_spCrossover",
    crossoverProbability = 0.8,
    selectionOperator = "gabin_tourSelection",
    mutationOperator = "gabin_raMutation",
    mutationProbability = 0.1,
    seed = 1234,
    pcaAlpha = 0.2,
    pcaSize = 1.4,
    oddsDefaultValues = character(0)
){


  #### CHECKING REQUIRED LIBRARIES ####
  if (!require("GA", character.only = TRUE)) {
    return("The package GA is not installed")
  }

  if (!require("ranger", character.only = TRUE)) {
    return("The package ranger is not installed")
  }

  if (!require("glmnet", character.only = TRUE)) {
    return("The package glmnet is not installed")
  }

  if (!require("caret", character.only = TRUE)) {
    return("The package caret is not installed")
  }

  if (!require("doParallel", character.only = TRUE)) {
    return("The package doParallel is not installed")
  }

  if (!require("parallel", character.only = TRUE)) {
    return("The package parallel is not installed")
  }

  if (!require("iterators", character.only = TRUE)) {
    return("The package iterators is not installed")
  }

  if (!require("foreach", character.only = TRUE)) {
    return("The package foreach is not installed")
  }

  if (!require("dplyr", character.only = TRUE)) {
    return("The package dplyr is not installed")
  }

  if (!require("plotly", character.only = TRUE)) {
    return("The package plotly is not installed")
  }

  if (!require("lazyeval", character.only = TRUE)) {
    return("The package lazyeval is not installed")
  }

  if (!require("ggtext", character.only = TRUE)) {
    return("The package ggtext is not installed")
  }

  if (!require("DT", character.only = TRUE)) {
    return("The package DT is not installed")
  }

  if (!require("jsmodule", character.only = TRUE)) {
    return("The package DT is not installed")
  }

  if (!require("ggpubr", character.only = TRUE)) {
    return("The package DT is not installed")
  }

  #### CHECKING INITIAL STEP ####

  if(justAnalysis & solutionPath == ""){
    return("If you want to perform the analysis only, you must indicate the path to the final solution of the genetic algorithm.")
  }

  if(justAnalysis & geneticPath == ""){
    return("If you want to perform the analysis only, you must indicate the path to the genetic algorithm object.")
  }

  if(justAnalysis & bestModelAfterDetectionPath == ""){
    return("If you want to perform the analysis only, you must indicate the path to the best model obtained after the detection.")
  }

  if(justAnalysis & mlAlgorithm == "Lasso" & lassoPredictorsPath == ""){
    return("If you want to perform the analysis of a Lasso model, you must indicate the path to the mean number of predictors selected in each generation.")
  }

  #### CHECKING ML ALGORITHM ####
  if(mlAlgorithm != "Lasso" & mlAlgorithm != "RF"){
    return("The Machine Learning algorithms that can be applied are: Lasso or RF (Random Forest)")
  }


  #### DATA READING ####

  if (omicDataPath == ""){
    omicRoute <- system.file("data", "omicData.tsv", package = "MLASDO")
    classVariable <- "Ca.Co.Last"
    idColumn <- "Trial"
  } else {
    omicRoute <- omicDataPath
  }

  if (clinicDataPath == ""){
    clinicRoute <- system.file("data", "clinicData.tsv", package = "MLASDO")
    classVariable <- "Ca.Co.Last"
    idColumn <- "Trial"
  } else {
    clinicRoute <- clinicDataPath
  }

  if (grepl("\\.tsv", omicRoute)) {

      omicData <- read.table(omicRoute, header = TRUE, sep = "\t", row.names = 1)

  } else if (grepl("\\.csv", omicRoute)){
      omicData <- read.csv(omicRoute)
  } else {

    return("The only valid formats for omic data are: .tsv and .csv")
  }

  if (grepl("\\.tsv", clinicRoute)) {

    clinicData <- read.table(clinicRoute, header = TRUE, sep = "\t", row.names = 1)

  } else if (grepl("\\.csv", clinicRoute)){
    clinicData <- read.csv(clinicRoute)
  } else {

    return("The only valid formats for clinic data are: .tsv and .csv")
  }

  #### RELEVEL CLINIC DATA ####

  # Variable to store the features with a default value
  defaultOddsPredictors <- character(0)

  # Function to convert a dataframe's (df) feature (predictor) to a factor with the first value as the default value (default_value)
  convert_to_factor <- function(df, predictor, default_value) {
    df[[predictor]] <- factor(df[[predictor]], levels = c(default_value, setdiff(unique(df[[predictor]]), default_value)))
    return(df)
  }

  # From each string, separate by ':', with the left side being the predictor and the right side being the desired default value (future first level of the feature)
  for (entry in oddsDefaultValues) {
    parts <- unlist(strsplit(entry, ":"))
    predictor <- parts[1]
    default_value <- parts[2]

    # Checking feature is not a numeric feature
    if(!predictor %in% numericActivePredictors){

      # Updeting variable of features with default predictors
      defaultOddsPredictors <- c(defaultOddsPredictors, predictor)

      clinicData <- convert_to_factor(clinicData, predictor, default_value)

    } else {
      print("You have attempted to convert a numerical predictor to have a default value")
    }
  }

  #### RELEVEL CLINIC DATA ####

  #### CHECKING PARAMETERS ####

  if(!is.null(omicPredictorsToIgnore)){

    validColumnsToRemove <- names(omicData) %in% omicPredictorsToIgnore

    if(any(validColumnsToRemove)){
      omicData <- omicData[, !validColumnsToRemove]
    }

  }

  if(!is.null(clinicPredictorsToIgnore)){

    validColumnsToRemove <- names(clinicData) %in% clinicPredictorsToIgnore

    if(any(validColumnsToRemove)){
      clinicData <- clinicData[, !validColumnsToRemove]
    }

  }


  if(is.null(activePredictors) & !is.null(numericActivePredictors) & !is.null(categoricActivePredictors)){

    activePredictors <- c(numericActivePredictors, categoricActivePredictors)

  } else if(is.null(activePredictors)){

    activePredictors <- names(clinicData)

    activePredictors <- activePredictors[activePredictors != classVariable]

    activePredictors <- activePredictors[activePredictors != idColumn]
  }

  # Checking if the patients in both datasets are the same
  omicIds <- sort(omicData[[idColumn]])
  clinicIDs <- sort(clinicData[[idColumn]])

  if(!identical(omicIds, clinicIDs)){
    return("The IDs of the patients don't match in both datasets")
  }

  omicClass <- omicData[c(idColumn, classVariable)]
  clinicClass <- clinicData[c(idColumn, classVariable)]

  omicClass <- omicClass[order(omicClass[[idColumn]]), ]
  clinicClass <- clinicClass[order(clinicClass[[idColumn]]), ]

  rownames(omicClass) <- NULL
  rownames(clinicClass) <- NULL

  if(!identical(omicClass, clinicClass)){
    return("The diagnostic of the patients don't match in both datasets")
  }

  if(predictorsToSelect == 0 | predictorsToSelect > ncol(omicData)){
    return("The number of predictors to select must be a number between 1 and the total number of predictors")
  }

  # First, check if the selected variable as a predictor is included in the variables to be studied
  if (classVariable %in% activePredictors) {
    return("The selected variable for prediction cannot be included in the variables chosen for analysis")
  }

  if (idColumn %in% activePredictors) {
    return("The variable that identifies the patients cannot be included in the variables chosen for analysis")
  }

  if (idColumn %in% omicPredictorsToIgnore) {
    return("The variable that identifies the patients cannot be included in the variables chosen to ignore")
  }

  if (idColumn %in% clinicPredictorsToIgnore) {
    return("The variable that identifies the patients cannot be included in the variables chosen to ignore")
  }

  # After this, check if the selected variable as the class variable for prediction is valid
  if (!classVariable %in% names(omicData)) {
    return("The selected variable for prediction does not exist in the omic dataframe")
  }

  if (!classVariable %in% names(clinicData)) {
    return("The selected variable for prediction does not exist in the clinic dataframe")
  }

  if (!idColumn %in% names(clinicData)) {
    return("The variable that identifies the patients does not exist in the clinic dataframe")
  }

  if (!idColumn %in% names(omicData)) {
    return("The variable that identifies the patients does not exist in the clinic omic dataframe")
  }

  # Next, check if the selected variable as the class variable for prediction is binary
  if (length(unique(omicData[[classVariable]])) > 2) {
    return("The selected variable for prediction is not binary")
  }

  if (length(unique(clinicData[[classVariable]])) > 2) {
    return("The selected variable for prediction is not binary")
  }

  # Checking if the class variable is well coded
  if(all(!(omicData[[classVariable]] %in% c(0, 1))) & all(!(omicData[[classVariable]] %in% c("Case", "Control")))){
    return("The class variable is not well coded")
  }

  if(all(!(clinicData[[classVariable]] %in% c(0, 1))) & all(!(clinicData[[classVariable]] %in% c("Case", "Control")))){
    return("The class variable is not well coded")
  }

  # Checking if the coding of the class variable in the omic dataset is correct
  if (0 %in% unique(omicData[[classVariable]])) {

    cat("The class variable is coded with 0, those values are changed to Control.\n")

    omicData[[classVariable]][omicData[[classVariable]] == 0] <- "Control"

    clinicData[[classVariable]][clinicData[[classVariable]] == 0] <- "Control"

  }

  # Checking if the coding of the class variable in the omic dataset is correct
  if (1 %in% unique(omicData[[classVariable]])) {

    cat("The class variable is coded with 1, those values are changed to Case.\n")

    omicData[[classVariable]][omicData[[classVariable]] == 1] <- "Case"

    clinicData[[classVariable]][clinicData[[classVariable]] == 1] <- "Case"

  }

  # Checking if the coding of the class variable in the clinic dataset is correct
  if (0 %in% unique(clinicData[[classVariable]])) {

    cat("The class variable is coded with 0, those values are changed to Control.\n")

    omicData[[classVariable]][omicData[[classVariable]] == 0] <- "Control"

    clinicData[[classVariable]][clinicData[[classVariable]] == 0] <- "Control"

  }


  # Checking if the coding of the class variable in the clinic dataset is correct
  if (1 %in% unique(clinicData[[classVariable]])) {

    cat("The class variable is coded with 1, those values are changed to Case.\n")

    omicData[[classVariable]][omicData[[classVariable]] == 1] <- "Case"

    clinicData[[classVariable]][clinicData[[classVariable]] == 1] <- "Case"

  }

  # Finally, check if the variables selected for ratio analysis are valid
  invalidPredictors <- activePredictors[!(activePredictors %in% names(clinicData))]

  if (length(invalidPredictors) > 0) {
    cat("The following variables selected for analysis do not exist in the dataframe:\n")
    cat(paste(invalidPredictors, collapse = ", "))
    cat("\n")

    activePredictors <- activePredictors[!(activePredictors %in% invalidPredictors)]
  }

  if(is.null(categoricActivePredictors) & !is.null(numericActivePredictors)){
    cat("You have only indicated which active predictors are numeric, all other active predictors will be treated as categorical.\n")

    restOfPredictors <- activePredictors[!(activePredictors %in% numericActivePredictors)]

    categoricActivePredictors <- restOfPredictors

    invalidPredictors <- numericActivePredictors[!(numericActivePredictors %in% names(clinicData))]

    numericActivePredictors <- numericActivePredictors[!(numericActivePredictors %in% invalidPredictors)]



    cat("Active numerical predictors: ")
    cat(paste(numericActivePredictors, collapse = ", "))
    cat("\n")
    cat("Active categoric predictors: ")
    cat(paste(categoricActivePredictors, collapse = ", "))
    cat("\n")
  }

  if(!is.null(categoricActivePredictors) & is.null(numericActivePredictors)){
    cat("You have only indicated which active predictors are categorical, all other active predictors will be treated as numeric.\n")

    restOfPredictors <- activePredictors[!(activePredictors %in% categoricActivePredictors)]

    numericActivePredictors <- restOfPredictors

    invalidPredictors <- categoricActivePredictors[!(categoricActivePredictors %in% names(clinicData))]

    categoricActivePredictors <- categoricActivePredictors[!(categoricActivePredictors %in% invalidPredictors)]

    cat("Active numerical predictors: ")
    cat(paste(numericActivePredictors, collapse = ", "))
    cat("\n")
    cat("Active categoric predictors: ")
    cat(paste(categoricActivePredictors, collapse = ", "))
    cat("\n")
  }


  if(!is.null(categoricActivePredictors) & !is.null(numericActivePredictors)){

    repeatedPredictors <- intersect(numericActivePredictors, categoricActivePredictors)

    if (length(repeatedPredictors) > 0) {
      cat("The following variables selected for analysis are listed as categorical and numerical variables:\n")
      cat(paste(repeatedPredictors, collapse = ", "))
      return("Please remove repeated variables.")
    }

  }

  # If the name is empty, we need to establish a default one
  if (savingName == ""){
    savingName <- format(Sys.time(), "%d-%m-%H-%M")
  }


  #### DATA PROCESSING ####
  omicData <- omicData[,!(names(omicData) %in% activePredictors)]

  omicsFiltered <- omicData[,!(names(omicData) %in% classVariable)]

  validClinicData <- merge(clinicData, omicsFiltered, by = idColumn)

  validClinicData <- validClinicData[, (names(validClinicData) %in% c(activePredictors, idColumn, classVariable))]

  dir.create(savingName)
  dir.create(paste(savingName, "analysisData", sep = "/"))

  set.seed(seed)

  subsetTrain <- createDataPartition(y = omicData[[classVariable]], p = partitionPercentage, list = FALSE)[,1]

  # Generating Baseline model
  omic <- omicData

  # We need to remove the column that indicates the patient id
  omic[[idColumn]] <- NULL

  omic[[classVariable]] <- ifelse(omic[[classVariable]] == "Case", 1, 0)

  #### CREATION OF TRAIN AND TEST SETS ####
  omicTrain <- omic[subsetTrain,]

  omicTrainDiagnosis <- omicTrain[[classVariable]]
  omicTrain[[classVariable]] <- NULL

  omicTest <- omic[-subsetTrain,]

  omicTestDiagnosis <- omicTest[[classVariable]]
  omicTest[[classVariable]] <- NULL

  baselinePrecision <- 0
  baselinePredictors <- 0

  bestBaselineModel <- NULL
  bestBaselineModelBA <- 0

  if(mlAlgorithm == "Lasso"){

    # This vector will store the balanced means of the numModelExecutions executions
    balancedAccValues <- vector(length=numModelExecutions)

    numPredictors <- list()

    set.seed(seed)

    # Train the Lasso model numModelExecutions times
    for (i in 1:numModelExecutions) {

      # Train the Lasso model
      model <- cv.glmnet(as.matrix(omicTrain), as.matrix(omicTrainDiagnosis), alpha = 1, family = "binomial", type.measure = "class", nfolds = 10)

      coefficients <- coef(model, s = model$lambda.min)

      indexes <- coefficients@i

      posZero <- which(indexes == 0)

      indexes <- indexes[-posZero]

      for(j in 1:length(indexes)){

        predName <- names(omicTrain)[[indexes[[j]]]]

        if (!exists(predName, where = numPredictors)) {

          numPredictors[[predName]] <- 1

        }

      }

      # Use the model to predict on the test set
      modelPrediction <- predict(model, newx = as.matrix(omicTest), alpha = 1, s = "lambda.min", type = "class")

      # Get the confusion matrix
      cfModel <- confusionMatrix(as.factor(as.integer(modelPrediction)), as.factor(omicTestDiagnosis))

      specificity <- ifelse(is.na(as.numeric(cfModel$byClass["Specificity"])), 0,  as.numeric(cfModel$byClass["Specificity"]))

      sensitivity <- ifelse(is.na(as.numeric(cfModel$byClass["Sensitivity"])), 0,  as.numeric(cfModel$byClass["Sensitivity"]))

      # Save the balanced mean obtained in this iteration
      balancedAccValues[i] <- (specificity + sensitivity) / 2

      if(balancedAccValues[i] > bestBaselineModelBA){

        bestBaselineModelBA <- balancedAccValues[i]
        bestBaselineModel <- model
      }

    }

    baselinePrecision <- mean(balancedAccValues)

    baselinePredictors <- length(numPredictors)

  } else if(mlAlgorithm == "RF"){

    # This vector will store the balanced means of the numModelExecutions executions
    balancedAccValues <- vector(length=numModelExecutions)

    # Train the Lasso model numModelExecutions times
    for (i in 1:numModelExecutions) {

      # Creating the ranger model
      model <- ranger(
        x = omicTrain,
        y = omicTrainDiagnosis,
        num.trees = numTrees,
        mtry = mtry,
        splitrule = splitRule,
        importance = "impurity", # In order to obtain the important variables in the prediction
        sample.fraction = sampleFraction,
        max.depth = maxDepth,
        min.node.size = minNodeSize,
        classification = TRUE,
        seed = seed
      )

      # Use the model to predict on the test set
      modelPrediction <- predict(model, omicTest)$predictions

      # Get the confusion matrix
      cfModel <-
        confusionMatrix(as.factor(as.integer(modelPrediction)), as.factor(omicTestDiagnosis))

      specificity <- ifelse(is.na(as.numeric(cfModel$byClass["Specificity"])), 0,  as.numeric(cfModel$byClass["Specificity"]))

      sensitivity <- ifelse(is.na(as.numeric(cfModel$byClass["Sensitivity"])), 0,  as.numeric(cfModel$byClass["Sensitivity"]))

      balancedAccValues[i] <- (specificity + sensitivity) / 2

      if(balancedAccValues[i] > bestBaselineModelBA){

        bestBaselineModelBA <- balancedAccValues[i]
        bestBaselineModel <- model
      }
    }

    baselinePrecision <- mean(balancedAccValues)

  }

  if(!justAnalysis){

    dir.create(paste(savingName, "geneticAlgorithm", sep = "/"))

    cat("\nExecuting the genetic algorithm.\n")
    MLASDO::executeGA(
      savingName = savingName,
      omicData = omicData,
      subsetTrain = subsetTrain,
      classVariable = classVariable,
      idColumn = idColumn,
      mlAlgorithm = mlAlgorithm,
      numModelExecutions = numModelExecutions,
      bestLambda = ifelse(bestLambda == TRUE, FALSE, model$lambda.min),
      predictorsToSelect = predictorsToSelect,
      numTrees = numTrees,
      mtry = mtry,
      splitRule = splitRule,
      sampleFraction = sampleFraction,
      maxDepth = maxDepth,
      minNodeSize = minNodeSize,
      nIterations = nIterations,
      nStopIter = nStopIter,
      populationSize = populationSize,
      diagnosticChangeProbability = diagnosticChangeProbability,
      crossoverOperator = crossoverOperator,
      crossoverProbability = crossoverProbability,
      selectionOperator = selectionOperator,
      mutationOperator = mutationOperator,
      mutationProbability = mutationProbability,
      maximumChangePercentage = maximumChangePercentage,
      nCores = nCores,
      seed = seed
      )
  }

  # Reading the GA solution
  name <- paste("GA", savingName, sep="_")

  if(justAnalysis){

    geneticAlgorithm <- readRDS(geneticPath)
    solutionGA <- readRDS(solutionPath)
    bestModelAfterDetection <- readRDS(bestModelAfterDetectionPath)

  } else {

    dirPath <- paste(savingName, "geneticAlgorithm", name, sep = "/")

    bestModelAfterDetectionPath <- paste(dirPath, "Best_Model_After_Detection.rds", sep="_")
    bestModelAfterDetection <- readRDS(bestModelAfterDetectionPath)

    solutionPath <- paste(dirPath, "Solution.rds", sep="_")
    solutionGA <- readRDS(solutionPath)

    gaPath <- paste(dirPath, ".rds", sep="")
    geneticAlgorithm <- readRDS(gaPath)

  }

  # Creating a copy of the original diagnosis
  changedDiagnoses <- omicData[[classVariable]]

  # I create two text strings to represent the changed diagnoses
  firstGroup = paste("Case", "Control", sep = "2")
  secondGroup = paste("Control", "Case", sep = "2")

  # Converting original diagnoses from strings to integers
  changedDiagnoses <- ifelse(changedDiagnoses == "Case", 1, 0)

  # I apply the solution from the genetic algorithm to the original diagnoses
  # I use an XOR function for this purpose
  changedDiagnoses[subsetTrain] <- bitwXor(changedDiagnoses[subsetTrain], solutionGA)

  # I obtain the indices where the genetic algorithm solution indicates changes in someone's diagnosis
  changeIndices <- which(solutionGA == 1)

  # Iterating through the indices that have been modified
  for(i in changeIndices){

    # If the gene, after modification, had the value 1
    # It means it had the value 0 (Control)
    if(changedDiagnoses[subsetTrain[[i]]] == 1){

      # Therefore, I set it to 2
      changedDiagnoses[subsetTrain[[i]]] <- secondGroup

      # If not, it means it had the value 1 (Case)
    } else{

      # Therefore, I set it to 3
      changedDiagnoses[subsetTrain[[i]]] <- firstGroup
    }
  }

  changedDiagnoses <- ifelse(changedDiagnoses == 0, "Control", changedDiagnoses)
  changedDiagnoses <- ifelse(changedDiagnoses == 1, "Case", changedDiagnoses)

  # After this, I have the following values:
  #     0 -> Gene remained with the value 0 of the variable being predicted
  #     1 -> Gene remained with the value 1 of the variable being predicted
  #     2 -> Gene changed from the value 0 to the value 1
  #     3 -> Gene changed from the value 1 to the value 0

  changedOmicData <- omicData

  changedClinicData <- validClinicData


  changedClinicData[[classVariable]] <- changedDiagnoses
  changedOmicData[[classVariable]] <- changedDiagnoses

  dirPath <- paste(savingName, "geneticAlgorithm", name, sep = "/")
  selectedPredictorsPath <- paste(dirPath, "Predictors_Importance.tsv", sep="_")

  selectedPredictors <- read.table(selectedPredictorsPath, header = TRUE, sep = "\t", row.names = 1)

  baselineModelPath <- paste(dirPath, "Best_Model_Before_Detection.rds", sep="_")
  saveRDS(bestBaselineModel, file = baselineModelPath)

  omicDataShow <- omicData
  omicDataShow[[classVariable]] <- changedDiagnoses

  selectedOmicPredictors <- omicDataShow[, c(idColumn, rownames(selectedPredictors), classVariable)]

  dirPath <- paste(savingName, "analysisData", name, sep = "/")
  selectedDataPath <- paste(dirPath, "Clinic_Selected_Predictors.tsv", sep="_")

  write.table(selectedOmicPredictors, selectedDataPath, row.names = T, col.names = T, sep =  '\t')


  if(justAnalysis){
    cat("\nPerforming PCA analysis.\n")
  } else {
    cat("Performing PCA analysis.\n")
  }

  MLASDO::performPCAAnalysis(
    savingName = savingName,
    mlAlgorithm = mlAlgorithm,
    changedOmicData = changedOmicData,
    selectedData = selectedOmicPredictors,
    idColumn = idColumn,
    classVariable = classVariable
  )

  if(!is.null(categoricActivePredictors) & !is.null(numericActivePredictors)){

    cat("Performing ratio analysis with the variable classification specified by the user.\n")
    MLASDO::performRatioAnalysisUserVariableClassification(
      savingName = savingName,
      changedClinicData = changedClinicData,
      firstGroup = firstGroup,
      secondGroup = secondGroup,
      activePredictors = activePredictors,
      categoricActivePredictors = categoricActivePredictors,
      numericActivePredictors = numericActivePredictors,
      classVariable = classVariable
    )

  } else {

    cat("Performing ratio analysis with automated variable classification.\n")
    MLASDO::performRatioAnalysisAutomatedVariableClassification(
      savingName = savingName,
      changedClinicData = changedClinicData,
      firstGroup = firstGroup,
      secondGroup = secondGroup,
      activePredictors = activePredictors,
      classVariable = classVariable
    )

  }

  if(!is.null(categoricActivePredictors) & !is.null(numericActivePredictors)){

    cat("Performing stadistic analysis with the variable classification specified by the user.\n")
    MLASDO::performStadisticAnalysisUserVariableClassification(
      savingName = savingName,
      changedClinicData = changedClinicData,
      firstGroup = firstGroup,
      secondGroup = secondGroup,
      activePredictors = activePredictors,
      categoricActivePredictors = categoricActivePredictors,
      numericActivePredictors = numericActivePredictors,
      classVariable = classVariable,
      defaultOddsPredictors = defaultOddsPredictors
    )

  } else {

    cat("Performing stadistic analysis with automated variable classification.\n")
    MLASDO::performStadisticAnalysisAutomatedVariableClassification(
      savingName = savingName,
      changedClinicData = changedClinicData,
      firstGroup = firstGroup,
      secondGroup = secondGroup,
      activePredictors = activePredictors,
      classVariable = classVariable,
      defaultOddsPredictors = defaultOddsPredictors
    )

  }


  bestAfterDetectionCM <- NULL
  bestBaselineModelCM <- NULL



  if(mlAlgorithm == "Lasso"){

    bestModelAfterDetectionPrediction <- predict(bestModelAfterDetection, newx = as.matrix(omicTest), alpha = 1, s = "lambda.min", type = "class")
    bestAfterDetectionCM <- confusionMatrix(as.factor(as.integer(bestModelAfterDetectionPrediction)), as.factor(omicTestDiagnosis))

    bestBaselineModelPrediction <- predict(bestBaselineModel, newx = as.matrix(omicTest), alpha = 1, s = "lambda.min", type = "class")
    bestBaselineModelCM <- confusionMatrix(as.factor(as.integer(bestBaselineModelPrediction)), as.factor(omicTestDiagnosis))

  } else if (mlAlgorithm == "RF"){

    bestModelAfterDetectionPrediction <- predict(bestModelAfterDetection, omicTest)$predictions
    bestAfterDetectionCM <-
      confusionMatrix(as.factor(as.integer(bestModelAfterDetectionPrediction)), as.factor(omicTestDiagnosis))

    bestBaselineModelPrediction <- predict(bestBaselineModel, omicTest)$predictions
    bestBaselineModelCM <-
      confusionMatrix(as.factor(as.integer(bestBaselineModelPrediction)), as.factor(omicTestDiagnosis))

  }

  rownames(selectedOmicPredictors) <- selectedOmicPredictors[[idColumn]]
  selectedOmicPredictors[[idColumn]] <- NULL

  rownames(changedClinicData) <- changedClinicData[[idColumn]]
  changedClinicData[[idColumn]] <- NULL

  cat("Compiling Markdown file.\n")
  MLASDO::compileMarkdown(
    savingName = savingName,
    justAnalysis = justAnalysis,
    geneticAlgorithm = geneticAlgorithm,
    mlAlgorithm = mlAlgorithm,
    numModelExecutions = numModelExecutions,
    lassoPredictorsPath = lassoPredictorsPath,
    baselinePrecision = baselinePrecision,
    baselinePredictors = baselinePredictors,
    originalDiagnosis = omicData[[classVariable]],
    clinicData = changedClinicData,
    categoricActivePredictors = categoricActivePredictors,
    numericActivePredictors = numericActivePredictors,
    selectedData = selectedOmicPredictors,
    classVariable = classVariable,
    idColumn = idColumn,
    predictorsToSelect = predictorsToSelect,
    numTrees = numTrees,
    mtry = mtry,
    splitRule = splitRule,
    sampleFraction = sampleFraction,
    maxDepth = maxDepth,
    minNodeSize = minNodeSize,
    bestLambda = bestLambda,
    nIterations = nIterations,
    nStopIter = nStopIter,
    populationSize = populationSize,
    diagnosticChangeProbability = diagnosticChangeProbability,
    crossoverOperator = crossoverOperator,
    crossoverProbability = crossoverProbability,
    selectionOperator = selectionOperator,
    mutationOperator = mutationOperator,
    mutationProbability = mutationProbability,
    nCores = nCores,
    seed = seed,
    bestAfterDetectionCM = bestAfterDetectionCM,
    bestBaselineModelCM = bestBaselineModelCM,
    pcaAlpha = pcaAlpha,
    pcaSize = pcaSize
    )
}

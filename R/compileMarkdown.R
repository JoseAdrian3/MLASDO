#' compileMarkdown
#'
#' This function compiles the markdown file with the analysis results
#'
#' @param savingName String | Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#'
#' @param justAnalysis Bool | Indicates whether to perform the analysis directly (TRUE) or to run the genetic algorithm (FALSE). Default value: FALSE.
#' @param geneticAlgorithm GA | Genetic algorithm object.
#' @param mlAlgorithm String | Machine Learning algorithm to be applied, the options are: Lasso or RF (Random Forest).
#' @param numModelExecutions Integer | Number of times the Lasso algorithm is executed. Default value: 5.
#'
#' @param lassoPredictorsPath String | Path to the mean number of predictors selected by Lasso in each generation.
#' @param baselinePrecision Decimal | Baseline precision obtained with the model before the detection.
#' @param baselinePredictors Integer | Baseline predictors selected with the model before the detection.
#'
#' @param originalDiagnosis Array of Strings | Original diagnostics of the patients.
#' @param clinicData Dataset | Dataset of clinic data that will be used.
#' @param omicPredictorsToIgnore Array of Strings | Variables to be removed from the omic dataset. These will not be taken into account in the execution.
#' @param clinicPredictorsToIgnore Array of Strings | Variables to be removed from the clinic dataset. These will not be taken into account in the execution.
#' @param selectedData Dataset | Dataset of omic data with only the predictors selected by the Lasso model.
#'
#' @param classVariable String | Target variable, which must be binary, meaning it has two possible values. If the user does not specify a path to his own data, the value for the sample data, Ca.Co.Last, will be used.
#' @param idColumn String | Variable that indicates the identifier of each patient in both datasets. If the user does not specify a path to his own data, the value for the sample data, Trial, will be used.
#'
#' @param predictorsToSelect Integer | Number of predictors to be selected from the most important predictors ranked by the RF model. This parameter is a integer number between 1 and the total number of predictors in the data. Default value: 15.
#' @param numTrees Integer | Number of trees of the Random Forest model. Default value: 100.
#' @param mtry Integer | Number of predictors that are evaluated at each partition (node) of each tree. Default value: 225.
#' @param splitRule String | This is the rule used by the algorithm to select the predictor and the optimal value to separate a node into two branches during the tree construction. Default value: gini.
#' @param sampleFraction Decimal | Fraction of the training data that will be used to create each of the trees in the forest. Default value: 1.
#' @param maxDepth Integer | Maximum height of each tree in the forest. Default value: 4.
#' @param minNodeSize Integer | Minimum number of observations required in a node to be able to split it. Default value: 30.
#'
#' @param bestLambda Bool | It indicates when to perform cv to find the best lambda for the Lasso model (TRUE) and when not (FALSE). If the value is FALSE, the best lambda found for the Baseline Lasso model will be used for all the Lasso models.
#'
#' @param nIterations Integer | Number of iterations (generations) the genetic algorithm will perform. Default value: 200.
#' @param nStopIter Integer | Number of iterations after which the algorithm will stop if all of them have the same fitness value. Default value: 25.
#' @param populationSize Integer | Number of solutions that will be part of the initial population. Default value: 150.
#' @param diagnosticChangeProbability Decimal | Percentage (expressed as a fraction) indicating the probability of each gene in the solutions to be changed. Default value: 0.1 (10%).
#' @param crossoverOperator String | Crossover operator used in the genetic algorithm. Default value: Single Point Crossover.
#' @param crossoverProbability Decimal | Percentage (expressed as a fraction) indicating the probability of crossover occurrence. Default value: 0.8 (80%).
#' @param selectionOperator String | Selection operator used in the genetic algorithm. Default value: Tournament Selection.
#' @param mutationOperator String | Mutation operator used in the genetic algorithm. Default value: Random Mutation.
#' @param mutationProbability Decimal | Percentage (expressed as a fraction) indicating the probability of mutation occurrence. Default value: 0.1 (10%).
#'
#' @param nCores Integer | Number of cores to be used in parallelization. Default value: 6.
#'
#' @param seed Integer | Seed used for the creation of training and test sets. Default value: 1234.
#'
#' @param bestAfterDetectionCM Confusion Matrix | Confusion matrix of the best model obtained after the detection.
#' @param bestBaselineModelCM Confusion Matrix | Confusion matrix of the best model obtained before the detection.
#'
#' @param pcaAlpha Decimal | Alpha used for the points that don't change in the pca plot. Default value: 0.2.
#' @param pcaSize Decimal | Size used for the points that change in the pca plot. Default value: 1.3.
#'
#' @export
#'
#' @examples
#'
#' MLASDO::compileMarkdown(savingName = savingName, justAnalysis = justAnalysis, geneticAlgorithm = geneticAlgorithm, mlAlgorithm = mlAlgorithm, numModelExecutions = numModelExecutions, lassoPredictorsPath = lassoPredictorsPath, baselinePrecision = baselinePrecision, baselinePredictors = baselinePredictors, originalDiagnosis = originalDiagnosis, clinicData = clinicData, categoricActivePredictors = categoricActivePredictors, numericActivePredictors = numericActivePredictors, selectedData = selectedData, classVariable = classVariable, idColumn = idColumn, predictorsToSelect = predictorsToSelect, numTrees = numTrees, mtry = mtry, splitRule = splitRule, sampleFraction = sampleFraction, maxDepth = maxDepth, minNodeSize = minNodeSize, bestLambda = bestLambda, nIterations = nIterations, nStopiter = nStopiter, populationSize = populationSize, diagnosticChangeProbability = diagnosticChangeProbability, crossoverOperator = crossoverOperator, crossoverProbability = crossoverProbability, selectionOperator = selectionOperator, mutationOperator = mutationOperator, mutationProbability = mutationProbability, nCores = nCores, seed = seed, bestAfterDetectionCM = bestAfterDetectionCM, bestBaselineModelCM = bestBaselineModelCM, pcaAlpha = pcaAlpha, pcaSize = pcaSize)
#'

compileMarkdown <- function(
    savingName,
    justAnalysis,
    geneticAlgorithm,
    mlAlgorithm,
    numModelExecutions,
    lassoPredictorsPath,
    baselinePrecision,
    baselinePredictors,
    originalDiagnosis,
    clinicData,
    categoricActivePredictors,
    numericActivePredictors,
    selectedData,
    classVariable,
    idColumn,
    predictorsToSelect,
    numTrees,
    mtry,
    splitRule,
    sampleFraction,
    maxDepth,
    minNodeSize,
    bestLambda,
    nIterations,
    nStopIter,
    populationSize,
    diagnosticChangeProbability,
    crossoverOperator,
    crossoverProbability,
    selectionOperator,
    mutationOperator,
    mutationProbability,
    nCores,
    seed,
    bestAfterDetectionCM,
    bestBaselineModelCM,
    pcaAlpha,
    pcaSize
    ){

  name <- paste("GA", savingName, sep="_")

  dirPath <- paste(savingName, "analysisData", name, sep = "/")

  pcaPath <- paste(dirPath, "PCA.tsv", sep="_")
  pca <- read.table(pcaPath, header = TRUE, sep = "\t", row.names = 1)

  pcaPath <- paste(dirPath, "PCA_Selected.tsv", sep="_")
  pcaSelected <- read.table(pcaPath, header = TRUE, sep = "\t", row.names = 1)

  pcaPath <- paste(dirPath, "PCA_Info.tsv", sep="_")
  pcaInfo <- read.table(pcaPath, header = TRUE, sep = "\t", row.names = 1)

  gaPathNumeric <- paste(dirPath, "NumericTable.tsv", sep="_")
  gaPathTotal <- paste(dirPath, "TotalTable.tsv", sep="_")

  numeric <- read.table(gaPathNumeric, header = TRUE, sep = "\t", row.names = 1)
  totalTable <- read.table(gaPathTotal, header = TRUE, sep = "\t", row.names = 1)

  wilcoxPathControls <- paste(dirPath, "WilcoxControls.tsv", sep="_")
  totalWilcoxControls <- read.table(wilcoxPathControls, header = TRUE, sep = "\t", row.names = 1)

  wilcoxPathCases <- paste(dirPath, "WilcoxCases.tsv", sep="_")
  totalWilcoxCases <- read.table(wilcoxPathCases, header = TRUE, sep = "\t", row.names = 1)

  OddsPathControls <- paste(dirPath, "OddsRatiosControls.tsv", sep="_")
  totalOddsControls <- read.table(OddsPathControls, header = TRUE, sep = "\t", row.names = 1)

  OddsPathCases <- paste(dirPath, "OddsRatiosCases.tsv", sep="_")
  totalOddsCases <- read.table(OddsPathCases, header = TRUE, sep = "\t", row.names = 1)

  outputName <- paste("analysisResult_", name, ".html", sep = "")

  outputPath <- paste("./", savingName, "/", sep = "")

  dirPath <- paste(savingName, "geneticAlgorithm", name, sep = "/")

  impPath <- paste(dirPath, "Predictors_Importance.tsv", sep="_")
  predictorsImp <- read.table(impPath, header = TRUE, sep = "\t", row.names = 1)

  lassoPredictors <- NULL

  if(mlAlgorithm == "Lasso"){

    if(justAnalysis){

      lassoPredictors <- readRDS(lassoPredictorsPath)

    } else {

      gaPath <- paste(dirPath, "Predictors.rds", sep="_")
      lassoPredictors <- readRDS(gaPath)
    }
  }

  rmarkdown::render(input = system.file("data", "analysisResult.Rmd", package = "MLASDO"),
                   params = list(
                                 geneticAlgorithm = geneticAlgorithm,
                                 lassoPredictors = lassoPredictors,
                                 baselinePrecision = baselinePrecision,
                                 baselinePredictors = baselinePredictors,
                                 originalDiagnosis = originalDiagnosis,
                                 mlAlgorithm = mlAlgorithm,
                                 clinicData = clinicData,
                                 categoricActivePredictors = categoricActivePredictors,
                                 numericActivePredictors = numericActivePredictors,
                                 selectedData = selectedData,
                                 pcaAnalysis = pca,
                                 pcaAnalysisSelected = pcaSelected,
                                 pcaInfo = pcaInfo,
                                 numericTable = numeric,
                                 totalTable = totalTable,
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
                                 numModelExecutions = numModelExecutions,
                                 bestAfterDetectionCM = bestAfterDetectionCM,
                                 bestBaselineModelCM = bestBaselineModelCM,
                                 predictorsImp = predictorsImp,
                                 pcaAlpha = pcaAlpha,
                                 pcaSize = pcaSize,
                                 totalWilcoxControls = totalWilcoxControls,
                                 totalWilcoxCases = totalWilcoxCases,
                                 totalOddsControls = totalOddsControls,
                                 totalOddsCases = totalOddsCases
                                 ),
                   output_file = outputName,
                   output_dir = outputPath
                   )
}

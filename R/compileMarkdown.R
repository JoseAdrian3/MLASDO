#' compileMarkdown
#'
#' This function compiles the markdown file with the analysis results
#'
#' @param savingName String | Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#'
#' @param geneticAlgorithm GA | Genetic algorithm object.
#' @param numModelExecutions Integer | Number of times the Lasso algorithm is executed. Default value: 5.
#'
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
#' @param numTrees Integer | Number of trees of the Random Forest model. Default value: 100.
#' @param mtry Integer | Number of predictors that are evaluated at each partition (node) of each tree. Default value: 225.
#' @param splitRule String | This is the rule used by the algorithm to select the predictor and the optimal value to separate a node into two branches during the tree construction. Default value: gini.
#' @param sampleFraction Decimal | Fraction of the training data that will be used to create each of the trees in the forest. Default value: 1.
#' @param maxDepth Integer | Maximum height of each tree in the forest. Default value: 4.
#' @param minNodeSize Integer | Minimum number of observations required in a node to be able to split it. Default value: 30.
#'
#'
#' @param nIterations Integer | Number of iterations (generations) the genetic algorithm will perform. Default value: 200.
#' @param diagnosticChangeProbability Decimal | Percentage (expressed as a fraction) indicating the probability of each gene in the solutions to be changed. Default value: 0.1 (10%).
#'
#' @param nCores Integer | Number of cores to be used in parallelization. Default value: 6.
#'
#' @param seed Integer | Seed used for the creation of training and test sets. Default value: 1234.
#'
#' @param bestBaselineModelCM Confusion Matrix | Confusion matrix of the best model obtained before the detection.
#'
#' @param pcaAlpha Decimal | Alpha used for the points that don't change in the pca plot. Default value: 0.2.
#' @param pcaSize Decimal | Size used for the points that change in the pca plot. Default value: 1.3.
#'
#' @export
#'
#' @examples
#'
#' MLASDO::compileMarkdown(savingName = savingName, geneticAlgorithm = geneticAlgorithm, numModelExecutions = numModelExecutions, originalDiagnosis = originalDiagnosis, clinicData = clinicData, categoricActivePredictors = categoricActivePredictors, numericActivePredictors = numericActivePredictors, selectedData = selectedData, classVariable = classVariable, idColumn = idColumn, numTrees = numTrees, mtry = mtry, splitRule = splitRule, sampleFraction = sampleFraction, maxDepth = maxDepth, minNodeSize = minNodeSize, diagnosticChangeProbability = diagnosticChangeProbability, nCores = nCores, seed = seed, bestBaselineModelCM = bestBaselineModelCM, pcaAlpha = pcaAlpha, pcaSize = pcaSize)
#'

compileMarkdown <- function(
    classVariable,
    savingName,
    originalDiagnosis,
    idColumn,
    numTrees,
    mtrysvm,
    splitRule,
    sampleFraction,
    maxDepth,
    minNodeSize,
    diagnosticChangeProbability,
    nCores,
    seed,
    clinicDataSVM,
    originalClinicData,
    categoricActivePredictors,
    numericActivePredictors,
    lambdas,
    C,
    gamma,
    kernel,
    partitionPercentage,
    enrichment,
    absolutePath,
    omicData,
    omic,
    domainTable
){
  ### SVM INFORMATION

  dirPath <- paste(savingName, "svm_data", sep = "/")
  dirPath <- paste(dirPath, savingName, sep = "/")

  svmData <- readRDS(paste(dirPath, "svm_data.rds", sep = "_"))
  svmDataAll <- readRDS(paste(dirPath, "svm_data_all.rds", sep = "_"))
  #svmSolution <- readRDS(paste(dirPath, "svm_solution.rds", sep = "_"))

  ### OUTLIERS INFORMATION ###

  dirPath <- paste(savingName, "bootstrapping_shuffle", sep = "/")
  dirPath <- paste(dirPath, savingName, sep = "/")

  outliersSolution <- readRDS(paste(dirPath, "solution_bootstrapping.rds", sep = "_"))
  thresholdControl <- readRDS(paste(dirPath, "threshold_control_bootstrapping.rds", sep = "_"))
  thresholdCase <- readRDS(paste(dirPath, "threshold_case_bootstrapping.rds", sep = "_"))
  mergeData <- readRDS(paste(dirPath, "merge_data_bootstrapping.rds", sep = "_"))

  ### MODELS INFORMATION SVM

  dirPath <- paste(savingName, "modelssvm", sep = "/")
  dirPath <- paste(dirPath, savingName, sep = "/")

  # LASSO Models
  svm_lasso_after_cf <- readRDS(paste(dirPath, "lasso_after_cf.rds", sep = "_"))
  svm_lasso_after_model <- readRDS(paste(dirPath, "lasso_after_model.rds", sep = "_"))
  svm_lasso_before_cf <- readRDS(paste(dirPath, "lasso_before_cf.rds", sep = "_"))
  svm_lasso_before_model <- readRDS(paste(dirPath, "lasso_before_model.rds", sep = "_"))
  svm_lasso_before_ci <- readRDS(paste(dirPath, "lasso_before_ci.rds", sep = "_"))

  # Random Forest (RF) Models
  svm_RF_after_cf <- readRDS(paste(dirPath, "RF_after_cf.rds", sep = "_"))
  svm_RF_after_model <- readRDS(paste(dirPath, "RF_after_model.rds", sep = "_"))
  svm_RF_before_cf <- readRDS(paste(dirPath, "RF_before_cf.rds", sep = "_"))
  svm_RF_before_model <- readRDS(paste(dirPath, "RF_before_model.rds", sep = "_"))
  svm_RF_before_ci <- readRDS(paste(dirPath, "RF_before_ci.rds", sep = "_"))

  # SVM Linear models
  svm_svmLinear_after_cf <- readRDS(paste(dirPath, "SVMLinear_after_cf.rds", sep = "_"))
  svm_svmLinear_after_model <- readRDS(paste(dirPath, "SVMLinear_after_model.rds", sep = "_"))
  svm_svmLinear_before_cf <- readRDS(paste(dirPath, "SVMLinear_before_cf.rds", sep = "_"))
  svm_svmLinear_before_model <- readRDS(paste(dirPath, "SVMLinear_before_model.rds", sep = "_"))
  svm_svmLinear_before_ci <- readRDS(paste(dirPath, "SVMLinear_before_ci.rds", sep = "_"))

  # SVM Radial models
  svm_svmRadial_after_cf <- readRDS(paste(dirPath, "SVMRadial_after_cf.rds", sep = "_"))
  svm_svmRadial_after_model <- readRDS(paste(dirPath, "SVMRadial_after_model.rds", sep = "_"))
  svm_svmRadial_before_cf <- readRDS(paste(dirPath, "SVMRadial_before_cf.rds", sep = "_"))
  svm_svmRadial_before_model <- readRDS(paste(dirPath, "SVMRadial_before_model.rds", sep = "_"))
  svm_svmRadial_before_ci <- readRDS(paste(dirPath, "SVMRadial_before_ci.rds", sep = "_"))

  ### STADISTIC SVM

  dirPath <- paste(savingName, "stadisticssvm", sep = "/")
  dirPath <- paste(dirPath, savingName, sep = "/")

  wilcoxPathControls <- paste(dirPath, "WilcoxControls.rds", sep="_")
  svmTotalWilcoxControls <- readRDS(wilcoxPathControls)

  OddsPathControls <- paste(dirPath, "OddsRatiosControls.rds", sep="_")
  svmTotalOddsControls <- readRDS(OddsPathControls)

  wilcoxPathCases <- paste(dirPath, "WilcoxCases.rds", sep="_")
  svmTotalWilcoxCases <- readRDS(wilcoxPathCases)

  OddsPathCases <- paste(dirPath, "OddsRatiosCases.rds", sep="_")
  svmTotalOddsCases <- readRDS(OddsPathCases)


  outputName <- paste0(savingName, "/MLASDO_", savingName, ".html")
  outputPath <- savingName

  rmarkdown::render(input = system.file("data", "analysisResult.Rmd", package = "MLASDO"),
                    params = list(
                      classVariable = classVariable,
                      idColumn = idColumn,
                      numTrees = numTrees,
                      mtrysvm = mtrysvm,
                      diagnosticChangeProbability = diagnosticChangeProbability,
                      nCores = nCores,
                      seed = seed,
                      svmData = svmData,
                      svmDataAll = svmDataAll,
                      svm_lasso_after_cf = svm_lasso_after_cf,
                      svm_lasso_after_model = svm_lasso_after_model,
                      svm_lasso_before_cf = svm_lasso_before_cf,
                      svm_lasso_before_model = svm_lasso_before_model,
                      svm_RF_after_cf = svm_RF_after_cf,
                      svm_RF_after_model = svm_RF_after_model,
                      svm_RF_before_cf = svm_RF_before_cf,
                      svm_RF_before_model = svm_RF_before_model,
                      svmTotalWilcoxControls = svmTotalWilcoxControls,
                      svmTotalOddsControls = svmTotalOddsControls,
                      svmTotalWilcoxCases = svmTotalWilcoxCases,
                      svmTotalOddsCases = svmTotalOddsCases,
                      clinicDataSVM = clinicDataSVM,
                      originalDiagnosis = originalDiagnosis,
                      categoricActivePredictors = categoricActivePredictors,
                      numericActivePredictors = numericActivePredictors,
                      svm_lasso_before_ci = svm_lasso_before_ci,
                      svm_RF_before_ci = svm_RF_before_ci,
                      lambdas = lambdas,
                      C = C,
                      gamma = gamma,
                      svm_svmLinear_after_cf = svm_svmLinear_after_cf,
                      svm_svmLinear_after_model = svm_svmLinear_after_model,
                      svm_svmLinear_before_cf = svm_svmLinear_before_cf,
                      svm_svmLinear_before_model = svm_svmLinear_before_model,
                      svm_svmLinear_before_ci = svm_svmLinear_before_ci,
                      svm_svmRadial_after_cf = svm_svmRadial_after_cf,
                      svm_svmRadial_after_model = svm_svmRadial_after_model,
                      svm_svmRadial_before_cf = svm_svmRadial_before_cf,
                      svm_svmRadial_before_model = svm_svmRadial_before_model,
                      svm_svmRadial_before_ci = svm_svmRadial_before_ci,
                      savingName = savingName,
                      partitionPercentage = partitionPercentage,
                      kernel = kernel,
                      enrichment = enrichment,
                      absolutePath = absolute_path,
                      omicData = omicData,
                      omic = omic,
                      domainTable = domainTable,
                      outliersSolution = outliersSolution,
                      thresholdControl = thresholdControl,
                      thresholdCase = thresholdCase,
                      mergeData = mergeData
                    ),
                    output_file = outputName,
                    output_dir = outputPath
  )
}

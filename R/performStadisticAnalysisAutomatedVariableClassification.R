#' performStadisticAnalysisUserVariableClassification
#'
#' This function performs ratio-based analysis on the active activePredictors passed as parameters.
#'
#' @param savingName String | Name under which the model and solution will be saved after execution. If the user does not set any name, it will create a string with the current date.
#'
#' @param changedClinicData Dataset | Dataset of clinic data that will be used.
#'
#' @param firstGroup String | String that contains the name of the class created for the firstGroup. Case -> Control
#' @param secondGroup String | String that contains the name of the class created for the secondGroup. Control -> Case
#'
#' @param activePredictors Array of Strings | Predictors on which the study of the ratios will be conducted after the genetic algorithm has been performed. Default value: All the predictors in clinic data, except classVariable and idColumn.
#' @param numericActivePredictors Array of Strings | Numerical predictors of the active predictor list. Default value: All predictors in the active predictor list that return TRUE in the is.numeric function or return TRUE in the is.integer function and have 7 or more distinct values.
#' @param categoricActivePredictors Array of Strings | Categoric predictors of the active predictor list. Default value: All predictors in the active predictor list that return FALSE in the is.numeric function or return TRUE in the is.integer function and have less than 7 distinct values.
#'
#' @param classVariable String | Target variable, which must be binary, meaning it has two possible values. If the user does not specify a path to his own data, the value for the sample data, Ca.Co.Last, will be used.
#'
#'
#'
#' @export
#'
#' @examples
#'
#' MLASDO::performStadisticAnalysisUserVariableClassification(savingName = savingName, changedClinicData = changedClinicData, firstGroup = firstGroup, secondGroup = secondGroup, activePredictors = activePredictors, categoricActivePredictors = categoricActivePredictors, numericActivePredictors = numericActivePredictors, classVariable = classVariable)


performStadisticAnalysisUserVariableClassification <- function(
    savingName,
    changedClinicData,
    firstGroup,
    secondGroup,
    activePredictors,
    categoricActivePredictors,
    numericActivePredictors,
    classVariable,
    oddPredictorsToIgnore
){

  #### DATA READING ####
  clinic <- changedClinicData

  #### GENETIC ALGORITHM SOLUTION READING ####

  #### DATA PROCESSING ####

  # Iterating through the categoric selected activePredictors for analysis
  for (activePredictor in categoricActivePredictors){

    # Checking if the activePredictor is categorical or discrete
    if(!is.numeric(clinic[[activePredictor]]) | (is.integer(clinic[[activePredictor]]) & length(unique(clinic[[activePredictor]])) < 7)){

      # Replacing NA with the string "Empty"
      clinic[[activePredictor]] <- replace(clinic[[activePredictor]], is.na(clinic[[activePredictor]]), "Empty")

    }

  }

  activePredictors <- activePredictors[!activePredictors %in% oddPredictorsToIgnore]


  #### PATIENT GROUP CREATION ####
  # Subset of Control's that remained as Control's
  cPrima <- subset(clinic,  clinic[[classVariable]] == "Control")

  # Subset of Case that changed to Control
  cPrimaPrima <- subset(clinic,  clinic[[classVariable]] == firstGroup)

  # Subset of Case's that remained as Case's
  ePrima <- subset(clinic, clinic[[classVariable]] == "Case")

  # Subset of Control that changed to Case
  ePrimaPrima <- subset(clinic,  clinic[[classVariable]] == secondGroup)

  # Total Control subset (before genetic)
  c <- rbind(cPrima, ePrimaPrima)

  # Total Case subset (before genetic)
  e <- rbind(ePrima, cPrimaPrima)

  ### STADISTICAL ANALYSIS ###
  totalOdds <- data.frame(OddRatio = numeric(0), LowerCI = numeric(0), UpperCI = numeric(0))

  totalWilcox <- data.frame(value = numeric(0), group = numeric(0), predictor = character(0))

  namesCategoricOdds <- character(0)
  namesNumericWilcox <- character(0)

  ### ODDS RATIO ###

  # Iterating through the selected activePredictors for analysis
  for(activePredictor in activePredictors){

    # Checking if the activePredictor is numerical
    if((is.numeric(clinic[[activePredictor]]) & !is.integer(clinic[[activePredictor]])) | (is.integer(clinic[[activePredictor]]) & length(unique(clinic[[activePredictor]])) >= 7)){

      # Obtaining the number of patients for the numerator and the denominator
      # Since the activePredictor is numeric, there will be empty columns
      NumPatientsCC <- nrow(ePrimaPrima)

      if(NumPatientsCC == 0){
        CC <- 0
      } else {
        CC <- ePrimaPrima[[activePredictor]]
      }

      NumPatientsCM <- nrow(cPrima)

      if(NumPatientsCM == 0){
        CM <- 0
      } else {
        CM <- cPrima[[activePredictor]]
      }


      NumPatientsEC <- nrow(cPrimaPrima)

      if(NumPatientsEC == 0){
        EC <- 0
      } else {
        EC <- cPrimaPrima[[activePredictor]]
      }


      NumPatientsEM <- nrow(ePrima)

      if(NumPatientsEM == 0){
        EM <- 0
      } else {
        EM <- ePrima[[activePredictor]]
      }


      NumPatientsC <- nrow(c)

      if(NumPatientsEM == 0){
        C <- 0
      } else {
        C <- c[[activePredictor]]
      }


      NumPatientsE <- nrow(e)

      if(NumPatientsE == 0){
        E <- 0
      } else {
        E <- e[[activePredictor]]
      }

      ### AÑADIMOS CONTROLS CAMBIADOS
      activePredictor_C <- paste0(activePredictor, "_Controls")

      new_rows <- data.frame(
        value = CC,
        group = "Changed",
        predictor = activePredictor_C
      )

      totalWilcox <- rbind(totalWilcox, new_rows)

      ### AÑADIMOS CONTROLES MANTENIDOS
      new_rows <- data.frame(
        value = CM,
        group = "Maintained",
        predictor = activePredictor_C
      )

      totalWilcox <- rbind(totalWilcox, new_rows)

      ###AÑADIMOS CASOS CAMBIADOS
      activePredictor_E <- paste0(activePredictor, "_Cases")

      new_rows <- data.frame(
        value = EC,
        group = "Changed",
        predictor = activePredictor_E
      )

      totalWilcox <- rbind(totalWilcox, new_rows)

      ### AÑADIMOS CONTROLES MANTENIDOS
      new_rows <- data.frame(
        value = EM,
        group = "Maintained",
        predictor = activePredictor_E
      )

      totalWilcox <- rbind(totalWilcox, new_rows)


    } else {

      factor <- as.factor(clinic[[activePredictor]])

      if (length(levels(factor)) > 2) {

        # We obtain a string formed by the two values separated by /

        ### ODDS RATIO CC ###

        clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

        clinic_odds$predictor <- as.factor(clinic[[activePredictor]])
        clinic_odds$odds <- ifelse(clinic[[classVariable]] == secondGroup, 1, 0)

        model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

        logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

        odd_ratio_table <- logistic_info$table

        defaultValue <-  sub(".*=([^=]+)$", "\\1", rownames(odd_ratio_table)[1])

        for (i in 2:(nrow(odd_ratio_table)-1)) {

          odd_ratio_str <- odd_ratio_table[i]

          combinationName <- paste(activePredictor, "/", sep = "")
          combinationName <- paste(combinationName, defaultValue, sep = "")
          combinationName <- paste(combinationName, trimws(rownames(odd_ratio_table)[i]), sep = " vs ")

          OddRatioCC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
          LowerCICC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
          UpperCICC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

          newRow <- data.frame(OddRatio = OddRatioCC, LowerCI = LowerCICC, UpperCI = UpperCICC)
          totalOdds <- rbind(totalOdds, newRow)
          combinationName_C <- paste(combinationName, "/Controls", sep = "")
          namesCategoricOdds <- c(namesCategoricOdds, combinationName_C)
        }

        ### ODDS RATIO CC ###


        ### ODDS RATIO EC ###

        clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

        clinic_odds$predictor <- as.factor(clinic[[activePredictor]])
        clinic_odds$odds <- ifelse(clinic[[classVariable]] == firstGroup, 1, 0)

        model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

        logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

        odd_ratio_table <- logistic_info$table

        defaultValue <-  sub(".*=([^=]+)$", "\\1", rownames(odd_ratio_table)[1])

        for (i in 2:(nrow(odd_ratio_table)-1)) {

          odd_ratio_str <- odd_ratio_table[i]

          combinationName <- paste(activePredictor, "/", sep = "")
          combinationName <- paste(combinationName, defaultValue, sep = "")
          combinationName <- paste(combinationName, trimws(rownames(odd_ratio_table)[i]), sep = " vs ")

          OddRatioEC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
          LowerCIEC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
          UpperCIEC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

          newRow <- data.frame(OddRatio = OddRatioEC, LowerCI = LowerCIEC, UpperCI = UpperCIEC)
          totalOdds <- rbind(totalOdds, newRow)
          combinationName_E <- paste(combinationName, "/Cases", sep="")
          namesCategoricOdds <- c(namesCategoricOdds, combinationName_E)

        }

      } else {

        # We obtain a string formed by the two values separated by /

        ### ODDS RATIO CC ###

        clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

        clinic_odds$predictor <- as.factor(clinic[[activePredictor]])
        clinic_odds$odds <- ifelse(clinic[[classVariable]] == secondGroup, 1, 0)

        model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

        logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

        odd_ratio_table <- logistic_info$table

        defaultValue <-  sub(".*=([^=]+)$", "\\1", rownames(odd_ratio_table)[1])

        combinationName <- paste(activePredictor, "/", sep = "")
        combinationName <- paste(combinationName, defaultValue, sep = "")

        odd_ratio_str <- odd_ratio_table[1]

        OddRatioCC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
        LowerCICC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
        UpperCICC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

        newRow <- data.frame(OddRatio = OddRatioCC, LowerCI = LowerCICC, UpperCI = UpperCICC)
        totalOdds <- rbind(totalOdds, newRow)
        combinationName_C <- paste(combinationName, "/Controls", sep = "")
        namesCategoricOdds <- c(namesCategoricOdds, combinationName_C)

        ### ODDS RATIO CC ###


        ### ODDS RATIO EC ###

        clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

        clinic_odds$predictor <- as.factor(clinic[[activePredictor]])
        clinic_odds$odds <- ifelse(clinic[[classVariable]] == firstGroup, 1, 0)

        model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

        logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

        odd_ratio_table <- logistic_info$table

        defaultValue <-  sub(".*=([^=]+)$", "\\1", rownames(odd_ratio_table)[1])

        combinationName <- paste(activePredictor, "/", sep = "")
        combinationName <- paste(combinationName, defaultValue, sep = "")

        odd_ratio_str <- odd_ratio_table[1]


        OddRatioEC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
        LowerCIEC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
        UpperCIEC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

        newRow <- data.frame(OddRatio = OddRatioEC, LowerCI = LowerCIEC, UpperCI = UpperCIEC)
        totalOdds <- rbind(totalOdds, newRow)
        combinationName_E <- paste(combinationName, "/Cases", sep="")
        namesCategoricOdds <- c(namesCategoricOdds, combinationName_E)

      }

    }

  }

  rownames(totalOdds) <- namesCategoricOdds
  totalOdds$y <- namesCategoricOdds

  na_rows <- totalOdds[is.na(totalOdds$UpperCI), ]

  if (nrow(na_rows) > 0) {
    cat("The following comparisons have obtained an infinite value in the upper limit of the Odds Ratio:\n")
    print(na_rows$y)
  }


  name <- paste("GA", savingName, sep="_")

  dirPath <- paste(savingName, "analysisData", name, sep = "/")

  # Save the stadistical analysis
  gaPathOdds <- paste(dirPath, "OddsRatios.tsv", sep="_")
  gaPathWilcox <- paste(dirPath, "Wilcox.tsv", sep="_")

  write.table(totalOdds, gaPathOdds, row.names = T, col.names = T, sep =  '\t')
  write.table(totalWilcox, gaPathWilcox, row.names = T, col.names = T, sep =  '\t')
}

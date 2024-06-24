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
    oddPredictorsToIgnore,
    oddsDefaultValues
){

  #### DATA READING ####
  clinic <- changedClinicData

  #### GENETIC ALGORITHM SOLUTION READING ####

  #### DATA PROCESSING ####

  # Iterating through the categoric selected activePredictors for analysis
  for (activePredictor in categoricActivePredictors){

    # Replacing NA with the string "Empty"
    clinic[[activePredictor]] <- replace(clinic[[activePredictor]], is.na(clinic[[activePredictor]]), "Empty")

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
  totalOddsControls <- data.frame(OddRatio = numeric(0), LowerCI = numeric(0), UpperCI = numeric(0), SumCCFirst = numeric(0), SumCCSecond = numeric(0), SumCMFirst = numeric(0), SumCMSecond = numeric(0))

  totalOddsCases <- data.frame(OddRatio = numeric(0), LowerCI = numeric(0), UpperCI = numeric(0), SumECFirst = numeric(0), SumECSecond = numeric(0), SumEMFirst = numeric(0), SumEMSecond = numeric(0))

  totalWilcoxControls <- data.frame(value = numeric(0), group = numeric(0), predictor = character(0))
  totalWilcoxCases <- data.frame(value = numeric(0), group = numeric(0), predictor = character(0))

  namesCategoricOddsControls <- character(0)
  namesCategoricOddsCases <- character(0)

  namesNumericWilcoxControls <- character(0)
  namesNumericWilcoxCases <- character(0)

  ### ODDS RATIO ###

  # Iterating through the selected activePredictors for analysis
  for(activePredictor in activePredictors){

    # Checking if the activePredictor is numerical
    if(activePredictor %in% numericActivePredictors){

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

      new_rows <- data.frame(
        value = CC,
        group = "Changed",
        predictor = activePredictor
      )

      totalWilcoxControls <- rbind(totalWilcoxControls, new_rows)

      ### AÑADIMOS CONTROLES MANTENIDOS
      new_rows <- data.frame(
        value = CM,
        group = "Maintained",
        predictor = activePredictor
      )

      totalWilcoxControls <- rbind(totalWilcoxControls, new_rows)

      ###AÑADIMOS CASOS CAMBIADOS

      new_rows <- data.frame(
        value = EC,
        group = "Changed",
        predictor = activePredictor
      )

      totalWilcoxCases <- rbind(totalWilcoxCases, new_rows)

      ### AÑADIMOS CONTROLES MANTENIDOS
      new_rows <- data.frame(
        value = EM,
        group = "Maintained",
        predictor = activePredictor
      )

      totalWilcoxCases <- rbind(totalWilcoxCases, new_rows)


    } else {

        factor <- as.factor(clinic[[activePredictor]])

        length_factor <- length(levels(factor))

        if (length_factor > 2) {

          #### ODDS RATIO FOR PREDICTORS WITH DEFAULT VALUE ####

          if(activePredictor %in% defaultOddsPredictors) {

            # We obtain a string formed by the two values separated by /

            ### ODDS RATIO CC ###

            clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

            clinic_odds$predictor <- as.factor(clinic[[activePredictor]])
            clinic_odds$odds <- ifelse(clinic[[classVariable]] == secondGroup, 1, 0)

            model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

            logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

            odd_ratio_table <- logistic_info$table

            defaultValue <-  sub(".*=([^=]+)$", "\\1", rownames(odd_ratio_table)[1])

            DCM <- sum(cPrima[[activePredictor]] == defaultValue)

            DCC <- sum(ePrimaPrima[[activePredictor]] == defaultValue)

            for (i in 2:(nrow(odd_ratio_table)-1)) {

              odd_ratio_str <- odd_ratio_table[i]

              actualValue <- trimws(rownames(odd_ratio_table)[i])

              combinationName <- paste(activePredictor, "/", sep = "")
              combinationName <- paste(combinationName, defaultValue, sep = "")
              combinationName <- paste(combinationName, actualValue, sep = " vs ")

              OddRatioCC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
              LowerCICC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
              UpperCICC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

              ACM <- sum(cPrima[[actualValue]] == defaultValue)

              ACC <- sum(ePrimaPrima[[actualValue]] == defaultValue)

              newRow <- data.frame(OddRatio = OddRatioCC, LowerCI = LowerCICC, UpperCI = UpperCICC, SUMCCFirst = DCC, SUMCCSecond = ACC, SUMCMFirst = DCM, SUMCMSecond = ACM)
              totalOddsControls <- rbind(totalOddsControls, newRow)
              namesCategoricOddsControls <- c(namesCategoricOddsControls, combinationName)
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

            DEM <- sum(ePrima[[activePredictor]] == defaultValue)

            DEC <- sum(cPrimaPrima[[activePredictor]] == defaultValue)

            for (i in 2:(nrow(odd_ratio_table)-1)) {

              odd_ratio_str <- odd_ratio_table[i]

              combinationName <- paste(activePredictor, "/", sep = "")
              combinationName <- paste(combinationName, defaultValue, sep = "")
              combinationName <- paste(combinationName, trimws(rownames(odd_ratio_table)[i]), sep = " vs ")

              OddRatioEC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
              LowerCIEC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
              UpperCIEC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

              AEM <- sum(ePrima[[actualValue]] == defaultValue)

              AEC <- sum(cPrimaPrima[[actualValue]] == defaultValue)

              newRow <- data.frame(OddRatio = OddRatioCC, LowerCI = LowerCICC, UpperCI = UpperCICC, SUMCCFirst = DEC, SUMCCSecond = AEC, SUMCMFirst = DEM, SUMCMSecond = AEM)
              totalOddsCases <- rbind(totalOddsCases, newRow)
              namesCategoricOddsCases <- c(namesCategoricOddsCases, combinationName)

            }

          }

          for(i in 1:length_factor){

              # We obtain a string formed by the two values separated by /

              ### ODDS RATIO CC ###

              level <- levels(factor)[i]

              clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

              clinic_odds$predictor <- as.factor(clinic[[activePredictor]])
              clinic_odds$predictor <- ifelse(clinic_odds$predictor == level, level, "Rest")
              clinic_odds$predictor <- as.factor(clinic_odds$predictor)

              clinic_odds$odds <- ifelse(clinic[[classVariable]] == secondGroup, 1, 0)

              model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

              logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

              odd_ratio_table <- logistic_info$table

              combinationName <- paste(activePredictor, "/", sep = "")
              combinationName <- paste(combinationName, level, sep = "")
              combinationName <- paste(combinationName, " vs Rest", sep = "")

              odd_ratio_str <- odd_ratio_table[1]

              OddRatioCC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
              LowerCICC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
              UpperCICC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

              DCM <- sum(cPrima[[activePredictor]] == level)
              DCC <- sum(ePrimaPrima[[activePredictor]] == level)

              RCM <- sum(cPrima[[activePredictor]] != level)
              RCC <- sum(ePrimaPrima[[activePredictor]] != level)

              newRow <- data.frame(OddRatio = OddRatioCC, LowerCI = LowerCICC, UpperCI = UpperCICC, SUMCCFirst = DCC, SUMCCSecond = RCC, SUMCMFirst = DCM, SUMCMSecond = RCM)
              totalOddsControls <- rbind(totalOddsControls, newRow)
              namesCategoricOddsControls <- c(namesCategoricOddsControls, combinationName)

              ### ODDS RATIO CC ###


              ### ODDS RATIO EC ###

              clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

              clinic_odds$predictor <- as.factor(clinic[[activePredictor]])
              clinic_odds$predictor <- ifelse(clinic_odds$predictor == level, level, "Rest")
              clinic_odds$predictor <- as.factor(clinic_odds$predictor)


              clinic_odds$odds <- ifelse(clinic[[classVariable]] == firstGroup, 1, 0)

              model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

              logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

              odd_ratio_table <- logistic_info$table


              combinationName <- paste(activePredictor, "/", sep = "")
              combinationName <- paste(combinationName, level, sep = "")
              combinationName <- paste(combinationName, " vs Rest", sep = "")

              odd_ratio_str <- odd_ratio_table[1]

              OddRatioEC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
              LowerCIEC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
              UpperCIEC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

              DEM <- sum(ePrima[[activePredictor]] == level)
              DEC <- sum(cPrimaPrima[[activePredictor]] == level)

              REM <- sum(ePrima[[activePredictor]] != level)
              REC <- sum(cPrimaPrima[[activePredictor]] != level)

              newRow <- data.frame(OddRatio = OddRatioCC, LowerCI = LowerCICC, UpperCI = UpperCICC, SUMCCFirst = DEC, SUMCCSecond = REC, SUMCMFirst = DEM, SUMCMSecond = REM)
              totalOddsCases <- rbind(totalOddsCases, newRow)
              namesCategoricOddsCases <- c(namesCategoricOddsCases, combinationName)

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

            defaultValue <-  gsub("predictor: ", "", rownames(odd_ratio_table)[1])

            combinationName <- paste(activePredictor, "/", sep = "")
            combinationName <- paste(combinationName, defaultValue, sep = "")

            odd_ratio_str <- odd_ratio_table[1]

            OddRatioCC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
            LowerCICC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
            UpperCICC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

            levels <- levels(clinic_odds$predictor)

            CM1 <- sum(cPrima[[activePredictor]] == levels[1])
            CC1 <- sum(ePrimaPrima[[activePredictor]] == levels[1])

            CM2 <- sum(cPrima[[activePredictor]] == levels[2])
            CC2 <- sum(ePrimaPrima[[activePredictor]] == levels[2])

            newRow <- data.frame(OddRatio = OddRatioCC, LowerCI = LowerCICC, UpperCI = UpperCICC, SUMCCFirst = CC1, SUMCCSecond = CC2, SUMCMFirst = CM1, SUMCMSecond = CM2)
            totalOddsControls <- rbind(totalOddsControls, newRow)
            namesCategoricOddsControls <- c(namesCategoricOddsControls, combinationName)

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

            levels <- levels(clinic_odds$predictor)

            EM1 <- sum(ePrima[[activePredictor]] == levels[1])
            EC1 <- sum(cPrimaPrima[[activePredictor]] == levels[1])

            EM2 <- sum(ePrima[[activePredictor]] == levels[2])
            EC2 <- sum(cPrimaPrima[[activePredictor]] == levels[2])

            newRow <- data.frame(OddRatio = OddRatioCC, LowerCI = LowerCICC, UpperCI = UpperCICC, SUMCCFirst = EC1, SUMCCSecond = EC2, SUMCMFirst = EM1, SUMCMSecond = EM2)
            totalOddsCases <- rbind(totalOddsCases, newRow)
            namesCategoricOddsCases <- c(namesCategoricOddsCases, combinationName)


        }

      }

  }

  rownames(totalOddsControls) <- namesCategoricOddsControls
  totalOddsControls$y <- namesCategoricOddsControls

  rownames(totalOddsCases) <- namesCategoricOddsCases
  totalOddsCases$y <- namesCategoricOddsCases

  na_rows_controls <- totalOddsControls[is.na(totalOddsControls$UpperCI), ]

  if (nrow(na_rows_controls) > 0) {
    cat("The following comparisons have obtained an infinite value in the upper limit of the Odds Ratio from Control to Cases patients:\n")
    print(na_rows_controls$y)
  }

  na_rows_cases <- totalOddsCases[is.na(totalOddsCases$UpperCI), ]

  if (nrow(na_rows_cases) > 0) {
    cat("The following comparisons have obtained an infinite value in the upper limit of the Odds Ratio from Cases to Controls patients:\n")
    print(na_rows_cases$y)
  }


  name <- paste("GA", savingName, sep="_")

  dirPath <- paste(savingName, "analysisData", name, sep = "/")

  # Save the stadistical analysis
  gaPathOddsControls <- paste(dirPath, "OddsRatiosControls.tsv", sep="_")
  gaPathWilcoxControls <- paste(dirPath, "WilcoxControls.tsv", sep="_")

  write.table(totalOddsControls, gaPathOddsControls, row.names = T, col.names = T, sep =  '\t')
  write.table(totalWilcoxControls, gaPathWilcoxControls, row.names = T, col.names = T, sep =  '\t')

  gaPathOddsCases <- paste(dirPath, "OddsRatiosCases.tsv", sep="_")
  gaPathWilcoxCases <- paste(dirPath, "WilcoxCases.tsv", sep="_")

  write.table(totalOddsCases, gaPathOddsCases, row.names = T, col.names = T, sep =  '\t')
  write.table(totalWilcoxCases, gaPathWilcoxCases, row.names = T, col.names = T, sep =  '\t')
}

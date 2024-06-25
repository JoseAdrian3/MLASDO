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
#' @param defaultOddsPredictors String | Vector of names of categorical features that have a default value to compare this value against the rest of the values of that feature.
#'
#'
#' @export
#'
#' @examples
#'
#' MLASDO::performStadisticAnalysisAutomatedVariableClassification(savingName = savingName, changedClinicData = changedClinicData, firstGroup = firstGroup, secondGroup = secondGroup, activePredictors = activePredictors, categoricActivePredictors = categoricActivePredictors, numericActivePredictors = numericActivePredictors, classVariable = classVariable, defaultOddsPredictors = defaultOddsPredictors)


performStadisticAnalysisAutomatedVariableClassification <- function(
    savingName,
    changedClinicData,
    firstGroup,
    secondGroup,
    activePredictors,
    categoricActivePredictors,
    numericActivePredictors,
    classVariable,
    defaultOddsPredictors
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

  # Variables to store the information and names of the odds ratios between the categorical variables and the patients who have transitioned from Control to Case
  totalOddsControls <- data.frame(OddRatio = numeric(0), LowerCI = numeric(0), UpperCI = numeric(0), SumCCFirst = numeric(0), SumCCSecond = numeric(0), SumCMFirst = numeric(0), SumCMSecond = numeric(0))
  namesCategoricOddsControls <- character(0)

  # Variables to store the information and names of the odds ratios between the categorical variables and the patients who have transitioned from Case to Control
  totalOddsCases <- data.frame(OddRatio = numeric(0), LowerCI = numeric(0), UpperCI = numeric(0), SumECFirst = numeric(0), SumECSecond = numeric(0), SumEMFirst = numeric(0), SumEMSecond = numeric(0))
  namesCategoricOddsCases <- character(0)

  # Variables to store the information and names on the odds ratios between the categorical variables and the patients who have transitioned from Control to Case
  totalWilcoxControls <- data.frame(value = numeric(0), group = numeric(0), predictor = character(0))
  namesNumericWilcoxControls <- character(0)

  # Variables to store the information and names on the odds ratios between the categorical variables and the patients who have transitioned from Case to Control
  totalWilcoxCases <- data.frame(value = numeric(0), group = numeric(0), predictor = character(0))
  namesNumericWilcoxCases <- character(0)

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

      #### ADD CONTROLS CHANGED ####

      new_rows <- data.frame(
        value = CC,
        group = "Changed",
        predictor = activePredictor
      )

      totalWilcoxControls <- rbind(totalWilcoxControls, new_rows)

      #### ADD CONTROLS CHANGED ####

      #### ADD CONTROLS MAINTAINED ####

      new_rows <- data.frame(
        value = CM,
        group = "Maintained",
        predictor = activePredictor
      )

      totalWilcoxControls <- rbind(totalWilcoxControls, new_rows)

      #### ADD CONTROLS MAINTAINED ####

      #### ADD CASES CHANGED ####

      new_rows <- data.frame(
        value = EC,
        group = "Changed",
        predictor = activePredictor
      )

      totalWilcoxCases <- rbind(totalWilcoxCases, new_rows)

      #### ADD CASES CHANGED ####

      #### ADD CASES MAINTAINED ####

      new_rows <- data.frame(
        value = EM,
        group = "Maintained",
        predictor = activePredictor
      )

      totalWilcoxCases <- rbind(totalWilcoxCases, new_rows)

      #### ADD CASES MAINTAINED ####

      # If not, is a categorical feature
    } else {

      # Getting levels of the factor of the categorical feature
      factor <- as.factor(clinic[[activePredictor]])
      length_factor <- length(levels(factor))

      # Checking if it has more than 2 levels
      if (length_factor > 2) {

        #### ODDS RATIO FOR PREDICTORS WITH DEFAULT VALUE ####

        # Checking if the feature has default value
        if(activePredictor %in% defaultOddsPredictors) {

          ### DEFAULT ODDS RATIO CONTROLS ###

          # Dataframe for storing the odds information for the glm model
          clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

          # Active predictor (factor already)
          clinic_odds$predictor <- clinic[[activePredictor]]

          # Vector with the patients changed from Controls to Cases of the feature
          clinic_odds$odds <- ifelse(clinic[[classVariable]] == secondGroup, 1, 0)

          # glm predicting changes using the predictor
          model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

          # "logistic.display2" for the odds ratios of the default value (first level) versus the others values
          logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

          # Logistic display in case of having more than 2 levels returns a table where, form row 2, each row contains the odds ratio of
          # comparing the default value with each of the other values
          odd_ratio_table <- logistic_info$table

          # In the row 1 is the text of the default value
          defaultValue <-  gsub("predictor: ref.=", "", rownames(odd_ratio_table)[1])

          # Number of patients maintained in Controls of the default value
          DCM <- sum(cPrima[[activePredictor]] == defaultValue)

          # Number of patients changed from Controls to Cases of the default value
          DCC <- sum(ePrimaPrima[[activePredictor]] == defaultValue)

          # iterating from every odds ratio value from the second row to the penultimate (the table has one last empty row)
          for (i in 2:(nrow(odd_ratio_table)-1)) {

            # Process to obtain the value with which the default value is compared
            odd_ratio_str <- odd_ratio_table[i]
            actualValue <- trimws(rownames(odd_ratio_table)[i])

            # Process to obtain the text "default value vs actual value"
            combinationName <- paste(activePredictor, "/", sep = "")
            combinationName <- paste(combinationName, defaultValue, sep = "")
            combinationName <- paste(combinationName, actualValue, sep = " vs ")

            # Process to obtain the numeric values of the odd ratio, upper ci and lower ci
            OddRatioCC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
            LowerCICC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
            UpperCICC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

            # Number of patients maintained in Controls of the actual value
            ACM <- sum(cPrima[[actualValue]] == defaultValue)

            # Number of patients changed from Controls to Cases of the actual value
            ACC <- sum(ePrimaPrima[[actualValue]] == defaultValue)

            # Updeting variable with the odds information and names of the Controls patients
            newRow <- data.frame(OddRatio = OddRatioCC, LowerCI = LowerCICC, UpperCI = UpperCICC, SUMCCFirst = DCC, SUMCCSecond = ACC, SUMCMFirst = DCM, SUMCMSecond = ACM)
            totalOddsControls <- rbind(totalOddsControls, newRow)
            namesCategoricOddsControls <- c(namesCategoricOddsControls, combinationName)
          }

          ### DEFAULT ODDS RATIO CONTROLS ###

          ### DEFAULT ODDS RATIO CASES ###

          # Dataframe for storing the odds information for the glm model
          clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

          # Active predictor (factor already)
          clinic_odds$predictor <- clinic[[activePredictor]]

          # Vector with the patients changed from Cases to Controls of the feature
          clinic_odds$odds <- ifelse(clinic[[classVariable]] == firstGroup, 1, 0)

          # glm predicting changes using the predictor
          model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

          # "logistic.display2" for the odds ratios of the default value (first level) versus the others values
          logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

          # Logistic display in case of having more than 2 levels returns a table where, form row 2, each row contains the odds ratio of
          # comparing the default value with each of the other values
          odd_ratio_table <- logistic_info$table

          # In the row 1 is the text of the default value
          defaultValue <-  gsub("predictor: ref.=", "", rownames(odd_ratio_table)[1])

          # Number of patients maintained in Cases of the default value
          DEM <- sum(ePrima[[activePredictor]] == defaultValue)

          # Number of patients changed from Cases to Controls of the default value
          DEC <- sum(cPrimaPrima[[activePredictor]] == defaultValue)

          # iterating from every odds ratio value from the second row to the penultimate (the table has one last empty row)
          for (i in 2:(nrow(odd_ratio_table)-1)) {

            # Process to obtain the value with which the default value is compared
            odd_ratio_str <- odd_ratio_table[i]

            # Process to obtain the text "default value vs actual value"
            combinationName <- paste(activePredictor, "/", sep = "")
            combinationName <- paste(combinationName, defaultValue, sep = "")
            combinationName <- paste(combinationName, trimws(rownames(odd_ratio_table)[i]), sep = " vs ")

            # Process to obtain the numeric values of the odd ratio, upper ci and lower ci
            OddRatioEC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
            LowerCIEC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
            UpperCIEC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

            # Number of patients maintained in Cases of the actual value
            AEM <- sum(ePrima[[actualValue]] == defaultValue)

            # Number of patients changed from Cases to Controls of the actual value
            AEC <- sum(cPrimaPrima[[actualValue]] == defaultValue)

            # Updeting variable with the odds information and names of the Cases patients
            newRow <- data.frame(OddRatio = OddRatioEC, LowerCI = LowerCIEC, UpperCI = UpperCIEC, SUMCCFirst = DEC, SUMCCSecond = AEC, SUMCMFirst = DEM, SUMCMSecond = AEM)
            totalOddsCases <- rbind(totalOddsCases, newRow)
            namesCategoricOddsCases <- c(namesCategoricOddsCases, combinationName)

          }

          ### DEFAULT ODDS RATIO CASES ###

        }

        # Whether a feature has a default value or not, the odds ratio of each value against the rest will be obtained
        for(i in 1:length_factor){

          ### DEFAULT ODDS RATIO CONTROLS ###

          # The current level is obtained with which it will be compared against the rest to obtain the odds ratio
          level <- levels(factor)[i]

          # Dataframe for storing the odds information for the glm model
          clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

          # Active predictor as a factor
          clinic_odds$predictor <- as.factor(clinic[[activePredictor]])

          # Predictor as a factor with the current value and in the remaining positions with others values the value "Rest"
          clinic_odds$predictor <- ifelse(clinic_odds$predictor == level, level, "Rest")
          clinic_odds$predictor <- as.factor(clinic_odds$predictor)

          # Vector with the patients changed from Controls to Cases of the feature
          clinic_odds$odds <- ifelse(clinic[[classVariable]] == secondGroup, 1, 0)

          # glm predicting changes using the predictor
          model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

          # "logistic.display2" for the odd ratio between the actual value versus the others
          logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

          # Logistic display in case of having 2 levels returns a table where in the row 1 contains the odd ratio
          odd_ratio_table <- logistic_info$table

          # Process to obtain the text "actual value vs Rest"
          combinationName <- paste(activePredictor, "/", sep = "")
          combinationName <- paste(combinationName, level, sep = "")
          combinationName <- paste(combinationName, " vs Rest", sep = "")

          # Process to obtain the numeric values of the odd ratio, upper ci and lower ci
          odd_ratio_str <- odd_ratio_table[1]
          OddRatioCC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
          LowerCICC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
          UpperCICC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

          # Number of patients maintained in Controls of the actual value and the rest
          DCM <- sum(cPrima[[activePredictor]] == level)
          RCM <- sum(cPrima[[activePredictor]] != level)

          # Number of patients changed from Controls to Cases of the actual value
          DCC <- sum(ePrimaPrima[[activePredictor]] == level)
          RCC <- sum(ePrimaPrima[[activePredictor]] != level)

          # Updeting variable with the odds information and names of the Controls patients
          newRow <- data.frame(OddRatio = OddRatioCC, LowerCI = LowerCICC, UpperCI = UpperCICC, SUMCCFirst = DCC, SUMCCSecond = RCC, SUMCMFirst = DCM, SUMCMSecond = RCM)
          totalOddsControls <- rbind(totalOddsControls, newRow)
          namesCategoricOddsControls <- c(namesCategoricOddsControls, combinationName)

          ### DEFAULT ODDS RATIO CONTROLS ###

          ### DEFAULT ODDS RATIO CASES ###

          # Dataframe for storing the odds information for the glm model
          clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

          # Active predictor as a factor
          clinic_odds$predictor <- as.factor(clinic[[activePredictor]])

          # Predictor as a factor with the current value and in the remaining positions with others values the value "Rest"
          clinic_odds$predictor <- ifelse(clinic_odds$predictor == level, level, "Rest")
          clinic_odds$predictor <- as.factor(clinic_odds$predictor)

          # Vector with the patients changed from Cases to Controls of the feature
          clinic_odds$odds <- ifelse(clinic[[classVariable]] == firstGroup, 1, 0)

          # glm predicting changes using the predictor
          model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

          # "logistic.display2" for the odd ratio between the actual value versus the others
          logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

          # Logistic display in case of having 2 levels returns a table where in the row 1 contains the odd ratio
          odd_ratio_table <- logistic_info$table

          # Process to obtain the text "actual value vs Rest"
          combinationName <- paste(activePredictor, "/", sep = "")
          combinationName <- paste(combinationName, level, sep = "")
          combinationName <- paste(combinationName, " vs Rest", sep = "")

          # Process to obtain the numeric values of the odd ratio, upper ci and lower ci
          odd_ratio_str <- odd_ratio_table[1]
          OddRatioEC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
          LowerCIEC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
          UpperCIEC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

          # Number of patients maintained in Cases of the actual value and the rest
          DEM <- sum(ePrima[[activePredictor]] == level)
          DEC <- sum(cPrimaPrima[[activePredictor]] == level)

          # Number of patients changed from Cases to Controls of the actual value
          REM <- sum(ePrima[[activePredictor]] != level)
          REC <- sum(cPrimaPrima[[activePredictor]] != level)

          # Updeting variable with the odds information and names of the Cases patients
          newRow <- data.frame(OddRatio = OddRatioEC, LowerCI = LowerCIEC, UpperCI = UpperCIEC, SUMCCFirst = DEC, SUMCCSecond = REC, SUMCMFirst = DEM, SUMCMSecond = REM)
          totalOddsCases <- rbind(totalOddsCases, newRow)
          namesCategoricOddsCases <- c(namesCategoricOddsCases, combinationName)

          ### DEFAULT ODDS RATIO CASES ###
        }

        # if not, it has only 2 levels
      } else {

        ### DEFAULT ODDS RATIO CONTROLS ###

        # Dataframe for storing the odds information for the glm model
        clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

        # Active predictor as a factor
        clinic_odds$predictor <- as.factor(clinic[[activePredictor]])

        # Vector with the patients changed from Controls to Cases of the feature
        clinic_odds$odds <- ifelse(clinic[[classVariable]] == secondGroup, 1, 0)

        # glm predicting changes using the predictor
        model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

        # "logistic.display2" for the odd ratio between the actual value versus the others
        logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

        # Logistic display in case of having 2 levels returns a table where in the row 1 contains the odd ratio
        odd_ratio_table <- logistic_info$table

        # Getting the text of "Value 1 vs Value 2"
        defaultValue <-  gsub("predictor: ", "", rownames(odd_ratio_table)[1])

        # Getting the text of "feature/Value 1 vs Value 2
        combinationName <- paste(activePredictor, "/", sep = "")
        combinationName <- paste(combinationName, defaultValue, sep = "")

        # Process to obtain the numeric values of the odd ratio, upper ci and lower ci
        odd_ratio_str <- odd_ratio_table[1]
        OddRatioCC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
        LowerCICC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
        UpperCICC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

        # Getting the two levels
        levels <- levels(clinic_odds$predictor)

        # Number of patients maintained in Controls of the value 1 and value 2
        CM1 <- sum(cPrima[[activePredictor]] == levels[1])
        CM2 <- sum(cPrima[[activePredictor]] == levels[2])

        # Number of patients changed from Controls to Cases of the value 1 and value 2
        CC1 <- sum(ePrimaPrima[[activePredictor]] == levels[1])
        CC2 <- sum(ePrimaPrima[[activePredictor]] == levels[2])

        # Updeting variable with the odds information and names of the Controls patients
        newRow <- data.frame(OddRatio = OddRatioCC, LowerCI = LowerCICC, UpperCI = UpperCICC, SUMCCFirst = CC1, SUMCCSecond = CC2, SUMCMFirst = CM1, SUMCMSecond = CM2)
        totalOddsControls <- rbind(totalOddsControls, newRow)
        namesCategoricOddsControls <- c(namesCategoricOddsControls, combinationName)

        ### DEFAULT ODDS RATIO CONTROLS ###

        ### DEFAULT ODDS RATIO CASES ###

        # Dataframe for storing the odds information for the glm model
        clinic_odds <- data.frame(predictor = numeric(nrow(clinic)), odds = numeric(nrow(clinic)))

        # Active predictor as a factor
        clinic_odds$predictor <- as.factor(clinic[[activePredictor]])

        # Vector with the patients changed from Cases to Controls of the feature
        clinic_odds$odds <- ifelse(clinic[[classVariable]] == firstGroup, 1, 0)

        # glm predicting changes using the predictor
        model <- glm(odds ~ predictor, data = clinic_odds, family = binomial)

        # "logistic.display2" for the odd ratio between the actual value versus the others
        logistic_info <- logistic.display2(model, crude = TRUE, crude.p.value = TRUE, decimal = 3)

        # Logistic display in case of having 2 levels returns a table where in the row 1 contains the odd ratio
        odd_ratio_table <- logistic_info$table

        # Getting the text of "Value 1 vs Value 2"
        defaultValue <-  gsub("predictor: ", "", rownames(odd_ratio_table)[1])

        # Getting the text of "feature/Value 1 vs Value 2
        combinationName <- paste(activePredictor, "/", sep = "")
        combinationName <- paste(combinationName, defaultValue, sep = "")

        # Process to obtain the numeric values of the odd ratio, upper ci and lower ci
        odd_ratio_str <- odd_ratio_table[1]
        OddRatioEC <- as.numeric(sub("\\s*\\(.*", "", odd_ratio_str))
        LowerCIEC <- as.numeric(sub(".*\\((.*),.*", "\\1", odd_ratio_str))
        UpperCIEC <- as.numeric(gsub("^.*\\(\\d+\\.\\d+,\\s*(\\d+\\.\\d+).*", "\\1", odd_ratio_str))

        # Getting the two levels
        levels <- levels(clinic_odds$predictor)

        # Number of patients maintained in Cases of the value 1 and value 2
        EM1 <- sum(ePrima[[activePredictor]] == levels[1])
        EM2 <- sum(ePrima[[activePredictor]] == levels[2])

        # Number of patients changed from Cases to Controls of the value 1 and value 2
        EC1 <- sum(cPrimaPrima[[activePredictor]] == levels[1])
        EC2 <- sum(cPrimaPrima[[activePredictor]] == levels[2])

        # Updeting variable with the odds information and names of the Cases patients
        newRow <- data.frame(OddRatio = OddRatioEC, LowerCI = LowerCIEC, UpperCI = UpperCIEC, SUMCCFirst = EC1, SUMCCSecond = EC2, SUMCMFirst = EM1, SUMCMSecond = EM2)
        totalOddsCases <- rbind(totalOddsCases, newRow)
        namesCategoricOddsCases <- c(namesCategoricOddsCases, combinationName)

        ### DEFAULT ODDS RATIO CASES ###

      }

    }

  }

  # Adding the Odd Ratio names
  rownames(totalOddsControls) <- namesCategoricOddsControls
  rownames(totalOddsCases) <- namesCategoricOddsCases

  # Printing NA Controls to Cases Odds Ratios
  na_rows_controls <- totalOddsControls[is.na(totalOddsControls$UpperCI), ]

  if (nrow(na_rows_controls) > 0) {
    cat("The following comparisons have obtained an infinite value in the upper limit of the Odds Ratio from Control to Cases patients:\n")
    print(rownames(na_rows_controls))
  }

  # Removing Controls to Cases Odds Ratios with NA values
  totalOddsControls <- totalOddsControls[!is.na(totalOddsControls$UpperCI), ]

  # Printing NA Cases to Controls Odds Ratios
  na_rows_cases <- totalOddsCases[is.na(totalOddsCases$UpperCI), ]

  if (nrow(na_rows_cases) > 0) {
    cat("The following comparisons have obtained an infinite value in the upper limit of the Odds Ratio from Cases to Controls patients:\n")
    print(rownames(na_rows_cases))
  }

  # Removing Cases to Controls Odds Ratios with NA values
  totalOddsCases <- totalOddsCases[!is.na(totalOddsCases$UpperCI), ]

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

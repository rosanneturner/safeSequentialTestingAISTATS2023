#' Perform one stratified 2x2 contingency table confidence sequence simulation.
#' @author R.J. Turner
#'
#' @param realNumberOfDataBlocks vector with number of data blocks to simulate for.
#' @param realProbs list of the probabilities and batch sizes to simulate for, for structure
#' see the example below.
#' @param sharing cross-talk type, one of  c("logOdds", "nuisance", "none", "mixNuisance").
#' @param alpha desired type-I error guarantee, numeric.
#' @param numberForSeed seed to set for the simulation.
#' @param returnRawEValues boolean indicating to return E-values for individual data blocks
#' @param theta0GreaterThanDelta evaluate null hypothesis of the form \Theta_0 \geq \delta
#' @param theta0SmallerThanDelta evaluate null hypothesis of the form \Theta_0 \leq \delta
#' @param CIGridStart start of the CI grid, numeric
#' @param CIGridStop end of the CI grid, numeric
#' @param CIGridPrecision precision of the CI grid, numeric
#' @param pretendAllSameStratum boolean indicating to ignore strata
#'
#' @return CIResults, a tibble with E values or rejection results per delta for each stratum.
#'
#' @examples
#' alpha <- 0.05
#' globalSeed <- 1842781
#' realNumberOfDataBlocks <- c(1, 1, 1) * 5
#' realProbs <- list(
#'  c(
#'    "pa" = 0.1,
#'    "pb" = 0.15,
#'    "na" = 1,
#'    "nb" = 1
#'  ),
#'  c(
#'    "pa" = 0.2,
#'    "pb" = 0.6,
#'    "na" = 1,
#'    "nb" = 1
#'  ),
#'  c(
#'    "pa" = 0.8,
#'    "pb" = 0.2,
#'    "na" = 1,
#'    "nb" = 1
#'  )
#' )
#'
#' CIResult <- performStratified2x2CISimulation(
#'  realNumberOfDataBlocks = realNumberOfDataBlocks,
#'  realProbs = realProbs,
#'  sharing = "none",
#'  alpha = alpha,
#'  numberForSeed = globalSeed
#' )
performStratified2x2CISimulation <-
  function(realNumberOfDataBlocks,
           realProbs,
           sharing = c("logOdds", "nuisance", "none", "mixNuisance"),
           alpha,
           numberForSeed,
           returnRawEValues = FALSE,
           theta0GreaterThanDelta = FALSE,
           theta0SmallerThanDelta = FALSE,
           CIGridStart = -0.99,
           CIGridStop = 0.99,
           CIGridPrecision = 0.01,
           pretendAllSameStratum = FALSE) {
    sharing <- match.arg(sharing)
    K <-
      ifelse(pretendAllSameStratum, 1, length(realNumberOfDataBlocks))
    
    set.seed(numberForSeed)
    observationsAndEstimates <-
      createStratified2x2CIData(
        realNumberOfDataBlocks = realNumberOfDataBlocks,
        realProbs = realProbs,
        pretendAllSameStratum = pretendAllSameStratum
      )
    
    CIResults <-
      calculateEPerStratumPerDelta(
        observationsAndEstimates = observationsAndEstimates,
        sharing = sharing,
        CIGridStart = CIGridStart,
        CIGridStop = CIGridStop,
        CIGridPrecision = CIGridPrecision,
        theta0GreaterThanDelta = theta0GreaterThanDelta,
        theta0SmallerThanDelta = theta0SmallerThanDelta,
        returnRawEValues = returnRawEValues,
        alpha = alpha
      )
    
    return(CIResults)
  }

#' Simulate one sequence of stratified 2x2 contingency table data and calculate estimators used
#' for the calculations of E-values.
#' @author R.J. Turner
#'
#' @param realNumberOfDataBlocks vector with number of data blocks to simulate for.
#' @param realProbs list of the probabilities and batch sizes to simulate for, for structure
#' see the example below.
#' @param pretendAllSameStratum boolean indicating to ignore strata
#'
#' @return a tibble with observations for each batch, and estimates based on data up to and
#' including the previous batch.
#'
#' @examples
#' realNumberOfDataBlocks <- c(1, 1, 1) * 5
#' realProbs <- list(
#'  c(
#'    "pa" = 0.1,
#'    "pb" = 0.15,
#'    "na" = 1,
#'    "nb" = 1
#'  ),
#'  c(
#'    "pa" = 0.2,
#'    "pb" = 0.6,
#'    "na" = 1,
#'    "nb" = 1
#'  ),
#'  c(
#'    "pa" = 0.8,
#'    "pb" = 0.2,
#'    "na" = 1,
#'    "nb" = 1
#'  )
#' )
#'
#' set.seed(globalSeed)
#' simulatedData <- createStratified2x2CIData(realNumberOfDataBlocks, realProbs)
createStratified2x2CIData <-
  function(realNumberOfDataBlocks,
           realProbs,
           pretendAllSameStratum = FALSE) {
    K <- length(realNumberOfDataBlocks)
    observedData <- data.frame()
    
    for (k in seq_along(realProbs)) {
      ya <- rbinom(n = realNumberOfDataBlocks[k],
                   size = realProbs[[k]][["na"]],
                   prob = realProbs[[k]][["pa"]])
      yb <- rbinom(n = realNumberOfDataBlocks[k],
                   size = realProbs[[k]][["nb"]],
                   prob = realProbs[[k]][["pb"]])
      observedData <-
        rbind(
          observedData,
          data.frame(
            ya = ya,
            yb = yb,
            na = realProbs[[k]][["na"]],
            nb = realProbs[[k]][["nb"]],
            stratum = k
          )
        )
    }
    
    #we ignore the strata: pretend everything is from the same stratum
    if (pretendAllSameStratum) {
      observedData <- observedData %>%
        mutate(stratum = 1)
    }
    
    #randomly determine the order of complete data blocks in the groups
    observedData %<>%
      mutate(mTotal = sample(
        1:nrow(observedData),
        size = nrow(observedData),
        replace = FALSE
      )) %>%
      arrange(mTotal)
    
    #retrieve total and stratified counts and odds estiamtes per time point
    observationsAndEstimates <- observedData %>%
      retrieveTotalAndStratifiedCountsAndOddsEstimatesPerTimePoint()
    
    return(observationsAndEstimates)
  }

#' For a tibble of simulated stratified 2x2 data, calculate E values for each batch, per stratum
#' per delta in the CI grid.
#' @author R.J. Turner
#'
#' @inheritParams performStratified2x2CISimulation
#'
#' @return CIResults, a tibble with E values or rejection results per delta for each stratum.
#'
#' @examples
#' alpha <- 0.05
#' globalSeed <- 1842781
#' realNumberOfDataBlocks <- c(1, 1, 1) * 5
#' realProbs <- list(
#'  c(
#'    "pa" = 0.1,
#'    "pb" = 0.15,
#'    "na" = 1,
#'    "nb" = 1
#'  ),
#'  c(
#'    "pa" = 0.2,
#'    "pb" = 0.6,
#'    "na" = 1,
#'    "nb" = 1
#'  ),
#'  c(
#'    "pa" = 0.8,
#'    "pb" = 0.2,
#'    "na" = 1,
#'    "nb" = 1
#'  )
#' )
#'
#' set.seed(globalSeed)
#' simulatedData <- createStratified2x2CIData(realNumberOfDataBlocks, realProbs)
#' CIResult <- calculateEPerStratumPerDelta(simulatedData, sharing = "none", alpha = alpha)
calculateEPerStratumPerDelta <- function(observationsAndEstimates, sharing = "none",
                                         CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.01,
                                         theta0GreaterThanDelta = FALSE, theta0SmallerThanDelta = FALSE,
                                         returnRawEValues = FALSE, alpha = 0.05) {
  #calculate E-values for each stratum, for each delta by looping first over strata,
  #and therewithin over delta
  # loop over stratum 
  CIResults <- data.frame()
  
  if(sharing == "mixNuisance"){
    stratumEstimatesAllStrata <- createStratumEstimates(stratumData = observationsAndEstimates, sharing = sharing)
  }
  
  for (k in unique(observationsAndEstimates$stratum)) {
    
    if(sharing != "mixNuisance"){
      stratumData <- observationsAndEstimates %>%
        dplyr::filter(stratum == k)
      
      #estimates breveThetaA and breveThetaB are the same, whichever H0 we test against.
      #can do this before we loop over delta
      stratumEstimates <-
        createStratumEstimates(stratumData = stratumData, sharing = sharing)
    } else {
      stratumEstimates <- stratumEstimatesAllStrata %>% dplyr::filter(stratum == k)
    }
    
    #loop over H0 delta, determine RIPr and calculate E
    deltaGrid <- seq(CIGridStart, CIGridStop, by = CIGridPrecision)
    for (delta in deltaGrid) {
      #if min difference estimation, only need to find RIPr if breveThetaB - breveThetaA < delta
      eValueResultsForDelta <- stratumEstimates %>%
        mutate(
          thetaARIPr = case_when(
            is.na(breveThetaA) ~ NA_real_,
            #H1 coincides with H0 in the two following cases:
            theta0GreaterThanDelta & (breveThetaB - breveThetaA >= delta) ~ breveThetaA,
            theta0SmallerThanDelta & (breveThetaB - breveThetaA <= delta) ~ breveThetaA,
            TRUE ~ vRetrieveThetaARIPrRiskDifference(
              delta = delta,
              na = na,
              nb = nb,
              breveThetaA = breveThetaA,
              breveThetaB = breveThetaB
            )
          ),
          thetaBRIPr = case_when(
            theta0GreaterThanDelta & (breveThetaB - breveThetaA >= delta) ~ breveThetaB,
            theta0SmallerThanDelta & (breveThetaB - breveThetaA <= delta) ~ breveThetaB,
            TRUE ~ thetaARIPr + delta
          ),
          likelihoodRIPr = likelihoodTwoProportions(
            na1 = ya,
            na = na,
            nb1 = yb,
            nb = nb,
            thetaA = thetaARIPr,
            thetaB = thetaBRIPr
          ),
          mE = case_when(
            is.na(breveThetaA) | is.na(thetaARIPr) | is.na(likelihoodAlternative) ~ 1,
            TRUE ~ likelihoodAlternative / likelihoodRIPr
          ),
          E = cumprod(mE),
          delta = delta
        )
      
      
      if(returnRawEValues) {
        CIResults <- rbind(
          CIResults,
          eValueResultsForDelta %>%
            select(
              stratum, delta, mTotal, mStratum, mE, E
            )
        )
      }else {
        CIResults <- rbind(
          CIResults,
          eValueResultsForDelta %>%
            summarise(
              stratum = k,
              delta  = unique(delta),
              reject = any(E >= (1 / alpha)),
              mRejectStratum = ifelse(reject, min(which(E >= 1 / alpha)), NA_integer_),
              mRejectTotal = ifelse(reject, mTotal[min(which(E >= 1 / alpha))], NA_integer_)
            )
        )
      }
    }
  }
  return(CIResults)
}

#' Create the point estimates for calculating E-values 
#' @author R.J. Turner
#'
#' @param stratumData tibble with data for one stratum, generated with createStratified2x2CIData
#' @param sharing 
#' @param controlSuccessPrior hyperparameter for calculating the E-values, by default 0.18 (optimal)
#' in empirical experiments, see Turner, Ly and Grunwald in arXiv TODO). 
#' @param controlFailPrior hyperparameter for calculating the E-values, see controlSuccessPrior
#' @param interventionSuccessPrior hyperparameter for calculating the E-values, see controlSuccessPrior
#' @param interventionFailPrior hyperparameter for calculating the E-values, see controlSuccessPrior
#' @param gridSize grid size for the point estimate of \theta_a and \theta_b in the alternative
#' of the E-value when placing a restriction on the parameter space in terms of the log odds ratio.
#'
#' @return tibble enriched with estimates for calculating the E values per batch, per stratum, per delta.
#'
#' @examples
#' realNumberOfDataBlocks <- c(1, 1, 1) * 5
#' realProbs <- list(
#'  c(
#'    "pa" = 0.1,
#'    "pb" = 0.15,
#'    "na" = 1,
#'    "nb" = 1
#'  ),
#'  c(
#'    "pa" = 0.2,
#'    "pb" = 0.6,
#'    "na" = 1,
#'    "nb" = 1
#'  ),
#'  c(
#'    "pa" = 0.8,
#'    "pb" = 0.2,
#'    "na" = 1,
#'    "nb" = 1
#'  )
#' )
#'
#' set.seed(globalSeed)
#' simulatedData <- createStratified2x2CIData(realNumberOfDataBlocks, realProbs)
#' stratumData <- simulatedData %>% dplyr::filter(stratum == 1)
#' stratumEstimates <- createStratumEstimates(stratumData, sharing = "none")
createStratumEstimates <- function(stratumData,
                                   sharing = c("logOdds", "nuisance", "none", "riskDifference", "mixNuisance"),
                                   controlSuccessPrior = 0.18,
                                   controlFailPrior = 0.18,
                                   interventionSuccessPrior = 0.18,
                                   interventionFailPrior = 0.18,
                                   gridSize = 2000) {
  sharing <- match.arg(sharing)
  
  if (sharing == "logOdds") {
    stratumEstimates <- stratumData %>%
      mutate(
        breveThetaA = vCalculateOddsRatioBasedMarginalPredProbThetaA(
          gridSize = gridSize,
          betaA1 = controlSuccessPrior,
          betaA2 = controlFailPrior,
          logOddsEstimate = log(previousOddsEstimateTotal),
          na1 = previousSumYAStratum,
          na = na * (mStratum - 1),
          nb1 = previousSumYBStratum,
          nb = nb * (mStratum - 1)
        ),
        breveThetaB = safestats:::calculateThetaBFromThetaAAndLOR(breveThetaA, log(previousOddsEstimateTotal)),
        likelihoodAlternative = likelihoodTwoProportions(
          na1 = ya,
          na = na,
          nb1 = yb,
          nb = nb,
          thetaA = breveThetaA,
          thetaB = breveThetaB
        )
      ) %>%
      select(
        mTotal,
        stratum,
        mStratum,
        na,
        nb,
        ya,
        yb,
        breveThetaA,
        breveThetaB,
        likelihoodAlternative
      )
  } else if (sharing == "nuisance") {
    stratumEstimates <- stratumData %>%
      mutate(
        breveThetaA = calculateMarginalPredProb(
          totalSuccess = previousSumYATotal,
          totalFail = previousFailATotal,
          priorSuccess = controlSuccessPrior,
          priorFail = controlFailPrior
        ),
        breveThetaB = calculateMarginalPredProb(
          totalSuccess = previousSumYBStratum,
          totalFail = nb * (mStratum - 1) - previousSumYBStratum,
          priorSuccess = interventionSuccessPrior,
          priorFail = interventionFailPrior
        ),
        likelihoodAlternative = likelihoodTwoProportions(
          na1 = ya,
          na = na,
          nb1 = yb,
          nb = nb,
          thetaA = breveThetaA,
          thetaB = breveThetaB
        )
      ) %>%
      select(
        mTotal,
        stratum,
        mStratum,
        na,
        nb,
        ya,
        yb,
        breveThetaA,
        breveThetaB,
        likelihoodAlternative
      )
  } else if (sharing == "riskDifference") {
    stratumEstimates <- stratumData %>%
      #add total risk difference estimate
      mutate(previousTotalRiskDifference = previousSumYBTotal/(nb*(mTotal - 1)) - previousSumYATotal/(na*(mTotal - 1))) %>%
      mutate(
        breveThetaA = calculateMarginalPredProb(
          totalSuccess = previousSumYAStratum,
          totalFail = na * (mStratum - 1) - previousSumYAStratum,
          priorSuccess = controlSuccessPrior,
          priorFail = controlFailPrior
        ),
        breveThetaB = breveThetaA + previousTotalRiskDifference,
        likelihoodAlternative = likelihoodTwoProportions(
          na1 = ya,
          na = na,
          nb1 = yb,
          nb = nb,
          thetaA = breveThetaA,
          thetaB = breveThetaB
        )
      ) %>%
      select(
        mTotal,
        stratum,
        mStratum,
        na,
        nb,
        ya,
        yb,
        breveThetaA,
        breveThetaB,
        likelihoodAlternative
      )
  } else if (sharing == "mixNuisance") {
    observationsAndEstimatesWithMixParameters <- stratumData %>%
      mutate(
        breveThetaAStratum = calculateMarginalPredProb(
          totalSuccess = previousSumYAStratum,
          totalFail = na * (mStratum - 1) - previousSumYAStratum,
          priorSuccess = controlSuccessPrior,
          priorFail = controlFailPrior
        ),
        breveThetaATotal = calculateMarginalPredProb(
          totalSuccess = previousSumYATotal,
          totalFail = previousFailATotal,
          priorSuccess = controlSuccessPrior,
          priorFail = controlFailPrior
        ),
        breveThetaB = calculateMarginalPredProb(
          totalSuccess = previousSumYBStratum,
          totalFail = nb * (mStratum - 1) - previousSumYBStratum,
          priorSuccess = interventionSuccessPrior,
          priorFail = interventionFailPrior
        )
      ) %>%
      arrange(mTotal)
    
    optimalAlphaAtJ <- sapply(1:nrow(observationsAndEstimatesWithMixParameters), function(rowNumber){
      findOptimalMixingRateForMinusLogLikelihood(data = observationsAndEstimatesWithMixParameters[1:rowNumber,])
    })
    
    observationsAndEstimatesWithOptimalAlpha <- cbind(observationsAndEstimatesWithMixParameters, optimalAlphaAtJ)
    
    #suppress warning about assigning NA values within case_when
    suppressWarnings(    
      stratumEstimates <- observationsAndEstimatesWithOptimalAlpha %>%
        mutate(previousOptimalAlpha = lag(optimalAlphaAtJ),
               breveThetaA = as.numeric(previousOptimalAlpha * breveThetaATotal + (1 - alpha)*breveThetaAStratum),
               likelihoodAlternative = case_when(
                 is.na(breveThetaA) | is.na(breveThetaB) ~ NA_real_,
                 TRUE ~ likelihoodTwoProportions(
                   na1 = ya,
                   na = na,
                   nb1 = yb,
                   nb = nb,
                   thetaA = breveThetaA,
                   thetaB = breveThetaB
                 )
               )
        ) %>%
        select(
          mTotal,
          stratum,
          mStratum,
          na,
          nb,
          ya,
          yb,
          breveThetaA,
          breveThetaB,
          likelihoodAlternative
        )
    )
    
  } else {
    stratumEstimates <- stratumData %>%
      mutate(
        breveThetaA = calculateMarginalPredProb(
          totalSuccess = previousSumYAStratum,
          totalFail = na * (mStratum - 1) - previousSumYAStratum,
          priorSuccess = controlSuccessPrior,
          priorFail = controlFailPrior
        ),
        breveThetaB = calculateMarginalPredProb(
          totalSuccess = previousSumYBStratum,
          totalFail = nb * (mStratum - 1) - previousSumYBStratum,
          priorSuccess = interventionSuccessPrior,
          priorFail = interventionFailPrior
        ),
        likelihoodAlternative = likelihoodTwoProportions(
          na1 = ya,
          na = na,
          nb1 = yb,
          nb = nb,
          thetaA = breveThetaA,
          thetaB = breveThetaB
        )
      ) %>%
      select(
        mTotal,
        stratum,
        mStratum,
        na,
        nb,
        ya,
        yb,
        breveThetaA,
        breveThetaB,
        likelihoodAlternative
      )
  }
  
  return(stratumEstimates)
}

#' Simulate paths of confidence bounds for the minimal difference in stratified 2x2 data.
#' @author R.J. Turner
#'
#' @param realNumberOfDataBlocks vector with number of data blocks to simulate for.
#' @param realProbs list of the probabilities and batch sizes to simulate for, for structure
#' see the example below.
#' @param alpha desired type-I error guarantee, numeric.
#' @param globalSeed seed to set for the simulations.
#' @param exponent exponent to set for the learning rate in the pseudo-Bayesian approach
#'
#' @return a tibble with data that can be used to plot a path of a lower- and upperbound
#' of a CS for the minimum effct in the stratified 2x2 setting.
#' 
#' @examples
#' #figure S4a in AISTATS 2023 paper
#' alpha <- 0.05
#' globalSeed <- 12244
#' realNumberOfDataBlocks <- c(30, 30, 30)
#' realProbs <- list(c("pa" = 0.1, "pb" = 0.5, "na" = 1, "nb" = 1),
#'                   c("pa" = 0.2, "pb" = 0.6, "na" = 1, "nb" = 1),
#'                   c("pa" = 0.3, "pb" = 0.8, "na" = 1, "nb" = 1))
#' #figure 6 in AISTATS 2023 paper
#' # realProbs <- list(c("pa" = 0.1, "pb" = 0.15, "na" = 1, "nb" = 1),
#'                   c("pa" = 0.2, "pb" = 0.6, "na" = 1, "nb" = 1),
#'                   c("pa" = 0.3, "pb" = 0.8, "na" = 1, "nb" = 1))
#'  plotData <- analyzeAndCreatePlotDataConfSetMinDifference(realNumberOfDataBlocks, realProbs, globalSeed, alpha,
#'  exponent = 2, mSwitch = 5)                  
analyzeAndCreatePlotDataConfSetMinDifference <- function(realNumberOfDataBlocks,
                                                        realProbs,
                                                        globalSeed,
                                                        alpha,
                                                        exponent = 4) {
  #obtain corresponding E-values for each stratum, per time point
  EValueResultsTheta0GreaterThan <-
    performStratified2x2CISimulation(
      realNumberOfDataBlocks,
      realProbs,
      numberForSeed = globalSeed,
      sharing = "none",
      alpha = alpha,
      returnRawEValues = TRUE,
      theta0GreaterThanDelta = TRUE
    )
  
  #also test in the other direction for the lower bound
  EValueResultsTheta0SmallerThan <-
    performStratified2x2CISimulation(
      realNumberOfDataBlocks,
      realProbs,
      numberForSeed = globalSeed,
      sharing = "none",
      alpha = alpha,
      returnRawEValues = TRUE,
      theta0SmallerThanDelta = TRUE
    )

  #weighted averaging based on previous E
  totalResult <- computeWeightedAverageEValuesOverStrata(EValueResultsTheta0GreaterThan)
  
  weightedResult <- totalResult %>%
    dplyr::filter(ETotal >= 1/alpha) %>%
    group_by(mTotal) %>%
    summarise(confidenceBound = min(delta)) %>%
    mutate(method = "posterior",
           bound = "upper")
  
  #more extreme: introduce exponential weight
  totalResult <- computeWeightedAverageEValuesOverStrata(EValueResultsTheta0GreaterThan, exponent = exponent)
  
  weightedResultExponent <- totalResult %>%
    dplyr::filter(ETotal >= 1/alpha) %>%
    group_by(mTotal) %>%
    summarise(confidenceBound = min(delta)) %>%
    mutate(method = "posteriorPolynomial",
           bound = "upper")
  
  #multiplying E Values 
  totalResult <- multiplyEValuesOverStrata(EValueResultsTheta0GreaterThan)
  
  multiplicResults <- totalResult %>%
    dplyr::filter(ETotal >= 1/alpha) %>%
    group_by(mTotal) %>%
    summarise(confidenceBound = min(delta)) %>%
    mutate(method = "multiply",
           bound = "upper")
  
  #switch 
  totalResult <- computeSwitchStrategyWithPrior(
    EValueResultsTheta0GreaterThan, exponent = 1, switchPrior = "uniform",
    switchTimesMin = 5, switchTimesMax = 30, switchTimesStep = 1
  )
  
  switchResult <- totalResult %>%
    dplyr::filter(ETotal >= 1/alpha) %>%
    group_by(mTotal) %>%
    summarise(confidenceBound = min(delta)) %>%
    mutate(method = "switch",
           bound = "upper")
  
  #take the maximum over E-values for the greater than theta0 shape to obtain an lower bound on the minimum
  totalResult <- takeMinimumEValueOverStrata(EValueResultsTheta0SmallerThan)
  
  lowerBoundResults <- totalResult %>%
    dplyr::filter(ETotal >= 1/alpha) %>%
    group_by(mTotal) %>%
    summarise(confidenceBound = max(delta)) %>%
    mutate(bound = "lower",
           method = "minimum")
  
  plotData <- rbind(data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
                      left_join(lowerBoundResults, by = "mTotal"),
                    data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
                      left_join(weightedResult, by = "mTotal"),
                    data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
                      left_join(weightedResultExponent, by = "mTotal"),
                    data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
                      left_join(multiplicResults, by = "mTotal"),
                    data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
                      left_join(switchResult, by = "mTotal")) %>%
    drop_na()
  
  return(plotData)
}

#' Given a tibble with E values per stratum, per batch number, calculate the E values
#' for a minimal of maximal difference over all strata for a switch prior strategy
#' @author R.J. Turner
#'
#' @param EValueResults output from performStratified2x2CISimulation()
#' @param exponent exponent to compute averaged E-values with before switching
#' @param switchPrior prior type, one of c("uniform", "normal")
#' @param switchTimesMin minimum switch time to place uniform prior on
#' @param switchTimesMax maximal switch time to place uniform prior on
#' @param switchTimesStep precision of switch times grid to place uniform prior on
#' @param switchPriorMean mean of normal switch prior
#' @param switchPriorSd standard deviation of normal switch prior
#'
#' @return tibble with E values per delta per time point 
#'
#' @examples
#' #figure S2a in AISTATS 2023 paper
#' alpha <- 0.05
#' globalSeed <- 12242
#' realNumberOfDataBlocks <- c(30, 30, 30)
#' realProbs <- list(
#'   c(
#'     "pa" = 0.1,
#'     "pb" = 0.15,
#'     "na" = 1,
#'     "nb" = 1
#'   ),
#'   c(
#'     "pa" = 0.2,
#'     "pb" = 0.6,
#'     "na" = 1,
#'     "nb" = 1
#'   ),
#'   c(
#'     "pa" = 0.3,
#'     "pb" = 0.8,
#'     "na" = 1,
#'     "nb" = 1
#'   )
#' )
#' EValueResultsTheta0GreaterThan <-
#'   performStratified2x2CISimulation(
#'     realNumberOfDataBlocks,
#'     realProbs,
#'     numberForSeed = globalSeed,
#'     sharing = "none",
#'     alpha = alpha,
#'     returnRawEValues = TRUE,
#'     theta0GreaterThanDelta = TRUE
#'   )
#' uniformSwitchPriorResult <-
#'   computeSwitchStrategyWithPrior(EValueResultsTheta0GreaterThan,
#'                                  exponent = 1,
#'                                  switchPrior = "uniform", 
#'                                  switchTimesMin = 1, 
#'                                  switchTimesMax = 40, 
#'                                  switchTimesStep = 1) 
computeSwitchStrategyWithPrior <- function(
  EValueResults,
  exponent = 2,
  switchPrior = "uniform", #"normal"
  switchTimesMin = 5,
  switchTimesMax = 25,
  switchTimesStep = 1,
  switchPriorMean = NA,
  switchPriorSd = NA
) {
  #prepare possible switch times and their prior
  switchTimes <- seq(switchTimesMin, switchTimesMax, by = switchTimesStep)
  if(switchPrior == "uniform") {
    dSwitchTimes <- rep(1/length(switchTimes), times = length(switchTimes))
  } else if(switchPrior == "normal"){
    dSwitchTimes <-  dnorm(switchTimes, mean = switchPriorMean, sd = switchPriorSd)
  } else if(switchPrior == "stephen") {
    dSwitchTimes <- 1/(switchTimes * (switchTimes + 1))
  }
  
  normalizedSwitchTimes <- dSwitchTimes / sum(dSwitchTimes)
  switchPriorDF <- data.frame(
    mSwitch = switchTimes,
    switchPrior = normalizedSwitchTimes
  )
  
  #prepare E observations per time point and lagged product E variable weights
  FilledObservationsPerTimePoint <- EValueResults %>%
    #unobserved data at time point m: fill in E = 1
    pivot_wider(id_cols = c(delta, mTotal), names_from = stratum, values_from = mE, names_prefix = "mE") %>%
    mutate_at(vars(starts_with("mE")), replace_na, 1) %>%
    pivot_longer(cols = c(-delta, -mTotal), names_to = "stratum", values_to = "mE") %>%
    mutate(stratum = str_remove(stratum, "mE")) %>%
    #multiply E-variables within each stratum: running product
    group_by(delta, stratum) %>%
    arrange(rev(delta),stratum, mTotal) %>%
    mutate(E = cumprod(mE)) %>%
    #create weights based on the multiplied E-variables
    group_by(delta, mTotal) %>%
    mutate(weight = E^exponent / sum(E^exponent)) %>%
    #create lagged weights
    group_by(delta, stratum) %>%
    arrange(rev(delta), mTotal) %>%
    mutate(weight = lag(weight)) %>%
    ungroup()
  
  resultsForEachMSwitch <- data.frame()
  
  for(mSwitch in switchTimes) {
    winners <- FilledObservationsPerTimePoint %>%
      filter(mTotal == mSwitch) %>%
      group_by(delta) %>%
      filter(E == max(E)) %>%
      #take the winner stratum: if several strata gave the maximum, arbitrarily pick the first
      summarise(stratumWinner = unique(stratum)[1], .groups = "drop")
    
    resultsForEachMSwitch <- rbind(
      resultsForEachMSwitch,
      FilledObservationsPerTimePoint %>%
        left_join(winners, by = "delta") %>%
        #weights: if mTotal > mSwitch, replace weight by 1 if winner, 0 if loser
        mutate(weight = case_when(
          mTotal > mSwitch & stratum == stratumWinner ~ 1,
          mTotal > mSwitch ~ 0,
          TRUE ~ weight
        )) %>%
        mutate(mEWeighted = weight * mE) %>%
        #sum weighted mE per observation moment over the strata
        group_by(delta, mTotal) %>%
        summarise(mETotal = sum(mEWeighted),
                  .groups = "drop") %>%
        #calculate the total E Variable by multiplying over the timepoints
        mutate(mETotal = replace_na(mETotal, 1)) %>%
        group_by(delta) %>%
        arrange(mTotal) %>%
        mutate(ETotal = cumprod(mETotal)) %>%
        ungroup() %>%
        select(-mETotal) %>%
        mutate(mSwitch = mSwitch)
    )
  }
  
  resultsForEachMSwitch %>%
    left_join(switchPriorDF, by = "mSwitch") %>%
    mutate(weightedESwitch = switchPrior * ETotal) %>%
    group_by(delta, mTotal) %>%
    summarise(ETotal = sum(weightedESwitch), .groups = "drop")
}

#' Simulate stratified 2x2 contingency table data and compute a confidence sequence for the mean effect.
#' @author R.J. Turner
#' 
#' @param alpha desired type-I error guarantee, numeric.
#' @param globalSeed seed to set for the simulation.
#' @param realNumberOfDataBlocks vector with number of data blocks to simulate for.
#' @param realProbs list of the probabilities and batch sizes to simulate for, for structure
#' see the example below.
#' @param deltaStarStart minimum value of mean effect to compute E-values for
#' @param deltaStarStop maximum value of mean effect to compute E-values for
#' @param deltaStarGridPrecision grid precision of mean effects to compute E-values for
#' @param minimumEMaxSteps maximum number of steps taken during universal inference
#' @param minimumEGridPrecision precision of minimization for universal inference
#'
#' @return list with tibble with the mean effect confidence interval data, and the simulated data sequence.
simulateAndComputeMeanEffectCI <- function(alpha,
                                           globalSeed,
                                           realNumberOfDataBlocks,
                                           realProbs,
                                           deltaStarStart = -0.9,
                                           deltaStarStop = 0.9,
                                           deltaStarGridPrecision = 0.01,
                                           minimumEMaxSteps = 1e3,
                                           minimumEGridPrecision = 1e-3){
  deltaStars <- seq(deltaStarStart, deltaStarStop, by = deltaStarGridPrecision)
  
  set.seed(globalSeed)
  observationsAndEstimates <-
    createStratified2x2CIData(realNumberOfDataBlocks = realNumberOfDataBlocks,
                              realProbs = realProbs, pretendAllSameStratum = FALSE) 
  
  pointEstimates <- observationsAndEstimates %>%
    createStratumEstimates(sharing = "none")
  
  #at each time point, need YAStratum, YBStratum, breveThetaAStratum, breveThetaBStratum, likelihoodAlternative
  filledForTwoStrata <- pointEstimates %>%
    select(-mStratum, -na, -nb) %>%
    pivot_wider(names_from = stratum, 
                values_from = c(ya, yb, breveThetaA, breveThetaB, likelihoodAlternative)
    ) %>%
    tidyr::fill(everything())
  
  minimumEValuesPerMPerDelta <- matrix(nrow = nrow(filledForTwoStrata), ncol = length(deltaStars)) 
  
  pb <- txtProgressBar(max = nrow(filledForTwoStrata), style = 3)
  for(m in 1:nrow(filledForTwoStrata)) {
    setTxtProgressBar(pb, m)
    dataSeenSoFar <- filledForTwoStrata[1:m,]
    
    for(deltaStarNo in 1:length(deltaStars)) {
      deltaStar <- deltaStars[deltaStarNo]
      
      #given delta star, delta1 (and delta2) domains are limited
      constraintsDelta1 <- c(min = max((deltaStar - 0.5)*2 + 1e-4, -1 + 1e-4), max = min((deltaStar + 0.5)*2 - 1e-4, 1 - 1e-4))
      currentMinDelta <- constraintsDelta1[["min"]]
      currentMaxDelta <- constraintsDelta1[["max"]]
      newCandidates <- c(currentMinDelta, currentMaxDelta)
      
      newCandidatesLogE <- sapply(newCandidates, function(currentDelta){
        getLogE2Strata(currentDelta1 = currentDelta, 
                       dataSeenSoFar = dataSeenSoFar,
                       deltaStar = deltaStar)
      })
      
      #we want to avoid searching for too long
      stepCounter <- 0
      
      while(stepCounter <= minimumEMaxSteps) {
        stepCounter <- stepCounter + 1
        
        #we add a new candidate: the average
        candidates <- c(newCandidates, mean(newCandidates))
        
        #determine the E-variable for the new delta candidate
        
        newAvgLogE <- getLogE2Strata(currentDelta1 = mean(newCandidates), 
                                     dataSeenSoFar = dataSeenSoFar,
                                     deltaStar = deltaStar)
        logEForDeltas <- c(newCandidatesLogE, newAvgLogE)
        
        #we select the two candidates with the smallest E-values
        newCandidates <- candidates[order(logEForDeltas)][1:2]
        newCandidatesLogE <- logEForDeltas[order(logEForDeltas)][1:2]
        
        #have we reached desired step precision of grid?
        if (abs(newCandidates[2] - newCandidates[1]) <= minimumEGridPrecision) {
          minimumEValuesPerMPerDelta[m,deltaStarNo] <- exp(min(newCandidatesLogE))
          break()
          
        }
      }
    }
  }
  close(pb)
  
  colnames(minimumEValuesPerMPerDelta) <- deltaStars
  minimumEValuesPerMPerDelta <- as_tibble(cbind(mTotal = filledForTwoStrata$mTotal, minimumEValuesPerMPerDelta)) %>%
    pivot_longer(cols = -mTotal, names_to = "delta", values_to = "E")
  
  meanEffectCIData <- minimumEValuesPerMPerDelta  %>%
    filter(E <= 1/alpha | is.na(E)) %>%
    group_by(mTotal) %>%
    summarise(lowerBound = min(delta), upperBound = max(delta)) %>%
    ungroup() %>%
    mutate(lowerBound = cummax(lowerBound), upperBound = cummin(upperBound))
  
  return(list(
    meanEffectCIData = meanEffectCIData,
    simulatedData = observationsAndEstimates
  ))
}

#' Simulate power of the E-variable for stratified 2x2 contingence tables for the null hypothesis that 
#' the treatment effect across all strata equals 0.
#' @author R.J. Turner
#' 
#' @param realNumberOfDataBlocks vector with number of data blocks to simulate for.
#' @param realProbs list of the probabilities and batch sizes to simulate for, for structure
#' see the example below.
#' @param sharing cross-talk type
#' @param alpha desired type-I error guarantee, numeric.
#' @param globalSeed seed to set for the simulation.
#' @param M number of simulations, integer.
#' @param combinationFunction one of c("multiply", "average", "switch")
#' @param averageWeightExponent learning rate for Bayesian averaging
#' @param mSwitch if switching at a predetermined point, the batch number to switch after
#' @param switchPriorType prior type, one of c("uniform", "normal")
#' @param switchTimesMin minimum switch time to place uniform prior on
#' @param switchTimesMax maximal switch time to place uniform prior on
#' @param switchTimesStep precision of switch times grid to place uniform prior on
#' @param switchPriorMean mean of normal switch prior
#' @param switchPriorSd standard deviation of normal switch prior
#' @param pretendAllSameStratum boolean indicating if stratification should be ignored
#'
#' @return vector with estimated power for each batch number
getPowerForClassicH0Crosstalk <- function(realNumberOfDataBlocks,
                                          realProbs,
                                          sharing,
                                          alpha,
                                          globalSeed,
                                          M,
                                          combinationFunction = c("multiply", "average", "switch"),
                                          averageWeightExponent = 1,
                                          mSwitch = NA,
                                          switchPriorType = "point", #uniform, normal
                                          switchTimesMin = NA,
                                          switchTimesMax = NA,
                                          switchTimesStep = NA,
                                          switchPriorMean = NA,
                                          switchPriorSd = NA,
                                          pretendAllSameStratum = FALSE){
  combinationFunction <- match.arg(combinationFunction)
  rejectionMatrix <- matrix(nrow = sum(realNumberOfDataBlocks), ncol = M)
  
  for(m in 1:M){
    if(sharing == "mixCrossTalkNuisance") {
      #Obtain E's with cross-talk
      CIResultCrossTalk <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                            realProbs = realProbs, 
                                                            sharing = "nuisance",alpha = alpha, numberForSeed = globalSeed + m,
                                                            returnRawEValues = TRUE,
                                                            CIGridStart = 0, CIGridStop = 0, CIGridPrecision = NA,
                                                            pretendAllSameStratum = FALSE)
      
      
      #Obtain E's without cross-talk
      CIResultStratum <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                          realProbs = realProbs, 
                                                          sharing = "none",alpha = alpha, numberForSeed = globalSeed + m,
                                                          returnRawEValues = TRUE,
                                                          CIGridStart = 0, CIGridStop = 0, CIGridPrecision = NA,
                                                          pretendAllSameStratum = FALSE)
      
      #Are going to calculate:
      #\prod_{j=1}^m ( w(I | Y^(j-1) ) *  S_{j,I} + (1- w(I | Y^(j-1) }) *  S_{j,II} )  (**)
      # with w(I| Y^(j-1) ) = p(I) E^(j-1),I/\sum over teller
      
      EPerDeltaPerMCrossTalk <- CIResultCrossTalk %>%
        multiplyEValuesOverStrata() %>%
        rename(ETotalCrossTalk = ETotal) %>%
        select(delta, mTotal, ETotalCrossTalk)
      
      EPerDeltaPerMStratum <- CIResultStratum %>%
        multiplyEValuesOverStrata() %>%
        rename(ETotalStratum = ETotal) %>%
        select(delta, mTotal, ETotalStratum)
      
      CIResultTotal <- EPerDeltaPerMCrossTalk %>%
        left_join(EPerDeltaPerMStratum, by = c("delta", "mTotal")) %>%
        #create weights
        mutate(weightCrossTalk = ETotalCrossTalk/(ETotalCrossTalk + ETotalStratum)) %>%
        #create lagged weights per delta investigated
        group_by(delta) %>%
        arrange(delta, mTotal) %>%
        mutate(laggedWeightCrossTalk = lag(weightCrossTalk)) %>%
        #calculate the mixed E-value
        mutate(ETotal = laggedWeightCrossTalk * ETotalCrossTalk + (1 - laggedWeightCrossTalk) * ETotalStratum) %>%
        mutate(ETotal = replace_na(ETotal, 1)) %>%
        select(delta, mTotal, ETotal) 
      
      #have we rejected?
      rejectionResult <- CIResultTotal %>%
        mutate(cumulativeMaxE = cummax(ETotal)) %>%
        mutate(reject = cumulativeMaxE >= 1/ alpha)
      
    } else if (sharing == "mixAll") {
      #Obtain E's with cross-talk
      CIResultControlRate <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                              realProbs = realProbs, 
                                                              sharing = "nuisance",alpha = alpha, numberForSeed = globalSeed + m,
                                                              returnRawEValues = TRUE,
                                                              CIGridStart = 0, CIGridStop = 0, CIGridPrecision = NA,
                                                              pretendAllSameStratum = FALSE)
      
      
      #Obtain E's without cross-talk
      CIResultStratum <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                          realProbs = realProbs, 
                                                          sharing = "none",alpha = alpha, numberForSeed = globalSeed + m,
                                                          returnRawEValues = TRUE,
                                                          CIGridStart = 0, CIGridStop = 0, CIGridPrecision = NA,
                                                          pretendAllSameStratum = FALSE)
      
      CIResultOdds <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                       realProbs = realProbs, 
                                                       sharing = "logOdds",alpha = alpha, numberForSeed = globalSeed + m,
                                                       returnRawEValues = TRUE,
                                                       CIGridStart = 0, CIGridStop = 0, CIGridPrecision = NA,
                                                       pretendAllSameStratum = FALSE)
      
      #Are going to calculate:
      #\prod_{j=1}^m ( w(I | Y^(j-1) ) *  S_{j,I} + (1- w(I | Y^(j-1) }) *  S_{j,II} )  (**)
      # with w(I| Y^(j-1) ) = p(I) E^(j-1),I/\sum over teller
      
      EPerDeltaPerMControlRate <- CIResultControlRate %>%
        multiplyEValuesOverStrata() %>%
        rename(ETotalControlRate = ETotal) %>%
        select(delta, mTotal, ETotalControlRate)
      
      EPerDeltaPerMStratum <- CIResultStratum %>%
        multiplyEValuesOverStrata() %>%
        rename(ETotalStratum = ETotal) %>%
        select(delta, mTotal, ETotalStratum)
      
      EPerDeltaPerMOdds <- CIResultOdds %>%
        multiplyEValuesOverStrata() %>%
        rename(ETotalOdds = ETotal) %>%
        select(delta, mTotal, ETotalOdds)
      
      CIResultTotal <- EPerDeltaPerMControlRate %>%
        left_join(EPerDeltaPerMStratum, by = c("delta", "mTotal")) %>%
        left_join(EPerDeltaPerMOdds, by = c("delta", "mTotal")) %>%
        #create weights
        mutate(weightControlRate = ETotalControlRate/(ETotalControlRate + ETotalStratum + ETotalOdds),
               weightOdds= ETotalOdds/(ETotalControlRate + ETotalStratum + ETotalOdds)) %>%
        #create lagged weights per delta investigated
        group_by(delta) %>%
        arrange(delta, mTotal) %>%
        mutate(laggedWeightControlRate = lag(weightControlRate),
               laggedWeightOdds = lag(weightOdds)) %>%
        #calculate the mixed E-value
        mutate(ETotal = laggedWeightControlRate * ETotalControlRate + 
                 (1 - laggedWeightControlRate - laggedWeightOdds) * ETotalStratum +
                 laggedWeightOdds * ETotalOdds) %>%
        mutate(ETotal = replace_na(ETotal, 1)) %>%
        select(delta, mTotal, ETotal) 
      
      #have we rejected?
      rejectionResult <- CIResultTotal %>%
        mutate(cumulativeMaxE = cummax(ETotal)) %>%
        mutate(reject = cumulativeMaxE >= 1/ alpha)
    } else {
      CIResult <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                   realProbs = realProbs, 
                                                   sharing = sharing ,alpha = alpha, numberForSeed = globalSeed + m,
                                                   returnRawEValues = TRUE,
                                                   CIGridStart = 0, CIGridStop = 0, CIGridPrecision = NA,
                                                   pretendAllSameStratum = pretendAllSameStratum)
      
      #multiply, average etc.
      if (pretendAllSameStratum) {
        rejectionResult <- CIResult %>% 
          mutate(cumulativeMaxE = cummax(E)) %>%
          mutate(reject = cumulativeMaxE >= 1/ alpha)
      }else if (combinationFunction == "multiply") {
        rejectionResult <- CIResult %>% 
          multiplyEValuesOverStrata() %>%
          mutate(cumulativeMaxE = cummax(ETotal)) %>%
          mutate(reject = cumulativeMaxE >= 1/ alpha)
      } else if (combinationFunction == "average") {
        rejectionResult <- CIResult %>% 
          computeWeightedAverageEValuesOverStrata(exponent = averageWeightExponent) %>%
          mutate(cumulativeMaxE = cummax(ETotal)) %>%
          mutate(reject = cumulativeMaxE >= 1/ alpha)
      } else if(combinationFunction == "switch" & switchPriorType == "point") {
        rejectionResult <- CIResult %>% 
          computeAverageThenSwitchEValuesOverStrata(mSwitch = mSwitch, exponent = averageWeightExponent) %>%
          mutate(cumulativeMaxE = cummax(ETotal)) %>%
          mutate(reject = cumulativeMaxE >= 1/ alpha)
      } else {
        rejectionResult <- CIResult %>% 
          computeSwitchStrategyWithPrior(exponent = averageWeightExponent, switchPrior = switchPriorType,
                                         switchTimesMin = switchTimesMin, switchTimesMax = switchTimesMax,
                                         switchTimesStep = switchTimesStep, switchPriorMean = switchPriorMean,
                                         switchPriorSd = switchPriorSd) %>%
          mutate(cumulativeMaxE = cummax(ETotal)) %>%
          mutate(reject = cumulativeMaxE >= 1/ alpha)
      }
      
    }
    
    rejectionMatrix[,m] <- rejectionResult$reject
  }
  
  rowMeans(rejectionMatrix)
}

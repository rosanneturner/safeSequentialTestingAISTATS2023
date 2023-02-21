#HELPER FUNCTIONS ------------------------------------------------------------

# a) generic helper functions ------------------------------------------------
likelihoodTwoProportions <-
  function(na1, na, nb1, nb, thetaA, thetaB) {
    exp(
      na1 * log(thetaA) + (na - na1) * log(1 - thetaA) + nb1 * log(thetaB) + (nb -
                                                                                nb1) * log(1 - thetaB)
    )
  }

loglikelihoodTwoProportions <-
  function(na1, na, nb1, nb, thetaA, thetaB) {
    na1 * log(thetaA) + (na - na1) * log(1 - thetaA) + nb1 * log(thetaB) + (nb -
                                                                              nb1) * log(1 - thetaB)
  }

derivativeKLTwoProportionsLinear <-
  function(candidateThetaA,
           delta,
           na,
           nb,
           breveThetaA,
           breveThetaB,
           c = 1) {
    candidateThetaB <- c * candidateThetaA + delta
    
    na * ((1 - breveThetaA) / (1 - candidateThetaA) - breveThetaA / candidateThetaA) + nb *
      c * ((1 - breveThetaB) / (1 - candidateThetaB) - breveThetaB / candidateThetaB)
  }

retrieveThetaARIPrRiskDifference <-
  function(delta, na, nb, breveThetaA, breveThetaB) {
    safestats:::tryOrFailWithNA(
      uniroot(
        derivativeKLTwoProportionsLinear,
        interval = c(max(c(0,-delta)) + 1e-5, min(c(1, 1 - delta)) - 1e-5),
        delta = delta,
        c = 1,
        na = na,
        nb = nb,
        breveThetaA = breveThetaA,
        breveThetaB = breveThetaB
      )$root
    )
  }
vRetrieveThetaARIPrRiskDifference <-
  Vectorize(
    retrieveThetaARIPrRiskDifference,
    vectorize.args = c("na", "nb", "breveThetaA", "breveThetaB")
  )

calculateMarginalPredProb <-
  function(totalSuccess,
           totalFail,
           priorSuccess,
           priorFail) {
    (totalSuccess + priorSuccess) / (totalSuccess + totalFail + priorSuccess + priorFail)
  }

calculateMinusLogLikelihoodForMixingRate <- function(alpha, data){
  -sum(log(likelihoodTwoProportions(
    na1 = data$ya,
    na = data$na,
    nb1 = data$yb,
    nb = data$nb,
    thetaA = (1 - alpha)*data$breveThetaAStratum + alpha * data$breveThetaATotal,
    thetaB = data$breveThetaB
  )), na.rm = TRUE)
}
findOptimalMixingRateForMinusLogLikelihood <- function(data) {
  optim(calculateMinusLogLikelihoodForMixingRate, par = 0.5, data = data, lower = 0, upper = 1, method = "L-BFGS-B")$par
}

calculateOddsRatioBasedMarginalPredProbThetaA <-
  function(gridSize,
           betaA1,
           betaA2,
           logOddsEstimate,
           na1,
           na,
           nb1,
           nb) {
    rhoGrid <- seq(1 / gridSize, 1 - 1 / gridSize, length.out = gridSize)
    rhoGridDensity <-
      dbeta(x = rhoGrid,
            shape1 = betaA1,
            shape2 = betaA2)
    priorDensity <- rhoGridDensity / sum(rhoGridDensity)
    
    thetaAgrid <- rhoGrid
    thetaBgrid <-
      sapply(thetaAgrid,
             safestats:::calculateThetaBFromThetaAAndLOR,
             lOR = logOddsEstimate)
    
    likelihoodTimesPrior <-
      exp(
        na1 * log(thetaAgrid) + (na - na1) * log(1 - thetaAgrid) +
          nb1 * log(thetaBgrid) + (nb - nb1) * log(1 -
                                                     thetaBgrid) +
          log(priorDensity)
      )
    
    #normalize
    posteriorDensity <- likelihoodTimesPrior / sum(likelihoodTimesPrior)
    
    #calculate new marginal pred. probs
    thetaA <- as.numeric(thetaAgrid %*% posteriorDensity)
    
    return(thetaA)
  }
#vectorize for quick use with dplyr
vCalculateOddsRatioBasedMarginalPredProbThetaA <-
  Vectorize(
    calculateOddsRatioBasedMarginalPredProbThetaA,
    vectorize.args = c("logOddsEstimate", "na", "na1", "nb", "nb1")
  )

retrieveTotalAndStratifiedCountsAndOddsEstimatesPerTimePoint <- function(observedData){
  observedData %>%
    group_by(stratum) %>%
    mutate(
      mStratum = row_number(),
      sumYAStratum = cumsum(ya),
      previousSumYAStratum = lag(sumYAStratum),
      sumYBStratum = cumsum(yb),
      previousSumYBStratum = lag(sumYBStratum),
      oddsEstimateStratum = sumYBStratum / (nb * mStratum - sumYBStratum) * (na * mStratum - sumYAStratum) / sumYAStratum,
      oddsEstimateStratum = case_when(
        is.na(oddsEstimateStratum) ~ NA_real_,
        is.infinite(oddsEstimateStratum) ~ NA_real_,
        TRUE ~ oddsEstimateStratum
      ),
      previousOddsEstimateStratum = lag(oddsEstimateStratum)
    ) %>%
    ungroup() %>%
    mutate(
      sumYATotal = cumsum(ya),
      previousSumYATotal = lag(sumYATotal),
      failATotal = cumsum(na - ya),
      previousFailATotal = lag(failATotal),
      sumYBTotal = cumsum(yb),
      previousSumYBTotal = lag(sumYBTotal),
      oddsEstimateTotal = sumYBTotal / (nb * mTotal - sumYBTotal) * (na * mTotal - sumYATotal) / sumYATotal,
      oddsEstimateTotal = case_when(
        is.na(oddsEstimateTotal) ~ NA_real_,
        is.infinite(oddsEstimateTotal) ~ NA_real_,
        TRUE ~ oddsEstimateTotal
      ),
      previousOddsEstimateTotal = lag(oddsEstimateTotal)
    )
}

convertCIResultsToPlotData <- function(CIResults) {
  plotData <- data.frame()
  CIResults$stratum <- as.numeric(CIResults$stratum)
  K <- max(CIResults$stratum)
  for (k in 1:K) {
    stratumResults <- CIResults %>%
      dplyr::filter(stratum == k)
    
    plotData <-
      rbind(plotData, data.frame(t(
        sapply(1:realNumberOfDataBlocks[k], function(mStratum) {
          stratumResults %>%
            dplyr::filter((mRejectStratum > mStratum) | !reject) %>%
            summarise(
              mStratum = mStratum,
              lowerBound = min(delta),
              upperBound = max(delta)
            )
        })
      ), stratum = k)) %>%
      mutate_all(as.numeric)
  }
  
  return(plotData)
}

# b) compount measures helper functions -------------------------
getWidePivottedEValuesPerStratumPerDeltaPerTimepoint <- function(rawEValueResults){
  resultDataFrame <- rawEValueResults %>%
    select(-mE, -mStratum) %>%
    pivot_wider(values_from = E, names_from = stratum, names_prefix = "E_") %>%
    group_by(delta) %>%
    arrange(mTotal) %>%
    tidyr::fill(starts_with("E")) %>%
    mutate_at(vars(starts_with("E")), replace_na, replace = 1) %>%
    ungroup()
  
  return(resultDataFrame)
}

takeMinimumEValueOverStrata <- function(rawEValueResults){
  resultDataFrame <- getWidePivottedEValuesPerStratumPerDeltaPerTimepoint(rawEValueResults)
  resultDataFrame$ETotal <- apply(resultDataFrame[,3:ncol(resultDataFrame)], MARGIN = 1, function(rowData){
    return(min(rowData))
  })
  
  return(resultDataFrame)
}


multiplyEValuesOverStrata <- function(rawEValueResults){
  resultDataFrame <- getWidePivottedEValuesPerStratumPerDeltaPerTimepoint(rawEValueResults)
  resultDataFrame$ETotal <- apply(resultDataFrame[,3:ncol(resultDataFrame)], MARGIN = 1, Vectorize(prod))
  
  return(resultDataFrame)
}

computeWeightedAverageEValuesOverStrata <- function(EValueResults, exponent = 1){
  
  resultDataFrame <-  EValueResults %>%
    #unobserved data at time point m: fill in E = 1
    pivot_wider(id_cols = c(delta, mTotal), names_from = stratum, values_from = mE, names_prefix = "mE") %>%
    mutate_at(vars(starts_with("mE")), replace_na, 1) %>%
    pivot_longer(cols = c(-delta, -mTotal), names_to = "stratum", values_to = "mE") %>%
    mutate(stratum = str_remove(stratum, "mE")) %>%
    #multiply E-variables within each stratum: running product
    group_by(delta, stratum) %>%
    arrange(stratum, rev(delta), mTotal) %>%
    mutate(E = cumprod(mE)) %>%
    #create weights based on the multiplied E-variables
    group_by(delta, mTotal) %>%
    mutate(weight = E^exponent / sum(E^exponent)) %>%
    #create lagged weights
    group_by(delta, stratum) %>%
    arrange(rev(delta), mTotal) %>%
    mutate(weight = lag(weight)) %>%
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
    ungroup()
  
  return(resultDataFrame)
}

computeAverageThenSwitchEValuesOverStrata <- function(EValueResults, mSwitch = 10, exponent = 4){
  
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
  
  winners <- FilledObservationsPerTimePoint %>%
    filter(mTotal == mSwitch) %>%
    group_by(delta) %>%
    filter(E == max(E)) %>%
    #take the winner stratum: if several strata gave the maximum, arbitrarily pick the first
    summarise(stratumWinner = unique(stratum)[1], .groups = "drop")
  
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
    ungroup()
}

getLogE2Strata <- function(currentDelta1, 
                           dataSeenSoFar,
                           deltaStar) {
  #stratum 1
  thetaARIPr1 <- vRetrieveThetaARIPrRiskDifference(delta = currentDelta1, 
                                                   breveThetaA = dataSeenSoFar$breveThetaA_1, 
                                                   breveThetaB = dataSeenSoFar$breveThetaB_1, 
                                                   na = 1, nb = 1)
  thetaBRIPr1 <- currentDelta1 + thetaARIPr1
  
  #test the E-variable on the data corresponding to the current block
  logLikelihoodRIPr1 <- loglikelihoodTwoProportions(
    na1 = dataSeenSoFar$ya_1,
    na = 1,
    nb1 = dataSeenSoFar$yb_1,
    nb = 1,
    thetaA = thetaARIPr1,
    thetaB = thetaBRIPr1
  )
  
  #stratum 2
  currentDelta2 <- 2*(deltaStar - 0.5*currentDelta1)
  thetaARIPr2 <- vRetrieveThetaARIPrRiskDifference(delta = currentDelta2, 
                                                   breveThetaA = dataSeenSoFar$breveThetaA_2, 
                                                   breveThetaB = dataSeenSoFar$breveThetaB_2, 
                                                   na = 1, nb = 1)
  thetaBRIPr2 <- currentDelta2 + thetaARIPr2
  
  #test the E-variable on the data corresponding to the current block
  logLikelihoodRIPr2 <- loglikelihoodTwoProportions(
    na1 = dataSeenSoFar$ya_2,
    na = 1,
    nb1 = dataSeenSoFar$yb_2,
    nb = 1,
    thetaA = thetaARIPr2,
    thetaB = thetaBRIPr2
  )
  
  sum(replace_na(log(dataSeenSoFar$likelihoodAlternative_1) - logLikelihoodRIPr1, 0)) +
    sum(replace_na(log(dataSeenSoFar$likelihoodAlternative_2) - logLikelihoodRIPr2, 0))
  
}
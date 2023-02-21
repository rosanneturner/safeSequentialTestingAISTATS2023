require(tidyverse)
require(magrittr)

alpha <- 0.05
globalSeed <- 12242
realNumberOfDataBlocks <- c(30, 30, 30)
realProbs <- list(
  c(
    "pa" = 0.1,
    "pb" = 0.15,
    "na" = 1,
    "nb" = 1
  ),
  c(
    "pa" = 0.2,
    "pb" = 0.6,
    "na" = 1,
    "nb" = 1
  ),
  c(
    "pa" = 0.3,
    "pb" = 0.8,
    "na" = 1,
    "nb" = 1
  )
)
trueMinDiff <-
  min(sapply(realProbs, function(item) {
    item[["pb"]] - item[["pa"]]
  }))

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

uniformSwitchPriorResult <-
  computeSwitchStrategyWithPrior(EValueResultsTheta0GreaterThan,
                                 exponent = 1,
                                 switchPrior = "uniform", 
                                 switchTimesMin = 1, 
                                 switchTimesMax = 40, 
                                 switchTimesStep = 1) %>%
  dplyr::filter(ETotal >= 1 / alpha) %>%
  group_by(mTotal) %>%
  summarise(confidenceBound = min(delta), .groups = "drop") %>%
  mutate(method = "switchUniform")

uniformSwitchPriorResultEarly <-
  computeSwitchStrategyWithPrior(EValueResultsTheta0GreaterThan,
                                 exponent = 1,
                                 switchPrior = "uniform", 
                                 switchTimesMin = 5, 
                                 switchTimesMax = 10, 
                                 switchTimesStep = 1) %>%
  dplyr::filter(ETotal >= 1 / alpha) %>%
  group_by(mTotal) %>%
  summarise(confidenceBound = min(delta), .groups = "drop") %>%
  mutate(method = "switchUniformEarly")

normalSwitchPriorResult <-
  computeSwitchStrategyWithPrior(EValueResultsTheta0GreaterThan,
                                 exponent = 1,
                                 switchPrior = "normal", 
                                 switchTimesMin = 5, 
                                 switchTimesMax = 35, 
                                 switchTimesStep = 1,
                                 switchPriorMean = 20,
                                 switchPriorSd = 5) %>%
  dplyr::filter(ETotal >= 1 / alpha) %>%
  group_by(mTotal) %>%
  summarise(confidenceBound = min(delta)) %>%
  mutate(method = "switchNormal")

normalSwitchPriorResultEarly <-
  computeSwitchStrategyWithPrior(EValueResultsTheta0GreaterThan,
                                 exponent = 1,
                                 switchPrior = "normal", 
                                 switchTimesMin = 5, 
                                 switchTimesMax = 10, 
                                 switchTimesStep = 1,
                                 switchPriorMean = 7.5,
                                 switchPriorSd = 1) %>%
  dplyr::filter(ETotal >= 1 / alpha) %>%
  group_by(mTotal) %>%
  summarise(confidenceBound = min(delta)) %>%
  mutate(method = "switchNormalEarly")

pointSwitch20Result <- computeAverageThenSwitchEValuesOverStrata(EValueResultsTheta0GreaterThan, 
                                                                 mSwitch = 20, exponent = 1) %>%
  dplyr::filter(ETotal >= 1/alpha) %>%
  group_by(mTotal) %>%
  summarise(confidenceBound = min(delta)) %>%
  mutate(method = "switchAt20")

plotData <-
  rbind(
    data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
      left_join(uniformSwitchPriorResult, by = "mTotal"),
    data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
      left_join(uniformSwitchPriorResultEarly, by = "mTotal"),
    data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
      left_join(normalSwitchPriorResult, by = "mTotal"),
    data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
      left_join(normalSwitchPriorResultEarly, by = "mTotal"),
    data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
      left_join(pointSwitch20Result, by = "mTotal")
    
  )

#S2B ----------------
nSims <- 100
alpha <- 0.05
globalSeed <- 12242
realNumberOfDataBlocks <- c(30, 30, 30)
realProbs <- list(c("pa" = 0.1, "pb" = 0.15, "na" = 1, "nb" = 1),
                  c("pa" = 0.2, "pb" = 0.6, "na" = 1, "nb" = 1),
                  c("pa" = 0.3, "pb" = 0.8, "na" = 1, "nb" = 1))
trueMinDiff <- min(sapply(realProbs, function(item){item[["pb"]] - item[["pa"]]}))
resultDataFrame <- data.frame()

pb <- txtProgressBar(max = nSims, style = 3)
for(mSim in 1:nSims) {
  EValueResultsTheta0GreaterThan <-
    performStratified2x2CISimulation(
      realNumberOfDataBlocks,
      realProbs,
      numberForSeed = globalSeed + mSim,
      sharing = "none",
      alpha = alpha,
      returnRawEValues = TRUE,
      theta0GreaterThanDelta = TRUE
    )
  
  uniformSwitchPriorResult <-
    computeSwitchStrategyWithPrior(EValueResultsTheta0GreaterThan,
                                   exponent = 1,
                                   switchPrior = "uniform", 
                                   switchTimesMin = 1, 
                                   switchTimesMax = 40, 
                                   switchTimesStep = 1) %>%
    dplyr::filter(ETotal >= 1 / alpha) %>%
    group_by(mTotal) %>%
    summarise(confidenceBound = min(delta), .groups = "drop") %>%
    mutate(method = "switchUniform")
  
  uniformSwitchPriorResultEarly <-
    computeSwitchStrategyWithPrior(EValueResultsTheta0GreaterThan,
                                   exponent = 1,
                                   switchPrior = "uniform", 
                                   switchTimesMin = 5, 
                                   switchTimesMax = 10, 
                                   switchTimesStep = 1) %>%
    dplyr::filter(ETotal >= 1 / alpha) %>%
    group_by(mTotal) %>%
    summarise(confidenceBound = min(delta), .groups = "drop") %>%
    mutate(method = "switchUniformEarly")
  
  normalSwitchPriorResult <-
    computeSwitchStrategyWithPrior(EValueResultsTheta0GreaterThan,
                                   exponent = 1,
                                   switchPrior = "normal", 
                                   switchTimesMin = 5, 
                                   switchTimesMax = 35, 
                                   switchTimesStep = 1,
                                   switchPriorMean = 20,
                                   switchPriorSd = 5) %>%
    dplyr::filter(ETotal >= 1 / alpha) %>%
    group_by(mTotal) %>%
    summarise(confidenceBound = min(delta)) %>%
    mutate(method = "switchNormal")
  
  normalSwitchPriorResultEarly <-
    computeSwitchStrategyWithPrior(EValueResultsTheta0GreaterThan,
                                   exponent = 1,
                                   switchPrior = "normal", 
                                   switchTimesMin = 5, 
                                   switchTimesMax = 10, 
                                   switchTimesStep = 1,
                                   switchPriorMean = 7.5,
                                   switchPriorSd = 1) %>%
    dplyr::filter(ETotal >= 1 / alpha) %>%
    group_by(mTotal) %>%
    summarise(confidenceBound = min(delta)) %>%
    mutate(method = "switchNormalEarly")
  
  pointSwitch20Result <- computeAverageThenSwitchEValuesOverStrata(EValueResultsTheta0GreaterThan, 
                                            mSwitch = 20, exponent = 1) %>%
    dplyr::filter(ETotal >= 1/alpha) %>%
    group_by(mTotal) %>%
    summarise(confidenceBound = min(delta)) %>%
    mutate(method = "switchAt20")
  
  plotData <-
    rbind(
      data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
        left_join(uniformSwitchPriorResult, by = "mTotal"),
      data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
        left_join(uniformSwitchPriorResultEarly, by = "mTotal"),
      data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
        left_join(normalSwitchPriorResult, by = "mTotal"),
      data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
        left_join(normalSwitchPriorResultEarly, by = "mTotal"),
      data.frame(mTotal = 1:sum(realNumberOfDataBlocks)) %>%
        left_join(pointSwitch20Result, by = "mTotal")
      
    )
  
  toStore <- plotData %>%
    mutate(distanceToMin = confidenceBound - trueMinDiff) %>%
    drop_na() %>%
    select(mTotal, method, distanceToMin) %>%
    mutate(simulation = mSim)
  
  resultDataFrame <- rbind(resultDataFrame, toStore)
  setTxtProgressBar(pb, mSim)
}
close(pb)


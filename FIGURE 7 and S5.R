library(tidyverse)

alpha <- 0.05
globalSeed <- 8375368
realNumberOfDataBlocks <- c(25, 25)
# realProbs <- list(
#   c(
#     "pa" = 0.1,
#     "pb" = 0.5,
#     "na" = 1,
#     "nb" = 1
#   ),
#   c(
#     "pa" = 0.3,
#     "pb" = 0.7,
#     "na" = 1,
#     "nb" = 1
#   )
# )
realProbs <- list(
  c(
    "pa" = 0.1,
    "pb" = 0.3,
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
realMeanDifference <- mean(sapply(realProbs, function(probSet){probSet[["pb"]] - probSet[["pa"]]}))

CIresult <- simulateAndComputeMeanEffectCI(alpha = alpha,
                                           globalSeed = globalSeed,
                                           realNumberOfDataBlocks = realNumberOfDataBlocks,
                                           realProbs = realProbs,
                                           deltaStarStart = -0.9,
                                           deltaStarStop = 0.9,
                                           deltaStarGridPrecision = 0.01,
                                           minimumEMaxSteps = 1e3,
                                           minimumEGridPrecision = 1e-3)
  
#comparison to Klingenberg method
source(file="http://sites.williams.edu/bklingen/files/2013/06/stratMHRD.r")
stratMHRDData <- CIresult[["simulatedData"]] %>%
  group_by(stratum) %>%
  summarise(y1 = sum(ya), n1 = sum(na), y2 = sum(yb), n2 = sum(nb)) %>%
  select(y2, n2, y1, n1)

#comparison to minimal difference estimation
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

totalResult <- computeSwitchStrategyWithPrior(
  EValueResultsTheta0GreaterThan, exponent = 1, switchPrior = "uniform",
  switchTimesMin = 5, switchTimesMax = 45, switchTimesStep = 1
)

switchResult <- totalResult %>%
  dplyr::filter(ETotal >= 1/alpha) %>%
  group_by(mTotal) %>%
  summarise(upperBound = min(delta))

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

totalResult <- takeMinimumEValueOverStrata(EValueResultsTheta0SmallerThan)

lowerBoundResults <- totalResult %>%
  dplyr::filter(ETotal >= 1/alpha) %>%
  group_by(mTotal) %>%
  summarise(lowerBound = max(delta))
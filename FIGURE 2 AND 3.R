require(tidyverse)
require(magrittr)

#FIGURE 2--------------------------------------
M <- 1000
alpha <- 0.05
globalSeed <- 1842781
realNumberOfDataBlocks <- c(1, 1, 1) * 40

realProbs <- list(c("pa" = 0.1, "pb" = 0.15, "na" = 1, "nb" = 1),
                  c("pa" = 0.2, "pb" = 0.6, "na" = 1, "nb" = 1),
                  c("pa" = 0.8, "pb" = 0.2, "na" = 1, "nb" = 1))

powerResult <- data.frame(
  mTotal = 1:sum(realNumberOfDataBlocks),
  multiply = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                       realProbs = realProbs,
                                       sharing = "none",
                                       alpha = alpha,
                                       globalSeed = globalSeed,
                                       M = M, combinationFunction = "multiply"),
  average = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                           realProbs = realProbs,
                                           sharing = "none",
                                           alpha = alpha,
                                           globalSeed = globalSeed,
                                           M = M, combinationFunction = "average", averageWeightExponent = 1),
  average2 = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                          realProbs = realProbs,
                                          sharing = "none",
                                          alpha = alpha,
                                          globalSeed = globalSeed,
                                          M = M, combinationFunction = "average", averageWeightExponent = 2),
  switchPoint = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                          realProbs = realProbs,
                                          sharing = "none",
                                          alpha = alpha,
                                          globalSeed = globalSeed,
                                          M = M, combinationFunction = "switch", averageWeightExponent = 1,
                                          switchPriorType = "point", mSwitch = 10),
  switchUniform = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                          realProbs = realProbs,
                                          sharing = "none",
                                          alpha = alpha,
                                          globalSeed = globalSeed,
                                          M = M, combinationFunction = "switch", averageWeightExponent = 1,
                                          switchPriorType = "uniform", switchTimesMin = 5, switchTimesMax = sum(realNumberOfDataBlocks)-5, switchTimesStep = 1),
  unstratified = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                           realProbs = realProbs,
                                           sharing = "none",
                                           alpha = alpha,
                                           globalSeed = globalSeed,
                                           M = M, combinationFunction = "multiply", pretendAllSameStratum = TRUE)
  
)

#FIGURE 3------------------------------------------
M <- 100
alpha <- 0.05
globalSeed <- 19319831
realNumberOfDataBlocks <- c(1, 1, 1) * 40

realProbs <- list(c("pa" = 0.5, "pb" = 0.01, "na" = 1, "nb" = 1),
                  c("pa" = 0.49, "pb" = 0.4, "na" = 1, "nb" = 1),
                  c("pa" = 0.51, "pb" = 0.9, "na" = 1, "nb" = 1))
powerResultCrossTalk <- data.frame(
  mTotal = 1:sum(realNumberOfDataBlocks),
  none = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                              realProbs = realProbs,
                                              sharing = "none",
                                              alpha = alpha,
                                              globalSeed = globalSeed,
                                              M = M),
  nuisance = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                              realProbs = realProbs,
                                              sharing = "nuisance",
                                              alpha = alpha,
                                              globalSeed = globalSeed,
                                              M = M),
  odds = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                              realProbs = realProbs,
                                              sharing = "logOdds",
                                              alpha = alpha,
                                              globalSeed = globalSeed,
                                              M = M),
  mix = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                              realProbs = realProbs,
                                              sharing = "mixAll",
                                              alpha = alpha,
                                              globalSeed = globalSeed,
                                              M = M),
  unstratified = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                       realProbs = realProbs,
                                       sharing = "none",
                                       alpha = alpha,
                                       globalSeed = globalSeed,
                                       M = M,
                                       pretendAllSameStratum = TRUE)
  
)

#similar odds ratio
realProbs <- list(c("pa" = 0.2, "pb" = safestats:::calculateThetaBFromThetaAAndLOR(0.2, log(4)), "na" = 1, "nb" = 1),
                  c("pa" = 0.25, "pb" = safestats:::calculateThetaBFromThetaAAndLOR(0.25, log(4.01)), "na" = 1, "nb" = 1),
                  c("pa" = 0.85, "pb" = safestats:::calculateThetaBFromThetaAAndLOR(0.85, log(2.95)), "na" = 1, "nb" = 1))
powerResultCrossTalkOddsSimilar <- data.frame(
  mTotal = 1:sum(realNumberOfDataBlocks),
  none = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                       realProbs = realProbs,
                                       sharing = "none",
                                       alpha = alpha,
                                       globalSeed = globalSeed,
                                       M = M),
  nuisance = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                           realProbs = realProbs,
                                           sharing = "nuisance",
                                           alpha = alpha,
                                           globalSeed = globalSeed,
                                           M = M),
  odds = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                       realProbs = realProbs,
                                       sharing = "logOdds",
                                       alpha = alpha,
                                       globalSeed = globalSeed,
                                       M = M),
  mix = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                      realProbs = realProbs,
                                      sharing = "mixAll",
                                      alpha = alpha,
                                      globalSeed = globalSeed,
                                      M = M),
  unstratified = getPowerForClassicH0Crosstalk(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                               realProbs = realProbs,
                                               sharing = "none",
                                               alpha = alpha,
                                               globalSeed = globalSeed,
                                               M = M,
                                               pretendAllSameStratum = TRUE)
  
)
 
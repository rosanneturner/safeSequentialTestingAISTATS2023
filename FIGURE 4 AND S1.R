require(tidyverse)
require(magrittr)

#A-C -----------------------------------------------------------------------------
alpha <- 0.05
globalSeed <- 184278474
#globalSeed <- 184278472 # odds plot
realNumberOfDataBlocks <- c(30, 30, 30)
#realNumberOfDataBlocks <- c(10, 40, 40)
realProbs <- list(c("pa" = 0.1, "pb" = 0.15, "na" = 1, "nb" = 1),
                  c("pa" = 0.2, "pb" = 0.6, "na" = 1, "nb" = 1),
                  c("pa" = 0.8, "pb" = 0.2, "na" = 1, "nb" = 1))
# realProbs <- list(c("pa" = 0.5, "pb" = 0.01, "na" = 1, "nb" = 1),
#                   c("pa" = 0.5, "pb" = 0.25, "na" = 1, "nb" = 1),
#                   c("pa" = 0.5, "pb" = 0.6, "na" = 1, "nb" = 1))
# realProbs <- list(c("pa" = 0.2, "pb" = safestats:::calculateThetaBFromThetaAAndLOR(0.2, log(2)), "na" = 1, "nb" = 1),
#                   c("pa" = 0.25, "pb" = safestats:::calculateThetaBFromThetaAAndLOR(0.25, log(2)), "na" = 1, "nb" = 1),
#                   c("pa" = 0.85, "pb" = safestats:::calculateThetaBFromThetaAAndLOR(0.85, log(2)), "na" = 1, "nb" = 1))
trueDifferences <- data.frame(
  stratum = seq_along(realProbs),
  trueRiskDiff = sapply(realProbs, function(item){item[["pb"]] - item[["pa"]]})
)

CIResultNone <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                 realProbs = realProbs, 
                                                 sharing = "none",alpha = alpha, numberForSeed = globalSeed,
                                                 CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.005) %>%
  convertCIResultsToPlotData() %>%
  mutate(sharing = "none")

CIResultShareOdds <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                      realProbs = realProbs, 
                                                      sharing = "logOdds",alpha = alpha, numberForSeed = globalSeed,
                                                      CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.005) %>%
  convertCIResultsToPlotData() %>%
  mutate(sharing = "odds")

CIResultNuisance <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                     realProbs = realProbs, 
                                                     sharing = "nuisance",alpha = alpha, numberForSeed = globalSeed,
                                                     CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.005) %>%
  convertCIResultsToPlotData() %>%
  mutate(sharing = "nuisance")

#mixing
#Obtain E's with cross-talk
CIResultControlRateForMix <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                        realProbs = realProbs, 
                                                        sharing = "nuisance",alpha = alpha, numberForSeed = globalSeed,
                                                        returnRawEValues = TRUE,
                                                        CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.005,
                                                        pretendAllSameStratum = FALSE)


#Obtain E's without cross-talk
CIResultStratumForMix <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                    realProbs = realProbs, 
                                                    sharing = "none",alpha = alpha, numberForSeed = globalSeed,
                                                    returnRawEValues = TRUE,
                                                    CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.005,
                                                    pretendAllSameStratum = FALSE)

CIResultOddsForMix <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                 realProbs = realProbs, 
                                                 sharing = "logOdds",alpha = alpha, numberForSeed = globalSeed,
                                                 returnRawEValues = TRUE,
                                                 CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.005,
                                                 pretendAllSameStratum = FALSE)

#Are going to calculate:
#\prod_{j=1}^m ( w(I | Y^(j-1) ) *  S_{j,I} + (1- w(I | Y^(j-1) }) *  S_{j,II} )  (**)
# with w(I| Y^(j-1) ) = p(I) E^(j-1),I/\sum over teller

EPerDeltaPerMControlRate <- CIResultControlRateForMix %>%
  rename(EControlRate = E) %>%
  select(delta,stratum, mStratum, EControlRate)

EPerDeltaPerMStratum <- CIResultStratumForMix %>%
  rename(EStratum = E) %>%
  select(delta,stratum, mStratum, EStratum)

EPerDeltaPerMOdds <- CIResultOddsForMix %>%
  rename(EOdds = E) %>%
  select(delta, stratum, mStratum, EOdds)

CIResultMix <- EPerDeltaPerMControlRate %>%
  left_join(EPerDeltaPerMStratum, by = c("delta", "stratum", "mStratum")) %>%
  left_join(EPerDeltaPerMOdds, by = c("delta", "stratum", "mStratum")) %>%
  #create weights
  mutate(weightControlRate = EControlRate/(EControlRate + EStratum + EOdds),
         weightOdds= EOdds/(EControlRate + EStratum + EOdds)) %>%
  #create lagged weights per delta investigated
  group_by(delta, stratum) %>%
  arrange(delta, mStratum) %>%
  mutate(laggedWeightControlRate = lag(weightControlRate),
         laggedWeightOdds = lag(weightOdds)) %>%
  #calculate the mixed E-value
  mutate(EMixed = laggedWeightControlRate * EControlRate + 
           (1 - laggedWeightControlRate - laggedWeightOdds) * EStratum +
           laggedWeightOdds * EOdds) %>%
  mutate(EMixed = replace_na(EMixed, 1)) %>%
  select(delta, stratum, mStratum, EMixed) 

CIResultMixForPlot <- CIResultMix %>%
  filter(EMixed < 1/alpha) %>%
  group_by(stratum, mStratum) %>%
  summarise(lowerBound = min(delta), upperBound = max(delta)) %>%
  group_by(stratum) %>%
  arrange(mStratum) %>%
  mutate(lowerBound = cummax(lowerBound), upperBound = cummin(upperBound)) %>%
  ungroup() %>%
  mutate(sharing = "mix")

# D-F -------------------------------------------------
alpha <- 0.05
M <- 100
globalSeed <- 184278472
realNumberOfDataBlocks <- c(30, 30, 30)
#realNumberOfDataBlocks <- c(10, 40, 40)
realProbs <- list(c("pa" = 0.1, "pb" = 0.15, "na" = 1, "nb" = 1),
                  c("pa" = 0.2, "pb" = 0.6, "na" = 1, "nb" = 1),
                  c("pa" = 0.8, "pb" = 0.2, "na" = 1, "nb" = 1))
# realProbs <- list(c("pa" = 0.5, "pb" = 0.01, "na" = 1, "nb" = 1),
#                   c("pa" = 0.5, "pb" = 0.25, "na" = 1, "nb" = 1),
#                   c("pa" = 0.5, "pb" = 0.6, "na" = 1, "nb" = 1))
# realProbs <- list(c("pa" = 0.2, "pb" = safestats:::calculateThetaBFromThetaAAndLOR(0.2, log(2)), "na" = 1, "nb" = 1),
#                   c("pa" = 0.25, "pb" = safestats:::calculateThetaBFromThetaAAndLOR(0.25, log(2)), "na" = 1, "nb" = 1),
#                   c("pa" = 0.85, "pb" = safestats:::calculateThetaBFromThetaAAndLOR(0.85, log(2)), "na" = 1, "nb" = 1))
trueDifferences <- data.frame(
  stratum = seq_along(realProbs),
  trueRiskDiff = sapply(realProbs, function(item){item[["pb"]] - item[["pa"]]})
)

coveredMatrixNone <- coveredMatrixOdds <- coveredMatrixNuisance <- matrix(nrow = sum(realNumberOfDataBlocks), ncol = M)
widthMatrixNone <- widthMatrixOdds <- widthMatrixNuisance <- matrix(nrow = sum(realNumberOfDataBlocks), ncol = M)

pb <- txtProgressBar(max = M, style = 3)
for (m in 1:M) {
  CIMonteCarloResultNone <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                         realProbs = realProbs, 
                                                         sharing = "none",alpha = alpha, numberForSeed = globalSeed + m,
                                                         CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.005) %>%
    convertCIResultsToPlotData() %>%
    left_join(trueDifferences, by = "stratum") %>%
    mutate(covered = trueRiskDiff >= lowerBound & trueRiskDiff <= upperBound,
           width = upperBound - lowerBound)
  
  coveredMatrixNone[,m] <- CIMonteCarloResultNone$covered
  widthMatrixNone[,m] <- CIMonteCarloResultNone$width
  
  CIMonteCarloResultOR <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                         realProbs = realProbs, 
                                                         sharing = "logOdds",alpha = alpha, numberForSeed = globalSeed + m,
                                                         CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.005) %>%
    convertCIResultsToPlotData() %>%
    left_join(trueDifferences, by = "stratum") %>%
    mutate(covered = trueRiskDiff >= lowerBound & trueRiskDiff <= upperBound,
           width = upperBound - lowerBound)
  
  coveredMatrixOdds[,m] <- CIMonteCarloResultOR$covered
  widthMatrixOdds[,m] <- CIMonteCarloResultOR$width
  
  CIMonteCarloResultNuis <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                         realProbs = realProbs, 
                                                         sharing = "nuisance",alpha = alpha, numberForSeed = globalSeed + m,
                                                         CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.005) %>%
    convertCIResultsToPlotData() %>%
    left_join(trueDifferences, by = "stratum") %>%
    mutate(covered = trueRiskDiff >= lowerBound & trueRiskDiff <= upperBound,
           width = upperBound - lowerBound)
  
  coveredMatrixNuisance[,m] <- CIMonteCarloResultNuis$covered
  widthMatrixNuisance[,m] <- CIMonteCarloResultNuis$width
  setTxtProgressBar(pb, m)
}
close(pb)

averageWidthPerTimePoint <- rbind(
  data.frame(
    mStratum = CIMonteCarloResultNone$mStratum,
    stratum = CIMonteCarloResultNone$stratum,
    width = rowMeans(widthMatrixNone, na.rm = TRUE),
    crosstalk = "none"
  ),
  data.frame(
    mStratum = CIMonteCarloResultNuis$mStratum,
    stratum = CIMonteCarloResultNuis$stratum,
    width = rowMeans(widthMatrixNuisance, na.rm = TRUE),
    crosstalk = "control rate"
  ),
  data.frame(
    mStratum = CIMonteCarloResultOR$mStratum,
    stratum = CIMonteCarloResultOR$stratum,
    width = rowMeans(widthMatrixOdds, na.rm = TRUE),
    crosstalk = "odds"
  )
)

widthMatrixMix <- matrix(nrow = sum(realNumberOfDataBlocks), ncol = M)
pb <- txtProgressBar(max = M, style = 3)
for (m in 1:M) {
  setTxtProgressBar(pb, m)
  #mixing
  #Obtain E's with cross-talk
  CIResultControlRateForMix <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                                realProbs = realProbs, 
                                                                sharing = "nuisance",alpha = alpha, numberForSeed = globalSeed + m,
                                                                returnRawEValues = TRUE,
                                                                CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.005,
                                                                pretendAllSameStratum = FALSE)
  
  
  #Obtain E's without cross-talk
  CIResultStratumForMix <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                            realProbs = realProbs, 
                                                            sharing = "none",alpha = alpha, numberForSeed = globalSeed + m,
                                                            returnRawEValues = TRUE,
                                                            CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.005,
                                                            pretendAllSameStratum = FALSE)
  
  CIResultOddsForMix <- performStratified2x2CISimulation(realNumberOfDataBlocks = realNumberOfDataBlocks,
                                                         realProbs = realProbs, 
                                                         sharing = "logOdds",alpha = alpha, numberForSeed = globalSeed + m,
                                                         returnRawEValues = TRUE,
                                                         CIGridStart = -0.99, CIGridStop = 0.99, CIGridPrecision = 0.005,
                                                         pretendAllSameStratum = FALSE)
  
  #Are going to calculate:
  #\prod_{j=1}^m ( w(I | Y^(j-1) ) *  S_{j,I} + (1- w(I | Y^(j-1) }) *  S_{j,II} )  (**)
  # with w(I| Y^(j-1) ) = p(I) E^(j-1),I/\sum over teller
  
  EPerDeltaPerMControlRate <- CIResultControlRateForMix %>%
    rename(EControlRate = E) %>%
    select(delta,stratum, mStratum, EControlRate)
  
  EPerDeltaPerMStratum <- CIResultStratumForMix %>%
    rename(EStratum = E) %>%
    select(delta,stratum, mStratum, EStratum)
  
  EPerDeltaPerMOdds <- CIResultOddsForMix %>%
    rename(EOdds = E) %>%
    select(delta, stratum, mStratum, EOdds)
  
  CIResultMix <- EPerDeltaPerMControlRate %>%
    left_join(EPerDeltaPerMStratum, by = c("delta", "stratum", "mStratum")) %>%
    left_join(EPerDeltaPerMOdds, by = c("delta", "stratum", "mStratum")) %>%
    #create weights
    mutate(weightControlRate = EControlRate/(EControlRate + EStratum + EOdds),
           weightOdds= EOdds/(EControlRate + EStratum + EOdds)) %>%
    #create lagged weights per delta investigated
    group_by(delta, stratum) %>%
    arrange(delta, mStratum) %>%
    mutate(laggedWeightControlRate = lag(weightControlRate),
           laggedWeightOdds = lag(weightOdds)) %>%
    #calculate the mixed E-value
    mutate(EMixed = laggedWeightControlRate * EControlRate + 
             (1 - laggedWeightControlRate - laggedWeightOdds) * EStratum +
             laggedWeightOdds * EOdds) %>%
    mutate(EMixed = replace_na(EMixed, 1)) %>%
    select(delta, stratum, mStratum, EMixed) 
  
  CIResultMixForPlot <- CIResultMix %>%
    filter(EMixed < 1/alpha) %>%
    group_by(stratum, mStratum) %>%
    summarise(lowerBound = min(delta), upperBound = max(delta), .groups = "drop") %>%
    group_by(stratum) %>%
    arrange(mStratum) %>%
    mutate(lowerBound = cummax(lowerBound), upperBound = cummin(upperBound)) %>%
    ungroup() %>%
    mutate(sharing = "mix")
  
  widthMatrixMix[,m] <- CIResultMixForPlot$upperBound - CIResultMixForPlot$lowerBound
}
close(pb)


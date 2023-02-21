require(tidyverse)
require(magrittr)

#FIGURE 6 and S4A ----------------
alpha <- 0.05
globalSeed <- 12244
realNumberOfDataBlocks <- c(30, 30, 30)
realProbs <- list(c("pa" = 0.1, "pb" = 0.5, "na" = 1, "nb" = 1),
                  c("pa" = 0.2, "pb" = 0.6, "na" = 1, "nb" = 1),
                  c("pa" = 0.3, "pb" = 0.8, "na" = 1, "nb" = 1))
# realProbs <- list(c("pa" = 0.1, "pb" = 0.15, "na" = 1, "nb" = 1),
#                   c("pa" = 0.2, "pb" = 0.6, "na" = 1, "nb" = 1),
#                   c("pa" = 0.3, "pb" = 0.8, "na" = 1, "nb" = 1))
trueMinDiff <- min(sapply(realProbs, function(item){item[["pb"]] - item[["pa"]]}))

plotData <- analyzeAndCreatePlotDataConfSetMinDifference(realNumberOfDataBlocks, realProbs, globalSeed, alpha,
                                                         exponent = 2, mSwitch = 5)

#FIGURE S3 and S4B ----------------------
nSims <- 100
alpha <- 0.05
globalSeed <- 12242
realNumberOfDataBlocks <- c(30, 30, 30)
realProbs <- list(c("pa" = 0.1, "pb" = 0.5, "na" = 1, "nb" = 1),
                  c("pa" = 0.2, "pb" = 0.6, "na" = 1, "nb" = 1),
                  c("pa" = 0.3, "pb" = 0.8, "na" = 1, "nb" = 1))
# realProbs <- list(c("pa" = 0.1, "pb" = 0.15, "na" = 1, "nb" = 1),
#                   c("pa" = 0.2, "pb" = 0.6, "na" = 1, "nb" = 1),
#                   c("pa" = 0.3, "pb" = 0.8, "na" = 1, "nb" = 1))
trueMinDiff <- min(sapply(realProbs, function(item){item[["pb"]] - item[["pa"]]}))
resultDataFrame <- data.frame()

pb <- txtProgressBar(max = nSims, style = 3)
for(mSim in 1:nSims) {
  plotData <- analyzeAndCreatePlotDataConfSetMinDifference(realNumberOfDataBlocks, realProbs, globalSeed + mSim, alpha,
                                                           exponent = 2, mSwitch = 5)
  
  toStore <- plotData %>%
    unite(col = boundType, bound, method, sep = "_") %>%
    pivot_wider(id_cols = mTotal, names_from = boundType, values_from = confidenceBound) %>%
    pivot_longer(cols = starts_with("upper"), names_to = "method", values_to = "upper") %>%
    mutate(width = upper - lower_minimum) %>%
    drop_na() %>%
    #clean up
    mutate(method = str_remove(method, "upper_")) %>%
    select(mTotal, method, width) %>%
    mutate(simulation = mSim)
  
  resultDataFrame <- rbind(resultDataFrame, toStore)
  setTxtProgressBar(pb, mSim)
}
close(pb)

averageWidthPerTimePoint <- resultDataFrame %>%
  group_by(mTotal, method) %>%
  summarise(width = mean(width, na.rm = TRUE)) %>%
  ungroup()
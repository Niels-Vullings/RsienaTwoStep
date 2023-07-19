#' @title helper functions
#'
#' @description These functions calculate the respective network statistic for ego. When multiplied with the importance of each statistic (the 'parameters') this constitutes the network evaluation of ego. See: [`ts_eval()`].
#'
#' @param cov numeric, covariate scores


DyadCensus.sna <- function(i, data, sims, period, groupName, varName, levls=1:3){
  unloadNamespace("igraph") # to avoid package clashes
  require(network)
  require(sna)
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  if (network.edgecount(x) <= 0){x <- symmetrize(x)}
  # because else triad.census(x) will lead to an error
  tc <- sna::dyad.census(x)[levls]
  # names are transferred automatically
  tc
}

TriadCensus.sna <- function(i, data, sims, period, groupName, varName, levls=1:16){
  unloadNamespace("igraph") # to avoid package clashes
  require(network)
  require(sna)
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  if (network.edgecount(x) <= 0){x <- symmetrize(x)}

  # because else triad.census(x) will lead to an error
  tc <- sna::triad.census(x)[levls]
  # names are transferred automatically
  tc
}


Nacf.sna <- function(i, data, sims, period, groupName, varName, cov){
  unloadNamespace("igraph") # to avoid package clashes
  require(network)
  require(sna)
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  if (network.edgecount(x) <= 0){x <- symmetrize(x)}
  x <- as.sociomatrix.sna(x)
  snaM1 <- sna::nacf(x, cov, type = "moran", neighborhood.type = "out", demean = TRUE)[2]
  snaM1
}




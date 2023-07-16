#' @title helper functions
#'
#' @description These functions calculate the respective network statistic for ego. When multiplied with the importance of each statistic (the 'parameters') this constitutes the network evaluation of ego. See: [`ts_eval()`].
#'
#' @param cov numeric, covariate scores
ts_centering <- function(cov) {
  attr(cov, "mean") <- mean(cov)
  centered <- cov - attributes(cov)$mean
  return(centered)
}

ts_simij <- function(cov) {
  rv <- max(cov, na.rm=TRUE) - min(cov, na.rm=TRUE)
  mat1 <- matrix(cov, nrow = length(cov), ncol = length(cov), byrow = TRUE)
  mat2 <- t(mat1)
  simij <- 1 - (abs(mat1 - mat2)/rv)
  diag(simij) <- NA
  attr(cov, "simMean") <- mean(simij, na.rm=TRUE)
  attr(cov, "range") <- rv
  return(cov)
}

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




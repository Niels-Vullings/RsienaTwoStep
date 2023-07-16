#' @title Behavioral Statistics
#'
#' @description These functions calculate the respective behavioral statistic for ego. When multiplied with the importance of each statistic (the 'parameters') this constitutes the behavioral evaluation of ego. See: [`ts_eval_beh()`].
#'
#' @details For examples on how to use these statistics see: `vignette("Introduction_RsienaTwoStep", package="RsienaTwoStep")`.
#'
#' For the mathematical definition of these network statistics see chapter 12 of the RSiena manual \insertCite{ripley2022manual}{RsienaTwoStep}.
#' @family behaviorstatistics
#' @param net matrix, the adjacency matrix representing the relations between actors. Valid values are 0 and 1.
#' @param ego numeric, the ego for which we want to calculate the network statistic.
#' @param behavior numeric, behavioral scores of actors
#'
#' @references
#' \insertRef{ripley2022manual}{RsienaTwoStep}
#' @return numeric value
#' @seealso [`ts_eval()`]
#' @examples
#' ts_degree(net=net1, ego=3)
#'
#' @importFrom Rdpack reprompt
#' @rdname ts_linear
#' @export
#' centering
ts_centering <- function(actors) {
  centered <- actors - mean(actors)
  return(centered)
}

#' @rdname ts_linear
#' @export
#calculate similarity score
fsimij <- function(actors, min, max) {
  # rv <- max(actors) - min(actors)
  rv <- max - min
  mat1 <- matrix(actors, nrow = length(actors), ncol = length(actors), byrow = TRUE)
  mat2 <- t(mat1)
  simij <- 1 - (abs(mat1 - mat2)/rv)
  return(simij)
}

#' @rdname ts_linear
#' @export
#effects
flinear <- function(ego, alters, ...) {
  actors <- c(ego, alters)  #define the network
  beh_centered <- fcentering(actors)  #center behavior scores

  statistic <- beh_centered[1]  #the actual statistic

  return(statistic)
}

#' @rdname ts_linear
#' @export
fquad <- function(ego, alters, ...) {
  actors <- c(ego, alters)  #define the network
  beh_centered <- fcentering(actors)  #center behavior scores

  statistic <- (beh_centered[1])^2  #the actual statistic

  return(statistic)
}

#' @rdname ts_linear
#' @export
favAlt <- function(ego, alters, ...) {
  actors <- c(ego, alters)
  beh_centered <- fcentering(actors)

  statistic <- beh_centered[1] * (sum(beh_centered[-1], na.rm = TRUE)/length(alters))

  return(statistic)
}

#' @rdname ts_linear
#' @export
# this is the interaction between avAlt and linear shape
favAltZ <- function(ego, alters, ...) {
  actors <- c(ego, alters)
  beh_centered <- fcentering(actors)

  statistic <- beh_centered[1]*(sum(beh_centered[-1], na.rm = TRUE)/length(alters))
  statistic <- beh_centered[1]*statistic #multiply with linear shape, which is the centered behavior score of ego

  return(statistic)
}

#' @rdname ts_linear
#' @export
favSim <- function(ego, alters, min, max) {
  actors <- c(ego, alters)  #define the network
  beh_centered <- fcentering(actors)  #center behavior scores
  simij <- fsimij(beh_centered, min, max)  #calculate the similarity scores
  diag(simij) <- NA
  msimij <- mean(simij, na.rm = TRUE)  #calculate the mean similarity score. only calculate mean on non-diagonal cells??!!
  simij_c <- simij - msimij  #center the similarity scores

  statistic <- sum(simij_c[1, ], na.rm = TRUE)/length(alters)  #the actual statistic

  return(statistic)
}

#' @rdname ts_linear
#' @export
favSimZ <- function(ego, alters, min, max) {
  actors <- c(ego, alters)  #define the network
  beh_centered <- fcentering(actors)  #center behavior scores
  simij <- fsimij(beh_centered, min, max)  #calculate the similarity scores
  diag(simij) <- NA
  msimij <- mean(simij, na.rm = TRUE)  #calculate the mean similarity score. only calculate mean on non-diagonal cells??!!
  simij_c <- simij - msimij  #center the similarity scores

  statistic <- sum(simij_c[1, ], na.rm = TRUE)/length(alters)  #the actual statistic

  statistic <- beh_centered[1]*statistic #interact with (centered) z_i


  return(statistic)
}


#' @rdname ts_linear
#' @export
favAttLower <- function(ego, alters, min, max) {
  actors <- c(ego, alters)
  beh_centered <- fcentering(actors)
  simij <- fsimij(beh_centered, min, max)
  diag(simij) <- NA

  simijL <- simij[1, ]
  simijL[beh_centered >= beh_centered[1]] <- 1
  simijL[1] <- NA
  statistic <- sum(simijL, na.rm = TRUE)/length(alters)

  return(statistic)
}

#' @rdname ts_linear
#' @export
favAttLowerZ <- function(ego, alters, min, max) {
  actors <- c(ego, alters)
  beh_centered <- fcentering(actors)
  simij <- fsimij(beh_centered, min, max)
  diag(simij) <- NA

  simijL <- simij[1, ]
  simijL[beh_centered >= beh_centered[1]] <- 1
  simijL[1] <- NA
  statistic <- sum(simijL, na.rm = TRUE)/length(alters)

  statistic <- beh_centered[1] * statistic

  return(statistic)
}

#' @rdname ts_linear
#' @export
favAttHigher <- function(ego, alters, min, max) {
  actors <- c(ego, alters)
  beh_centered <- fcentering(actors)
  simij <- fsimij(beh_centered, min, max)
  diag(simij) <- NA

  simijH <- simij[1, ]
  simijH[beh_centered <= beh_centered[1]] <- 1
  simijH[1] <- NA
  statistic <- sum(simijH, na.rm = TRUE)/length(alters)

  return(statistic)
}


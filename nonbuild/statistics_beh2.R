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
#' @param min numeric, minimum value of behavioral scores of actors. If `NULL` the empirically observed minimum is used.
#' @param max numeric, maximum value of behavioral scores of actors. If `NULL` the empirically observed maximum is used.
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
ts_centering <- function(behavior) {
  centered <- behavior - mean(behavior)
  return(centered)
}

#' @rdname ts_linear
#' @export
ts_simij <- function(behavior, min=NULL, max=NULL) {
  if (is.null(min) & is.null(max)) rv <- max(behavior) - min(behavior)
  if (!is.null(min) & !is.null(max)) rv <- max - min
  mat1 <- matrix(behavior, nrow = length(behavior), ncol = length(behavior), byrow = TRUE)
  mat2 <- t(mat1)
  simij <- 1 - (abs(mat1 - mat2)/rv)
  return(simij)
}

#' @rdname ts_linear
#' @export
ts_linear <- function(ego, alters, ...) {
  actors <- c(ego, alters)  #define the network
  beh_centered <- ts_centering(actors)  #center behavior scores

  statistic <- beh_centered[1]  #the actual statistic

  return(statistic)
}

#' @rdname ts_linear
#' @export
ts_quad <- function(ego, alters, ...) {
  actors <- c(ego, alters)  #define the network
  beh_centered <- ts_centering(actors)  #center behavior scores

  statistic <- (beh_centered[1])^2  #the actual statistic

  return(statistic)
}

#' @rdname ts_linear
#' @export
ts_avAlt <- function(ego, alters, ...) {
  actors <- c(ego, alters)
  beh_centered <- ts_centering(actors)

  statistic <- beh_centered[1] * (sum(beh_centered[-1], na.rm = TRUE)/length(alters))

  return(statistic)
}

#' @rdname ts_linear
#' @export
ts_avSim <- function(ego, alters, min, max) {
  actors <- c(ego, alters)  #define the network
  beh_centered <- ts_centering(actors)  #center behavior scores
  simij <- ts_simij(beh_centered, min, max)  #calculate the similarity scores
  diag(simij) <- NA
  msimij <- mean(simij, na.rm = TRUE)  #calculate the mean similarity score. only calculate mean on non-diagonal cells??!!
  simij_c <- simij - msimij  #center the similarity scores

  statistic <- sum(simij_c[1, ], na.rm = TRUE)/length(alters)  #the actual statistic

  return(statistic)
}


#' @rdname ts_linear
#' @export
ts_avAttLower <- function(ego, alters, min, max) {
  actors <- c(ego, alters)
  beh_centered <- ts_centering(actors)
  simij <- ts_simij(beh_centered, min, max)
  diag(simij) <- NA

  simijL <- simij[1, ]
  simijL[beh_centered >= beh_centered[1]] <- 1
  simijL[1] <- NA
  statistic <- sum(simijL, na.rm = TRUE)/length(alters)

  return(statistic)
}

#' @rdname ts_linear
#' @export
ts_avAttHigher <- function(ego, alters, min, max) {
  actors <- c(ego, alters)
  beh_centered <- ts_centering(actors)
  simij <- ts_simij(beh_centered, min, max)
  diag(simij) <- NA

  simijH <- simij[1, ]
  simijH[beh_centered <= beh_centered[1]] <- 1
  simijH[1] <- NA
  statistic <- sum(simijH, na.rm = TRUE)/length(alters)

  return(statistic)
}


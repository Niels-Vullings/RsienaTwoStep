#' @title Behavior Statistics
#'
#' @description These functions calculate the respective behavior statistic for
#'   ego. When multiplied with the importance of each statistic (the
#'   'parameters') this constitutes the network evaluation of ego. See:
#'   [`ts_eval()`].
#'
#' @details For examples on how to use these statistics see:
#'   `vignette("Introduction_RsienaTwoStep", package="RsienaTwoStep")`.
#'
#'   For the mathematical definition of these network statistics see chapter 12
#'   of the RSiena manual \insertCite{ripley2022manual}{RsienaTwoStep}.
#' @family networkstatistics
#' @param beh behavioral dependent variable
#' @param net matrix, the adjacency matrix representing the relations between
#'   actors. Valid values are 0 and 1.
#' @param cov numeric, covariate scores
#' @param ego numeric, the ego for which we want to calculate the network
#'   statistic.
#'
#' @references \insertRef{ripley2022manual}{RsienaTwoStep}
#' @return numeric value
#' @seealso [`ts_eval()`]
#' @examples
#' ts_linear(df_ccovar1$cov2, ego=3)
#'
#' @importFrom Rdpack reprompt
#' @export
ts_linear <- function(beh, ego) {
  statistic <- beh[ego]
  return(statistic)
}
attr(ts_linear, "name") <- "linear"


#' @rdname ts_linear
#' @export
ts_quad <- function(beh, ego) {
  statistic <- (beh[ego])^2
  return(statistic)
}
attr(ts_quad, "name") <- "quad"


#' @rdname ts_linear
#' @export
ts_avAlt <- function(beh, net, ego, cov=NULL) {
  alters <- which(net[ego,]==1)
  if (length(alters) == 0) {
    statistic <- 0
    } else {
  statistic <- beh[ego]* mean(beh[alters])
    }
  return(statistic)
}
attr(ts_avAlt, "name") <- "avAlt"

#' @rdname ts_linear
#' @export
ts_effFrom <- function(beh, net=NULL, ego, cov) {
  statistic <- beh[ego] * cov[ego]
  return(statistic)
}
attr(ts_effFrom, "name") <- "effFrom"


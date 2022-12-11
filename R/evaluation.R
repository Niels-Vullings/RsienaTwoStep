#' @title Evaluation function
#' @details For the mathematical definition of the evaluation function see chapter 12 of the RSiena manual \insertCite{ripley2022manual}{RsienaTwoStep}.
#' @param net matrix, the adjacency matrix representing the relations between actors. Valid values are 0 and 1.
#' @param ego numeric, the ego for which we are going to calculat how its evaluates the network.
#' @param statistics list of names of statistic functions the same length as `parameters`. See `?f_degree` for a list of available statistic functions.
#' @param parameters vector of numeric values the same length as `statistics`.
#' @references
#' \insertRef{ripley2022manual}{RsienaTwoStep}
#' @seealso [RsienaTwoStep::f_alternatives_ministep()], `f_alternatives_twostep`, `f_sims`
#' @export
#' @examples
#' f_eval(net=net1, ego=6, statistics=list(f_degree, f_recip),
#' parameters=c(-2,1)) f_eval(net=net1, ego=10, statistics=list(f_degree,
#' f_recip, f_transTrip, f_transMedTrip), parameters=c(-2,2,7,7))
f_eval <- function(net, ego, statistics, parameters) {
  # calculuate value of evaluation function
  s <- 0
  for (j in 1:length(statistics)) {
    s <- s + parameters[j] * statistics[[j]](net, ego)
  }
  return(s)
}

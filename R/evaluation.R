#' @title Evaluation function
#' @details For the mathematical definition of the evaluation function see
#'   chapter 12 of the RSiena manual
#'   \insertCite{ripley2022manual}{RsienaTwoStep}.
#' @param net matrix, the adjacency matrix representing the relations between
#'   actors. Valid values are 0 and 1.
#' @param ego numeric, the ego for which we are going to calculate how its
#'   evaluates the network.
#' @param statistics list of names of statistic functions the same length as
#'   `parameters`. See e.g. [`ts_degree()`] for a list of available statistic
#'   functions.
#' @param ccovar data frame with named time-constant covariates
#' @param parameters vector of numeric values the same length as `statistics`.
#' @references \insertRef{ripley2022manual}{RsienaTwoStep}
#' @seealso [`ts_alternatives_ministep()`], [`ts_alternatives_twostep()`],
#'   [`ts_sims()`]
#' @export
#' @examples
#' ts_eval(net=net1, ego=6, statistics=list(ts_degree, ts_recip),
#' parameters=c(-2,1))
#' ts_eval(net=net1, ego=10, ccovar=df_ccovar1, statistics=list(ts_degree, ts_recip, ts_transTrip,
#' ts_transMedTrip, list(ts_egoX, "cov1")), parameters=c(-2,2,7,7,1))

ts_eval <- function(net, ego, statistics, ccovar=NULL, parameters) {
  # calculuate value of evaluation function
  s <- 0
  for (j in 1:length(statistics)) {
    if (length(statistics[[j]])==1) s <- s + parameters[j] * statistics[[j]](net, ego)
    if (length(statistics[[j]])==2) s <- s + parameters[j] * statistics[[j]][[1]](net, ego, ccovar[,statistics[[j]][[2]]])
  }
  return(s)
}

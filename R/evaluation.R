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
#' ts_eval(net=ts_net1, ego=6, statistics=list(ts_degree, ts_recip),
#' parameters=c(-2,1))
#' ts_eval(net=ts_net1, ego=10, ccovar=df_ccovar1, statistics=list(ts_degree, ts_recip, ts_transTrip,
#' ts_transMedTrip, list(ts_egoX, "cov1")), parameters=c(-2,2,7,7,1))

ts_eval <- function(net, ego, statistics, ccovar = NULL, parameters) {

  # prepare dataset
  ccovar <- ts_prepdata(ccovar)

  # Initialize the result
  s <- 0

  # Iterate over statistics
  for (j in seq_along(statistics)) {
    stat <- statistics[[j]]

    if (length(stat) == 1) {
      # Single argument statistic function
      s <- s + parameters[j] * stat(net, ego)
    #} else if (length(stat) == 2) { #we only have stats of length 1 or 2
    } else {
      # Two-argument statistic function
      if (stat[[2]] %in% names(ccovar)) { #hence it is a ccovar
      s <- s + parameters[j] * stat[[1]](net, ego, ccovar[, stat[[2]]])
      } else { #it is a depvar
          s <- s + parameters[j] * stat[[1]](net, ego, ccovar[, 1])
        }
    }
  }

  return(s)
}


ts_eval_beh <- function(beh, net, ego, statistics, ccovar = NULL, parameters) {

  # prepare dataset
  ccovar <- ts_prepdata(ccovar)

  # Initialize the result
  s <- 0

  # Iterate over behavior statistics
  for (j in seq_along(statistics)) {
    stat <- statistics[[j]]

    if (length(stat) == 1) {
      # Single argument statistic function
      s <- s + parameters[j] * stat(beh, ego)
      #} else if (length(stat) == 2) { #we only have stats of length 1 or 2
    } else {
      # Two-argument statistic function
      cov <- NULL
      if (stat[[2]] %in% names(ccovar)) cov <- ccovar[, stat[[2]]]
      s <- s + parameters[j] * stat[[1]](beh, net, ego, cov)
    }
  }

  return(s)
}


# #faster?
# #' @export
# ts_eval2 <- function(net, ego, statistics, ccovar = NULL, parameters) {
#
#   # prepare dataset
#   ccovar <- ts_prepdata(ccovar)
#
#   netstats <- sapply(statistics, FUN = function(stat) {
#     if (length(stat) == 1) {
#       # Single argument statistic function
#       stat(net, ego)
#       # } else if (length(stat) == 2) { #we only have stats of length 1 or 2
#     } else {
#       stat[[1]](net, ego, ccovar[, stat[[2]]])
#     }
#   })
#   return(sum(parameters * netstats))
# }

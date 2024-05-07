#' @title Possible behavior after ministep, simstep or twostep
#'
#' @param beh numerical, vector representing the behavioral scores
#'   actors.
#' @param net matrix, the adjacency matrix representing the relations between
#'   actors. Valid values are 0 and 1.
#' @param ego numeric, value indicating ego (row number of net)
#' @param dist1 numeric, minimal path length between ego1 and ego2 at time1 in
#'   order to be allowed to start a coordination. If `NULL` all dyads are
#'   allowed to start a coordination (i.e. **simultaneity**).
#' @param modet1 string indicating the type of ties being evaluated at time1.
#'   "`degree`" considers all ties as undirected. "`outdegree`" only allows
#'   directed paths starting from ego1 and ending at ego2. "`indegree`" only
#'   allows directed paths starting from ego2 and ending at ego1. See:
#'   `DETAILS`.
#' @description [`ts_alternatives_ministep_beh()`] constructs the possible future
#' behavior score after a ministep of `ego`.
#' [`ts_alternatives_twostep_beh()`] constructs the possible future behavioral
#' scores after two ministeps of two egos.
#' [`ts_alternatives_simstep_beh()`] constructs all possible future behavioral
#' scores (over the complete range of the behavioral variable) of one ego.Ego is
#' thus allowed to jump from one extreme to the other
#' @return list, a list of the alternative vector representing te behavioral
#' scores
#' @export
#' @seealso [`ts_select()`], [`ts_sims()`]
#' @examples
#' ccovar <- ts_prepdata(df_ccovar1)
#' ts_alternatives_ministep_beh(beh = ccovar[, "cov2"], ego = 3)
#' @export
ts_alternatives_ministep_beh <- function(beh, ego) {
  # save the alternative behavior scores
  list_alternatives <- list(beh, beh, beh)
  # set the boundaries of the behavior
  min <- round(attributes(beh)$range2[1] - attributes(beh)$mean, 4)
  max <- round(attributes(beh)$range2[2] - attributes(beh)$mean, 4)
  # change behavior with ministep
  # I need to round everything due to floating point imprecision in duplicated
  list_alternatives[[1]][ego] <-
    ifelse(list_alternatives[[1]][ego] - 1 < min, min,
           round(list_alternatives[[1]][ego] - 1,4))
  list_alternatives[[2]][ego] <- round(list_alternatives[[2]][ego],4)
  list_alternatives[[3]][ego] <-
    ifelse(list_alternatives[[3]][ego] + 1 > max, max,
           round(list_alternatives[[3]][ego] + 1,4))
  # only keep unique alternatives
  list_alternatives <-
    list_alternatives[!duplicated(list_alternatives)]
  return(list_alternatives)
}


#' @rdname ts_alternatives_ministep_beh
#' @export
ts_alternatives_simstep_beh <- function(beh, ego) {
  list_alternatives <- list()
  beh_alt <- c(attributes(beh)$range2[1]:attributes(beh)$range2[2]) -
    attributes(beh)$mean

  for (i in 1:length(beh_alt)) {
    beh_temp <- beh
    beh_temp[ego] <- beh_alt[i]
    list_alternatives[[i]] <- beh_temp
  }
  return(list_alternatives)
}



#' @rdname ts_alternatives_ministep_beh
#' @export
ts_alternatives_twostep_beh <- function(beh, net, dist1 = NULL, modet1 = "degree") {
  # dist1: minimal distance between ego1 and ego2 at time1
  # modet1: distance at time1 based on outdegree, indegree or degree

  if (is.null(dist1)) { # complete random selection
    egos <- ts_select(net = net, steps = 2)
    results <- ts_alternatives_ministep_beh(beh = beh, ego = egos[1])
    results2 <- lapply(results, ts_alternatives_ministep_beh, ego = egos[2])
    results2 <- unlist(results2, recursive = FALSE)
  } else if (!is.null(dist1)) {
    # ego1 and ego2 should be connected at t1
    if (sum(net, na.rm = TRUE) == 0) {
      # if empty network revert to complete random selection
      egos <- ts_select(net = net, steps = 2)
      results <- ts_alternatives_ministep_beh(beh = beh, ego = egos[1])
      results2 <- lapply(results, ts_alternatives_ministep_beh, ego = egos[2])
      results2 <- unlist(results2, recursive = FALSE)
    } else {
      egos <- ts_select(net = net, steps = 2, dist1 = dist1, modet1 = modet1)
      results <- ts_alternatives_ministep_beh(beh = beh, ego = egos[1])
      results2 <- lapply(results, ts_alternatives_ministep_beh, ego = egos[2])
      results2 <- unlist(results2, recursive = FALSE)
    }
  }
  return(list(egos, results2))
}

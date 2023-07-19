#' Simulation of Network evolution
#'
#' @description
#' `ts_sims` is the workhorse function of the package `RsienaTwoStep`.
#' It simulates the network evolution `nsims` times given an existing network `net`,
#'  the defined evaluation function ([`ts_eval()`]) by `statistics` and `parameters` and
#'  the average number of possible tie changes per actor as defined by `rate`.
#'
#' @details
#' For examples on how to use `ts_sims` see: `vignette("1.Introduction_RsienaTwoStep", package="RsienaTwoStep")` or the [package website](https://jochemtolsma.github.io/RsienaTwoStep/).
#' Before you set `parallel` to TRUE make sure to set-up a cluster with the package `doParallel` (see `Examples`).
#'

#'
#' @param nsims numeric, number of simulations.
#' @param parallel TRUE/FALSE
#' @param net matrix, the adjacency matrix representing the relations between
#'   actors. Valid values are 0 and 1.
#' @param ccovar data frame with named time-constant covariates
#' @param rate numeric, the average number of possible tie-changes per actor in
#'   the simulation.
#' @param statistics, list of names of statistic functions (see e.g.
#'   [`ts_degree()`] for a list of available functions)
#' @param parameters, numeric vector the same length as `parameters`
#' @param p2step numeric vector of length 3, setting the ratio of *ministep*,
#'   *twostep* and *twoministeps*.
#' @param chain TRUE/FALSE, set to `TRUE` if you want to save all the subsequent
#'   networks (after the ministep or twostep) during the simulation. If `FALSE`
#'   only the end network is saved.
#' @param dist1 numeric, minimal path length between ego1 and ego2 at time1 in
#'   order to be allowed to start a coorporation. If `NULL` all dyads are
#'   allowed to start a cooperation.
#' @param dist2 numeric, minimal path length between ego1 and ego2 at time2 in
#'   order for twostep to be counted as coorporation.
#' @param modet1 string, indicating the type of ties being evaluated at time1.
#'   "`degree`" considers all ties as undirected. "`outdegree`" only allows
#'   directed paths starting from ego1 and ending at ego2. "`indegree`" only
#'   allows directed paths starting from ego2 and ending at ego2.
#' @param modet2 string, indicating the type of ties being evaluated at time2.
#'   "`degree`" considers all ties as undirected. "`outdegree`" only allows
#'   directed paths starting from ego1 and ending at ego2. "`indegree`" only
#'   allows directed paths starting from ego2 and ending at ego2.
#'
#' @return If `chain=FALSE` a `list` (of length `nsims`) of adjacency matrices
#' representing the final network after the simulated evolution. If `chain=TRUE`
#' a `list` of lists of adjacency matrices. Each inner list represents the
#' complete network evolution of one simulation. The outer list refers to the
#' simulation run (with length `nsims`).
#' @export
#' @seealso [`ts_alternatives_ministep()`], [`ts_alternatives_twostep()`],
#'   [`ts_alternatives_simstep()`], [`ts_eval()`]
#' @examples
#' ts_sims(net=net2, nsims=2, rate=2, parallel=FALSE, statistics=list(ts_degree, ts_recip),
#' parameters=c(-2,1), p2step=c(0,1,0))
#' @importFrom foreach %dopar%
ts_sims <- function(nsims=1000, parallel=FALSE, net, ccovar=NULL, rate, statistics, parameters, p2step=c(0,1,0), chain=FALSE, dist1=NULL, dist2=NULL, modet1="degree", modet2="degree") {

  #prepare dataset
  if (!is.null(ccovar)) {
    for (i in 1:ncol(ccovar)) {
      ccovar[,i] <- ts_centering(ccovar[,i])
      ccovar[,i] <- ts_simij(ccovar[,i])
    }
  }


  if (parallel) {
    foreach::foreach(Nsim = 1:nsims) %dopar% {
      ts_sim(net=net, ccovar=ccovar, rate=rate, statistics=statistics, parameters=parameters, p2step=p2step, chain=chain, dist1=dist1, dist2=dist2, modet1=modet1, modet2=modet2)
    }
  } else {
    sims <- list()
    for (i in 1:nsims) {
      sims[[i]] <- ts_sim(net=net, ccovar=ccovar, rate=rate, statistics=statistics, parameters=parameters, p2step=p2step, chain=chain, dist1=dist1, dist2=dist2, modet1=modet1, modet2=modet2)
    }
    return(sims)
  }
}

ts_sim <- function(net, ccovar, rate, statistics=list(ts_degree, ts_recip), parameters=c(-1,2), p2step=c(0,1,0), chain=FALSE, dist1=NULL, dist2=NULL, modet1="degree", modet2="degree"){
  nministep <- rate*nrow(net) + 1
  net_n <- net
  nets <- list()
  ministep <- 1
  iteration <- 1
  while (ministep < nministep){
    #normal ministep or 2step?
    normal <- sample(c(1,2,3), 1, prob=p2step)

    #normal model
    if (normal==1) {
      #sample agent
      ego <- ts_select(net)
      #options
      options <- ts_alternatives_ministep(net=net_n, ego=ego)
      #evaluations
      eval <- sapply(options, FUN=ts_eval,  ccovar=ccovar, ego=ego, statistics=statistics, parameters=parameters)
      eval <- eval - max(eval)
      #pick new network
      net_n <- options[[sample(1:length(eval), size=1, prob=exp(eval)/sum(exp(eval)))]]
      if (chain) {nets[[iteration]] <- net_n}
      iteration <- iteration + 1
      ministep <- ministep + 1
    }

    #model with simultaneity
    if (normal==2) {
      results <- ts_alternatives_twostep(net=net_n, dist1=dist1, dist2=dist2, modet1=modet1, modet2=modet2)
      egos <- results[[1]] #sampled dyad
      options <- results[[2]] #all possible future networks after the twostep
      if (is.null(egos) & is.null(options)) { iteration <- iteration + 1 }

      if (!is.null(egos) & !is.null(options)) { #check if it was possible to sample agents
        #evaluations ego1
        eval1 <- sapply(options, FUN=ts_eval,  ccovar=ccovar, ego=egos[1], statistics=statistics, parameters=parameters)
        #evaluations ego2
        eval2 <- sapply(options, FUN=ts_eval,  ccovar=ccovar, ego=egos[2], statistics=statistics, parameters=parameters)
        #pick new network
        eval <- eval1 + eval2
        eval <- eval - max(eval)
        net_n <- options[[sample(1:length(eval), size=1, prob=exp(eval)/sum(exp(eval)))]]
        if (chain) {nets[[iteration]] <- net_n}
        iteration <- iteration + 1
        ministep <- ministep + 2
      }
    }

    #model with two simultaneous ministeps of the same ego
    if (normal==3) {
      #sample agent
      ego <- ts_select(net)
      #options
      options <- ts_alternatives_simstep(net=net_n, ego=ego)
      #evaluations
      eval <- sapply(options, FUN=ts_eval,  ccovar=ccovar, ego=ego, statistics=statistics, parameters=parameters)
      eval <- eval - max(eval)
      #pick new network
      net_n <- options[[sample(1:length(eval), size=1, prob=exp(eval)/sum(exp(eval)))]]
      if (chain) {nets[[iteration]] <- net_n}
      iteration <- iteration + 1
      ministep <- ministep + 2
    }

  }
  if (chain) { return(nets) } else { return(net_n)}
}

#' Simulation of Network evolution
#'
#' @description
#' `f_sims` is the workhorse function of the package `RsienaTwoStep`.
#' It simulates the network evolution `nsims` times given an existing network `net`,
#'  the defined evaluation function (`f_eval`) by `statistics` and `parameters` and
#'  the average number of possible tie changes per actor as defined by `rate`.
#'
#' @details
#' For examples on how to use `f_sims` see: vignette("ABMministep-vs-ABMtwostep", package="RsienaTwoStep").
#' Before you set `parallel` to TRUE make sure to set-up a cluster with the package `doParallel` (see `Examples`).
#'

#'
#' @param nsims numeric, number of simulations.
#' @param parallel TRUE/FALSE
#' @param net matrix, the adjacency matrix representing the relations between actors. Valid values are 0 and 1.
#' @param rate numeric, the average number of possible tie-changes per actor in the simulation.
#' @param statistics, list of names of statistic functions (see e.g. `?f_degree` for a list of available functions)
#' @param parameters, vector of numeric values the same length as `parameters`
#' @param p2step numeric, value between range `[0,1]` setting the probability that a twostep will occur relative to a ministep.
#' @param chain TRUE/FALSE, set to `TRUE` if you want to save all the subsequent networks (after the ministep or twostep) during the simulation. If `FALSE` only the end network is saved.
#' @param dist1 numeric, minimal path length between ego1 and ego2 at time1 in order to be allowed to start a coorporation. If `NULL` all dyads are allowed to start a cooperation.
#' @param dist2 numeric, minimal path length between ego1 and ego2 at time2 in order for twostep to be counted as coorporation.
#' @param modet1 string, indicating the type of ties being evaluated at time1. "`degree`" considers all ties as undirected. "`outdegree`" only allows directed paths starting from ego1 and ending at ego2. "`indegree`" only allows directed paths starting from ego2 and ending at ego2.
#' @param modet2 string, indicating the type of ties being evaluated at time2. "`degree`" considers all ties as undirected. "`outdegree`" only allows directed paths starting from ego1 and ending at ego2. "`indegree`" only allows directed paths starting from ego2 and ending at ego2.
#'
#' @return
#' If `chain=FALSE` a `list` (of length `nsims`) of adjacency matrices representing the final network after the simulated evolution.
#' If `chain=TRUE` a `list` of lists of adjacency matrices. Each inner list represents the complete network evolution of one simulation. The outer list refers to the simulation run (with length `nsims`).
#' @export
#' @seealso [f_alternatives_ministep], [f_alternatives_twostep], [f_eval]
#' @examples
#' #simulation with ministep only
#' f_sims(net=net1, rate=5, statistics=list(f_degree, f_recip), parameters=c(-3,1))
#'
#' #simulation with twosteps only
#' f_sims(net=net1, rate=5, statistics=list(f_degree, f_recip), parameters=c(-3,1), p2step=1)
#'
#' #running simulations in parallel
#' n.cores <- parallel::detectCores() - 1  #save one core for other work
#' #create the cluster
#' my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
#' #register it to be used by %dopar%
#' #doParallel::registerDoParallel(cl = my.cluster)
#' f_sims(net=net1, rate=5, parallel = TRUE, statistics=list(f_degree, f_recip), parameters=c(-3,1), p2step=1)
#'
#' @importFrom foreach %dopar%
f_sims <- function(nsims=1000, parallel=FALSE, net, rate, statistics, parameters, p2step=0, chain=FALSE, dist1=NULL, dist2=NULL, modet1="degree", modet2="degree") {
  if (parallel) {
    foreach::foreach(Nsim = 1:nsims) %dopar% {
      f_sim(net=net, rate=rate, statistics=statistics, parameters=parameters, p2step=p2step, chain=chain, dist1=dist1, dist2=dist2, modet1=modet1, modet2=modet2)
    }
  } else {
    sims <- list()
    for (i in 1:nsims) {
      sims[[i]] <- f_sim(net=net, rate=rate, statistics=statistics, parameters=parameters, p2step=p2step, chain=chain, dist1=dist1, dist2=dist2, modet1=modet1, modet2=modet2)
    }
    return(sims)
  }
}

f_sim <- function(net, rate, statistics=list(f_degree, f_recip), parameters=c(-1,2), p2step=0, chain=FALSE, dist1=NULL, dist2=NULL, modet1="degree", modet2="degree"){
  nministep <- rate*nrow(net)
  net_n <- net
  nets <- list()
  ministep <- 1
  iteration <- 1
  while (ministep < nministep){
    #normal ministep or 2step?
    normal <- sample(c(1,0), 1, prob=(c(1-p2step, p2step)))

    #normal model
    if (normal) {
      #sample agent
      ego <- f_select(net)
      #options
      options <- f_alternatives_ministep(net=net_n, ego=ego)
      #evaluations
      eval <- sapply(options, FUN=f_eval, ego=ego, statistics=statistics, parameters=parameters)
      eval <- eval - max(eval)
      #pick new network
      net_n <- options[[sample(1:length(eval), size=1, prob=exp(eval)/sum(exp(eval)))]]
      if (chain) {nets[[iteration]] <- net_n}
      iteration <- iteration + 1
      ministep <- ministep + 1
    }

    #model with simultaneity
    if (!normal) {
      results <- f_alternatives_twostep(net=net_n, dist1=dist1, dist2=dist2, modet1=modet1, modet2=modet2)
      egos <- results[[1]] #sampled dyad
      options <- results[[2]] #all possible future networks after the twostep
      #evaluations ego1
      eval1 <- sapply(options, FUN=f_eval, ego=egos[1], statistics=statistics, parameters=parameters)
      #evaluations ego2
      eval2 <- sapply(options, FUN=f_eval, ego=egos[2], statistics=statistics, parameters=parameters)
      #pick new network
      eval <- eval1 + eval2
      eval <- eval - max(eval)
      net_n <- options[[sample(1:length(eval), size=1, prob=exp(eval)/sum(exp(eval)))]]
      if (chain) {nets[[iteration]] <- net_n}
      iteration <- iteration + 1
      ministep <- ministep + 2
    }

  }
  if (chain) { return(nets) } else { return(net_n)}
}

#' Possible networks after ministep of ego
#'
#' @param net matrix, the adjacency matrix representing the relations between actors. Valid values are 0 and 1.
#' @param ego numeric, value indicating ego (row number of net)
#'
#' @return list, a list of the alternative adjacency matrices after all possible ministeps of ego
#' @export
#'
#' @examples
#' f_alternatives(net=net1, ego=3)

#' @importFrom foreach %dopar%
f_sims <- function(nsims=1000, parallel=FALSE, net, rate, statistics, estimates, p2step=0, chain=FALSE, dist1=NULL, dist2=NULL, modet1="degree", modet2="degree") {
  # Actual function
  if (parallel) {
    foreach::foreach(Nsim = 1:nsims) %dopar% {
      f_sim(net=net, rate=rate, statistics=statistics, estimates=estimates, p2step=p2step, chain=chain, dist1=dist1, dist2=dist2, modet1=modet1, modet2=modet2)
    }
  } else {
    sims <- list()
    for (i in 1:nsims) {
      sims[[i]] <- f_sim(net=net, rate=rate, statistics=statistics, estimates=estimates, p2step=p2step, chain=chain, dist1=dist1, dist2=dist2, modet1=modet1, modet2=modet2)
    }
    return(sims)
  }
}

f_sim <- function(net, rate, statistics=list(f_degree, f_recip), estimates=c(-1,2), p2step=0, chain=FALSE, dist1=NULL, dist2=NULL, modet1="degree", modet2="degree"){
  #with p2step we determine the probability whether it is just a ministep (p2step=0) or whether there is simultaneity. with p2step=1 only simultaneity.

  nministep <- rate*nrow(net)
  net_n <- net

  #let us save the steps
  nets <- list()
  ministep <- 1
  iteration <- 1
  while (ministep < nministep){
    #normal ministep of 2step?
    normal <- sample(c(1,0), 1, prob=(c(1-p2step, p2step)))

    #normal model
    if (normal) {
      #sample agent
      ego <- f_select(net)
      print(ego)
      #options
      options <- f_alternatives(net=net_n, ego=ego)
      #evaluations
      eval <- sapply(options, FUN=f_eval, ego=ego, stat=statistics, params=estimates)
      #pick new network
      net_n <- options[[sample(1:length(eval), size=1, prob=exp(eval)/sum(exp(eval)))]]
      print(net_n)
      if (chain) {nets[[iteration]] <- net_n}
      iteration <- iteration + 1
      ministep <- ministep + 1
    }

    #model with simultaneity
    if (!normal) {
      #options f_alternatives_2egos <- function(net, )
      results <- f_alternatives_2egos(net=net_n, dist1=dist1, dist2=dist2, modet1=modet1, modet2=modet2)
      egos <- results[[1]]
      options <- results[[2]]
      #evaluations ego1
      eval1 <- sapply(options, FUN=f_eval, ego=egos[1], stat=statistics, params=estimates)
      #evaluations ego2
      eval2 <- sapply(options, FUN=f_eval, ego=egos[2], stat=statistics, params=estimates)
      #pick new network
      eval <- eval1 + eval2
      net_n <- options[[sample(1:length(eval), size=1, prob=exp(eval)/sum(exp(eval)))]]
      if (chain) {nets[[iteration]] <- net_n}
      iteration <- iteration + 1
      ministep <- ministep + 2
    }

  }
  if (chain) { return(nets) } else { return(net_n)}
}

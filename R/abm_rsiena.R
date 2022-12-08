#' Select two actors
#'
#' Select two actors
#'
#' @details This function selects two actors at random that are part of the network

#' @param net matrix, the adjacency matrix representing the relations between actors. Valid values are 0 and 1.
#' @param steps numeric, do we want to sample one or two actors
#'
#' @return vector of length one or two, with the actors being sampled.
#'
#' @export
#' @examples
#' f_select(net=matrix(1:100, nrow=10, ncol=10), steps=2)
#'
f_select <- function(net, steps=1){
  actors <- NA
  #sample ego1
  actors[1] <- sample(1:nrow(net), 1)
  #if we allow simultaneity
  if (steps==2) {
    actors[2] <- sample(c(1:nrow(net))[-actors[1]], 1)
  }
  return(actors)
}

f_alternatives <- function(net, ego) {

  #list of all possible future networks
  list_alternatives <- list()
  for (alter in 1:ncol(net)){
    net_alt <- net
    #change tie (break or make)
    net_alt[ego,alter] <- ifelse(net_alt[ego,alter]==0, 1, 0) #this indexing will make the loop very slow. How to speed things up??
    list_alternatives[[alter]] <- net_alt
  }
  #we dont want ego to change relation to itself but want to include a non changed network
  list_alternatives[[ego]] <- net
  return(list_alternatives)
}


f_geodist <- function(net, ego1, ego2, degree="degree"){
  if (degree=="outdegree") {
    sna::geodist(dat=net, inf.replace=Inf, count.paths=FALSE, predecessors=FALSE, ignore.eval=TRUE, na.omit=TRUE)$gdist[ego1,ego2]
  } else if (degree=="indegree") {
    sna::geodist(dat=t(net), inf.replace=Inf, count.paths=FALSE, predecessors=FALSE, ignore.eval=TRUE, na.omit=TRUE)$gdist[ego1,ego2]
  } else {
    sna::geodist(dat=net + t(net), inf.replace=Inf, count.paths=FALSE, predecessors=FALSE, ignore.eval=TRUE, na.omit=TRUE)$gdist[ego1,ego2]
  }
}

f_alternatives_2egos <- function(net, dist1=NULL, dist2=NULL, modet1="degree", modet2="degree") {
  #dist1: minimal distance between ego1 and ego2 at time1
  #modet1: distance at time1 based on outdegree, indegree or degree
  #dist2: minimal distance between ego1 and ego2 at time2
  #modet2: distance at time2 based on outdegree, indegree or degree

  if (is.null(dist1) & is.null(dist2)) { #complete random selection
    egos <- f_select(net=net, steps=2)
    results <- f_alternatives(net=net, ego=egos[1])
    results2 <- lapply(results, f_alternatives, ego=egos[2])
    results2 <- unlist(results2, recursive = FALSE)
  }

  if (!is.null(dist1) & is.null(dist2)) { #ego1 and ego2 should be connected at t1
    succes <- FALSE
    repeat {
      egos <- f_select(net=net, steps=2)
      dist_t1 <- f_geodist(net=net, ego1=egos[1], ego2=egos[2], degree=modet1)
      if (dist_t1<=dist1) { #check if connected at t1, if not sample again, if true construct alternative nets
        succes <- TRUE
        results <- f_alternatives(net=net, ego=egos[1])
        results2 <- lapply(results, f_alternatives, ego=egos[2])
        results2 <- unlist(results2, recursive = FALSE)
        if (succes) break
      }
    }
  }

  if (!is.null(dist1) & !is.null(dist2)) { #ego1 and ego2 should be connected at t1 OR t2
    #if already connected at t1, easy:
    repeat {
      egos <- f_select(net=net, steps=2)
      dist_t1 <- f_geodist(net=net, ego1=egos[1], ego2=egos[2], degree=modet1)
      #we will contstruct all possible nets after the twosteps, even if we do not know beforehand if we need them.
      results <- f_alternatives(net=net, ego=egos[1])
      results2 <- lapply(results, f_alternatives, ego=egos[2])

      if (dist_t1<=dist1) { #if connected at t1, we can simply use all constructed alternative nets.
        results2 <- unlist(results2, recursive = FALSE)
        break
        } else {
          #only keep those alternative networks were distance between ego1 and 2 smaller than dist2
          keep <- rapply(results2, f=f_geodist, ego1=egos[1], ego2=egos[2], degree=modet2, how=c("list"))
          results2 <- unlist(results2, recursive=FALSE)
          results2 <- results2[unlist(keep, recursive=FALSE)<=dist2]
        }
        if (length(results2)>=1) break #check if we have at least one alternative network
      }
    }
  return(list(egos,results2))
}


f_degree <- function(net, ego) {
  statistic <- sum(net[ego,])
  return(statistic)
}

f_recip <- function(net, ego) {
  statistic <- sum(net[ego,]==1 & t(net)[ego,]==1)
  return(statistic)
}


f_outAct <- function(net, ego) {
  statistic <- sum(net[ego,])^2
  return(statistic)
}


f_inAct <- function(net, ego) {
  statistic <- sum(t(net)[ego,])*sum(net[ego,])
  return(statistic)
}

f_outPop <- function(net, ego) {
  outdegrees <- rowSums(net) #outdegrees of alters
  statistic <- sum(net[ego,] * outdegrees)
  return(statistic)
}

f_inPop <- function(net, ego) {
  indegrees <- colSums(net) #indegrees of alters
  statistic <- sum(net[ego,] * indegrees)
  return(statistic)
}

f_transTrip <- function(net, ego) {
  statistic <- 0
  alters <- which(net[ego,]==1)
  if (length(alters)>1) {
    #check if alters are connected
    for (alter1 in alters) {
      for (alter2 in alters) {
        statistic <- statistic + net[alter1, alter2]
      }
    }
  }
  return(statistic)
}

f_transMedTrip <- function(net, ego) {
  statistic <- 0
  alters1 <- which(net[ego,]==1) #ego connected to by outdegree
  alters2 <- which(t(net)[ego,]==1) #ego connected to by indegree
  if (length(alters1)>0 & length(alters2)>0) {
    #check if alters are connected
    for (alter1 in alters1) {
      for (alter2 in alters2) {
        statistic <- statistic + net[alter1, alter2]
      }
    }
  }
  return(statistic)
}

f_eval <- function(net, ego, stat, params) {
  # calculuate value of evaluation function
  s <- 0
  for (j in 1:length(stat)) {
    s <- s + params[j] * stat[[j]](net, ego)
  }
  return(s)
}


DyadCensus.sna <- function(i, data, sims, period, groupName, varName, levls=1:3){
  unloadNamespace("igraph") # to avoid package clashes
  x <- network::network.extraction(i, data, sims, period, groupName, varName)
  if (network::network.edgecount(x) <= 0){x <- sna::symmetrize(x)}
  # because else triad.census(x) will lead to an error
  tc <- sna::dyad.census(x)[levls]
  # names are transferred automatically
  tc
}

TriadCensus.sna <- function(i, data, sims, period, groupName, varName, levls=1:16){
  unloadNamespace("igraph") # to avoid package clashes
  x <- network::network.extraction(i, data, sims, period, groupName, varName)
  if (network::network.edgecount(x) <= 0){x <- sna::symmetrize(x)}
  # because else triad.census(x) will lead to an error
  tc <- sna::triad.census(x)[levls]
  # names are transferred automatically
  tc
}

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

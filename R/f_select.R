#' Select actor(s)
#'
#' Select one or two actors that will perform a (simultaneous) ministep.
#'
#' @details This function selects one actor (`steps=1`) or two actors (`steps=2`) at random that are part of the network. If it is not possible to select two actors to start a coordination (not connected dyads at time1), the function will tell you so.

#' @param net matrix, the adjacency matrix representing the relations between actors. Valid values are 0 and 1.
#' @param steps numeric, do we want to sample one or two actors
#' @param dist1 numeric, minimal path length between ego1 and ego2 at time1 in order to be allowed to start a coordination. If `NULL` just simultaneity. If not `NULL` `steps` is set to 2.
#' @param modet1 string indicating the type of ties being evaluated at time1. "`degree`" considers all ties as undirected. "`outdegree`" only allows directed paths starting from ego1 and ending at ego2. "`indegree`" only allows directed paths starting from ego2 and ending at ego2.

#'
#' @return vector of length one or two, with the actors being sampled.
#'
#' @export
#' @examples
#' ts_select(net=net1, steps=1)
#' ts_select(net=net1, steps=2)
#' ts_select(net=net1, dist1=2, modet1="degree")
#'
ts_select <- function(net, steps=1, dist1=NULL, modet1="degree"){
  actors <- NA
  if (is.null(dist1)) {
    #sample ego1
    actors[1] <- sample(1:nrow(net), 1)
    #if we allow simultaneity
    if (steps==2) {
      actors[2] <- sample(c(1:nrow(net))[-actors[1]], 1)
    }
  }
  if (!is.null(dist1)) {
    #determine distance matrix
    if (modet1=="degree") {
      net <- net + t(net)
    } else if (modet1=="indegree") {
      net <- t(net)
    }
    distmat <- sna::geodist(net, inf.replace=dist1 + 1)$gdist
    diag(distmat) <- dist1 + 1
    distmat <- distmat <= dist1
    #sample egos
    tryCatch({
      actors[1] <- sample(1:nrow(net), 1, prob = rowSums(distmat)>1) #only egos who have at least one alter fitting the criteria
      actors[2] <- sample(1:nrow(net), 1, prob = distmat[actors[1],]) #only alters that fit the criteria
      },
      error = function(e) {
               cat("We have problems sampling agents to start coordination. \n Probably your parameters for the ABM are ill-chosen and the density of the network became too low.")
      })
    }
  return(actors)
}


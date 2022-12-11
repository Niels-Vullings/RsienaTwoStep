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
#' f_select(net=net1, steps=2)
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


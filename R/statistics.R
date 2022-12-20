#' @title Network Statistics
#'
#' @description These functions calculate the respective network statistic for ego. When multiplied with the importance of each statistic (the 'parameters') this constitutues the network evaluation of ego. See: [`f_eval()`].
#'
#' @details For examples on how to use these statistics see: [vignette("Introduction_RsienaTwoStep", package="RsienaTwoStep")]:.
#'
#' For the mathematical definition of these network statistics see chapter 12 of the RSiena manual \insertCite{ripley2022manual}{RsienaTwoStep}.
#' @family networkstatistics
#' @param net matrix, the adjacency matrix representing the relations between actors. Valid values are 0 and 1.
#' @param ego numeric, the ego for which we want to calculate the network statistic.
#'
#' @references
#' \insertRef{ripley2022manual}{RsienaTwoStep}
#' @return numeric value
#' @seealso [`f_eval()`]
#' @examples
#' f_degree(net=net1, ego=3)
#'
#' @importFrom Rdpack reprompt
#' @export
f_degree <- function(net, ego) {
  statistic <- sum(net[ego,])
  return(statistic)
}

#' @rdname f_degree
#' @export
f_recip <- function(net, ego) {
  statistic <- sum(net[ego,]==1 & t(net)[ego,]==1)
  return(statistic)
}

#' @rdname f_degree
#' @export
f_outAct <- function(net, ego) {
  statistic <- sum(net[ego,])^2
  return(statistic)
}

#' @rdname f_degree
#' @export
f_inAct <- function(net, ego) {
  statistic <- sum(t(net)[ego,])*sum(net[ego,])
  return(statistic)
}

#' @rdname f_degree
#' @export
f_outPop <- function(net, ego) {
  outdegrees <- rowSums(net) #outdegrees of alters
  statistic <- sum(net[ego,] * outdegrees)
  return(statistic)
}

#' @rdname f_degree
#' @export
f_inPop <- function(net, ego) {
  indegrees <- colSums(net) #indegrees of alters
  statistic <- sum(net[ego,] * indegrees)
  return(statistic)
}

#' @rdname f_degree
#' @export
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

#' @rdname f_degree
#' @export
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

#' @title Network Statistics
#'
#' @description These functions calculate the respective network statistic for ego. When multiplied with the importance of each statistic (the 'parameters') this constitutues the network evaluation of ego. See: [`ts_eval()`].
#'
#' @details For examples on how to use these statistics see: `vignette("1.Introduction_RsienaTwoStep", package="RsienaTwoStep")`.
#'
#' For the mathematical definition of these network statistics see chapter 12 of the RSiena manual \insertCite{ripley2022manual}{RsienaTwoStep}.
#' @family networkstatistics
#' @param net matrix, the adjacency matrix representing the relations between actors. Valid values are 0 and 1.
#' @param ego numeric, the ego for which we want to calculate the network statistic.
#'
#' @references
#' \insertRef{ripley2022manual}{RsienaTwoStep}
#' @return numeric value
#' @seealso [`ts_eval()`]
#' @examples
#' ts_degree(net=net1, ego=3)
#'
#' @importFrom Rdpack reprompt
#' @export
ts_degree <- function(net, ego) {
  statistic <- sum(net[ego,])
  return(statistic)
}

#' @rdname ts_degree
#' @export
ts_recip <- function(net, ego) {
  statistic <- sum(net[ego,]==1 & t(net)[ego,]==1)
  return(statistic)
}

#' @rdname ts_degree
#' @export
ts_outAct <- function(net, ego) {
  statistic <- sum(net[ego,])^2
  return(statistic)
}

#' @rdname ts_degree
#' @export
ts_inAct <- function(net, ego) {
  statistic <- sum(t(net)[ego,])*sum(net[ego,])
  return(statistic)
}

#' @rdname ts_degree
#' @export
ts_outPop <- function(net, ego) {
  outdegrees <- rowSums(net) #outdegrees of alters
  statistic <- sum(net[ego,] * outdegrees)
  return(statistic)
}

#' @rdname ts_degree
#' @export
ts_inPop <- function(net, ego) {
  indegrees <- colSums(net) #indegrees of alters
  statistic <- sum(net[ego,] * indegrees)
  return(statistic)
}

#' @rdname ts_degree
#' @export
ts_transTrip <- function(net, ego) {
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

#' @rdname ts_degree
#' @export
ts_transMedTrip <- function(net, ego) {
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

#' @rdname ts_degree
#' @export
ts_transRecTrip <- function(net, ego) {
  # i<->j, i->h, h->j
  statistic <- 0
  alters <- which(net[ego,]==1)
  if (length(alters)>1) {
    #check if alters are connected
    for (alter1 in alters) {
      for (alter2 in alters) {
        #check first if alter is connected to ego (check for reciprocal tie i <-> j.
        if (net[alter1, ego]==1) {
          statistic <- statistic + net[alter2, alter1]
        }
      }
    }
  }
  return(statistic)
}

#' @rdname ts_degree
#' @export
ts_cycle3 <- function(net, ego) {
  # i->j, j->h, h->i
  statistic <- 0
  altersi <- which(net[ego,]==1) #identify alters of ego
  if (length(altersi)>0) {
    for (alter1 in altersi) {
      net_temp <- net
      net_temp[alter1, ego] <- 0
      altersj <- which(net_temp[alter1,]==1) #identify alters of alter but not including ego
      if (length(altersj)>0) {
        for (alter2 in altersj) {
          statistic <- statistic + net[alter2, ego]
        }
      }
    }
  }
  return(statistic)
}



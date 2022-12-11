#' @title Possible networks after ministep of ego
#'
#' @param net matrix, the adjacency matrix representing the relations between actors. Valid values are 0 and 1.
#' @param ego numeric, value indicating ego (row number of net)
#' @param dist1 numeric, minimal path length between ego1 and ego2 at time1 in order to be allowed to start a coorporation. If `NULL` all dyads are allowed to start a cooperation.
#' @param dist2 numeric, minimal path length between ego1 and ego2 at time2 in order for twostep to be counted as coorporation. See `DETAILS`.
#' @param modet1 string indicating the type of ties being evaluated at time1. "`degree`" considers all ties as undirected. "`outdegree`" only allows directed paths starting from ego1 and ending at ego2. "`indegree`" only allows directed paths starting from ego2 and ending at ego2. See: `DETAILS`.
#' @param modet2 string, indicating the type of ties being evaluated at time2. "`degree`" considers all ties as undirected. "`outdegree`" only allows directed paths starting from ego1 and ending at ego2. "`indegree`" only allows directed paths starting from ego2 and ending at ego2. See: `DETAILS`.
#'
#' @description
#' `f_alternatives_ministep` constructs the possible future networks at time2 after a ministep of `ego` given the network `net` at time1.
#' `f_alternatives_twostep` constructs the possible future networks at time2 after a twostep of two internally sampled egos (via [`f_select`]) given the network `net` at time1.
#' @details
#' `f_alternatives_ministep` mimics the ministep assumption as implemented in the SAOM of [`RSiena::siena07()`] \insertCite{ripley2022manual}{RsienaTwoStep}.
#' `f_alternatives_twostep` allows two actors to simultaneously make a ministep, that is a **twostep**.
#' Further restrictions can be set to which actors are allowed to make a twostep:
#' 1. Two random actors
#' 2. Two actors who are connected at time1 in a specific way (determined by `dist1` and `mode1`)
#'
#' A special case is to allow two random actors to make a twostep but only consider some of the alternative future networks as the result of cooperation, namely those networks that:
#' 1. result from actors who are connected at time1 in a specific way (determined by `dist1` and `mode1`), **or**
#' 2. where actors are connected at time2 in a specific way (determined by `dist2` and `mode2`)
#' @return list, a list of the alternative adjacency matrices after all possible ministeps of ego (`f_alternatives_ministep`) or after all possible twosteps of two egos (`f_alternatives_twostep`)
#' @export
#' @references
#' \insertRef{ripley2022manual}{RsienaTwoStep}
#'
#' @seealso f_select
#' @examples
#' f_alternatives_ministep(net=net1, ego=3)
f_alternatives_ministep <- function(net, ego) {

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

#' @rdname f_alternatives_ministep
#' @export
f_alternatives_twostep <- function(net, dist1=NULL, dist2=NULL, modet1="degree", modet2="degree") {
  #dist1: minimal distance between ego1 and ego2 at time1
  #modet1: distance at time1 based on outdegree, indegree or degree
  #dist2: minimal distance between ego1 and ego2 at time2
  #modet2: distance at time2 based on outdegree, indegree or degree

  if (is.null(dist1) & is.null(dist2)) { #complete random selection
    egos <- f_select(net=net, steps=2)
    results <- f_alternatives_ministep(net=net, ego=egos[1])
    results2 <- lapply(results, f_alternatives_ministep, ego=egos[2])
    results2 <- unlist(results2, recursive = FALSE)
  }

  if (!is.null(dist1) & is.null(dist2)) { #ego1 and ego2 should be connected at t1
    succes <- FALSE
    repeat {
      egos <- f_select(net=net, steps=2)
      dist_t1 <- f_geodist(net=net, ego1=egos[1], ego2=egos[2], degree=modet1)
      if (dist_t1<=dist1) { #check if connected at t1, if not sample again, if true construct alternative nets
        succes <- TRUE
        results <- f_alternatives_ministep(net=net, ego=egos[1])
        results2 <- lapply(results, f_alternatives_ministep, ego=egos[2])
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
      results <- f_alternatives_ministep(net=net, ego=egos[1])
      results2 <- lapply(results, f_alternatives_ministep, ego=egos[2])

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

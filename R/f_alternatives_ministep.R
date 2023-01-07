#' @title Possible networks after ministep, simstep or twostep
#'
#' @param net matrix, the adjacency matrix representing the relations between
#'   actors. Valid values are 0 and 1.
#' @param ego numeric, value indicating ego (row number of net)
#' @param dist1 numeric, minimal path length between ego1 and ego2 at time1 in
#'   order to be allowed to start a coordination. If `NULL` all dyads are
#'   allowed to start a coordination (i.e. **simultaneity**).
#' @param dist2 numeric, minimal path length between ego1 and ego2 at time2 in
#'   order for twostep to be counted as coordination. See `DETAILS`.
#' @param modet1 string indicating the type of ties being evaluated at time1.
#'   "`degree`" considers all ties as undirected. "`outdegree`" only allows
#'   directed paths starting from ego1 and ending at ego2. "`indegree`" only
#'   allows directed paths starting from ego2 and ending at ego1. See:
#'   `DETAILS`.
#' @param modet2 string, indicating the type of ties being evaluated at time2.
#'   "`degree`" considers all ties as undirected. "`outdegree`" only allows
#'   directed paths starting from ego1 and ending at ego2. "`indegree`" only
#'   allows directed paths starting from ego2 and ending at ego1. See:
#'   `DETAILS`.
#'
#' @description [`ts_alternatives_ministep()`] constructs the possible future
#' networks at time2 after a ministep of `ego` given the network `net` at time1.
#' [`ts_alternatives_twostep()`] constructs the possible future networks at time2
#' after a twostep of two internally sampled egos (via [`ts_select()`]) given the
#' network `net` at time1.
#' [`ts_alternatives_simstep()`] constructs the possible future
#' networks at time2 after two simultaneous ministeps of the same `ego` given the network `net` at time1.
#' @details [`ts_alternatives_ministep()`] mimics the ministep assumption as
#' implemented in the SAOM of [`RSiena::siena07()`]
#' \insertCite{ripley2022manual}{RsienaTwoStep}.
#' [`ts_alternatives_twostep()`] allows
#' two actors to simultaneously make a ministep, that is a **twostep**.
#' The function implements three types of coordination:
#'  1. **simultaneity**: when two actors are picked at random to simultaneously make
#' a ministep;
#' 2. ***weak* coordination**: two actors are picked at random to
#' simultaneously make a ministep but only specific possible future networks are
#' regarded as the result of coordination (as determined by `dist1`, `dist2`
#' `modet1` and `modet2`) and included in the choice set of the two actors;
#' 3. ***strict* coordination**: only actors are sampled to make a twostep who are
#' connected at time1 (as determined by `dist1` and `modet1`).
#'
#' [`ts_alternatives_simstep()`] allows one actor to make two subsequent ministeps and thus
#' opens the door for strategic actions. That is, the first ministep may not contribute to
#' increased satisfaction of the actor (the network after the first ministep is not evaluated
#' more favorably than the original network) but the subsequent ministep may make up for this.
#'
#' @return list, a list of the alternative adjacency matrices after all possible ministeps of ego (`ts_alternatives_ministep`) or after all possible twosteps of two egos (`ts_alternatives_twostep`)
#' @export
#' @references
#' \insertRef{ripley2022manual}{RsienaTwoStep}
#'
#' @seealso [`ts_select()`], [`ts_sims()`]
#' @examples
#' ts_alternatives_ministep(net=net1, ego=3)
ts_alternatives_ministep <- function(net, ego) {

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

#' @rdname ts_alternatives_ministep
#' @export
ts_alternatives_simstep <- function(net, ego) {
  #first ministep
  results <- ts_alternatives_ministep(net=net, ego=ego)
  #second ministep
  results2 <- lapply(results, ts_alternatives_ministep, ego=ego)
  results2 <- unlist(results2, recursive = FALSE)
  #remove duplicates
  results2 <- results2[!duplicated(results2)]
  return(results2)
}



#' @rdname ts_alternatives_ministep
#' @export
ts_alternatives_twostep <- function(net, dist1=NULL, dist2=NULL, modet1="degree", modet2="degree") {
  #dist1: minimal distance between ego1 and ego2 at time1
  #modet1: distance at time1 based on outdegree, indegree or degree
  #dist2: minimal distance between ego1 and ego2 at time2
  #modet2: distance at time2 based on outdegree, indegree or degree

  if (is.null(dist1) & is.null(dist2)) { #complete random selection
    egos <- ts_select(net=net, steps=2)
    results <- ts_alternatives_ministep(net=net, ego=egos[1])
    results2 <- lapply(results, ts_alternatives_ministep, ego=egos[2])
    results2 <- unlist(results2, recursive = FALSE)
  }

  if (!is.null(dist1) & is.null(dist2)) { #ego1 and ego2 should be connected at t1
    repeat {
      egos <- ts_select(net=net, steps=2) #need to optimize the ts_select function (only select ego1 among those who have at least alters that fit condition, than sample ego2 among those alters)
      dist_t1 <- ts_geodist(net=net, ego1=egos[1], ego2=egos[2], degree=modet1)
      if (dist_t1<=dist1) { #check if connected at t1, if not sample again, if true construct alternative nets
        results <- ts_alternatives_ministep(net=net, ego=egos[1])
        results2 <- lapply(results, ts_alternatives_ministep, ego=egos[2])
        results2 <- unlist(results2, recursive = FALSE)
        break
      }
    }
  }

  if (!is.null(dist1) & !is.null(dist2)) { #ego1 and ego2 should be connected at t1 OR t2
    #if already connected at t1, easy:
    repeat {
      egos <- ts_select(net=net, steps=2)
      dist_t1 <- ts_geodist(net=net, ego1=egos[1], ego2=egos[2], degree=modet1)
      #we will contstruct all possible nets after the twosteps, even if we do not know beforehand if we need them.
      results <- ts_alternatives_ministep(net=net, ego=egos[1])
      results2 <- lapply(results, ts_alternatives_ministep, ego=egos[2])

      if (dist_t1<=dist1) { #if connected at t1, we can simply use all constructed alternative nets.
        results2 <- unlist(results2, recursive = FALSE)
        break
      } else {
        #only keep those alternative networks were distance between ego1 and 2 smaller than dist2
        keep <- rapply(results2, f=ts_geodist, ego1=egos[1], ego2=egos[2], degree=modet2, how=c("list"))
        results2 <- unlist(results2, recursive=FALSE)
        results2 <- results2[unlist(keep, recursive=FALSE)<=dist2]
      }
      if (length(results2)>=1) break #check if we have at least one alternative network
    }
  }
  return(list(egos,results2))
}

#' Distance between alters
#'
#' Calculates the distance, or path length, between two actors part of the same network.
#'
#' @details This function is a wrapper around the function [`sna::geodist`].
#'
#' @return numeric vector of length one, with the distance between the two egos. If egos are not connected `ts_geodist` returns `Inf`.
#'
#' @param net matrix, the adjacency matrix representing the relations between actors. Valid values are 0 and 1.
#' @param ego1 numeric, value indicating ego1 (row number of net)
#' @param ego2 numeric, value indicating ego2 (row number of net)
#' @param degree string, the type of path to be considered: "`outdegree`"; "`indegree`"; "`degree`".
#' @export
#' @examples
#' ts_geodist(net=net1, ego1=6, ego2=8, degree="outdegree")
#' ts_geodist(net=net1, ego1=8, ego2=6, degree="outdegree")
#' ts_geodist(net=net1, ego1=8, ego2=6, degree="degree")
#' @export
ts_geodist <- function(net, ego1, ego2, degree="degree"){
  if (degree=="outdegree") {
    sna::geodist(dat=net, inf.replace=Inf, count.paths=FALSE, predecessors=FALSE, ignore.eval=TRUE, na.omit=TRUE)$gdist[ego1,ego2]
  } else if (degree=="indegree") {
    sna::geodist(dat=t(net), inf.replace=Inf, count.paths=FALSE, predecessors=FALSE, ignore.eval=TRUE, na.omit=TRUE)$gdist[ego1,ego2]
  } else {
    sna::geodist(dat=net + t(net), inf.replace=Inf, count.paths=FALSE, predecessors=FALSE, ignore.eval=TRUE, na.omit=TRUE)$gdist[ego1,ego2]
  }
}

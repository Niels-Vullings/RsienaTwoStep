ts_eval2 <- function(net, ego, statistics, ccovar=NULL, parameters) {
  # calculuate value of evaluation function

  s <- 0
  for (j in 1:length(statistics)) {
    if (length(statistics[[j]])==1) s <- s + parameters[j] * statistics[[j]](net, ego)
    if (length(statistics[[j]])==2) s <- s + parameters[j] * statistics[[j]][[1]](net, ego, ccovar[,statistics[[j]][[2]]])
  }


  # s <- foreach(j = 1:length(statistics), .combine = 'c') %dopar% {
  #   if (length(statistics[[j]])==1) { s <- parameters[j] * statistics[[j]](net, ego) }
  #   if (length(statistics[[j]])==2) { s <- parameters[j] * statistics[[j]][[1]](net, ego, ccovar[,statistics[[j]][[2]]]) }
  #   s
  #   }
  # s <- sum(s)
  #
  return(s)
}

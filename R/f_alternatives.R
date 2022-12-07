#' Title
#'
#' @param net
#' @param ego
#'
#' @return
#' @export
#'
#' @examples
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

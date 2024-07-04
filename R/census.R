#' @title Network census
#'
#' @description These functions calculate characteristics of the simulated
#'   networks.
#'
#' @details For examples on how to use these statistics see:
#'   vignette("Introduction_RsienaTwoStep").
#'
#' @family networkcensus
#' @param sims list, a list of (simulated) networks, the adjacency matrices
#'   representing the relations between actors. Valid values are 0 and 1. These
#'   simulated networks can be saved in objects that result from running
#'   [`ts_sims()`], [`ts_estim()`] or [`RSiena::siena07()`].
#' @param simtype string, name of the simulation type used (e.g. *ministep*,
#'   *twostep*).
#' @param forplot logical, if set to `FALSE` a dataframe is returned with in the
#'   column the network characteristic and each row represents a simulation
#'   outcome. If set to `TRUE` this dataframe is manipulated a bit, so that each
#'   row represents one specific network characteristic for each simulation
#'   outcome, this is useful for plotting.
#' @param cov numeric, covariate scores
#' @param mode Character string, “out” for out-degree, “in” for in-degree or “total” for the sum of the two. “all” is a synonym of “total”.
#' @param ans Results of class sienaFit, produced by a call to
#'   [`RSiena::siena07()`]
#' @importFrom foreach %dopar%
#' @importFrom iterators icount
#' @return `data.frame`
#' @seealso [`RSiena::sienaGOF()`], [`RSiena::sienaGOF-auxiliary()`]
#' @examples
#' \dontrun{
#' results_ministep <- ts_sims(net=net1, rate=5, statistics=list(ts_degree, ts_recip),
#' parameters=c(-3,1))
#' results_twostep <- ts_sims(net=net1, rate=5, statistics=list(ts_degree, ts_recip),
#' parameters=c(-3,1), p2step=1)
#'
#' dts_ms <- ts_dyads(sims=results_ministep, simtype="ministep")
#' dts_ts <- ts_dyads(sims=results_twostep, simtype="twostep")
#'
#' df <- rbind(dts_ms, dts_ts)
#' p <- ggplot(df, aes(x=x, y=y, fill=type)) +
#'  geom_violin(position=position_dodge(1)) +
#'  stat_summary(fun = mean,
#'               geom = "errorbar",
#'               fun.max = function(x) mean(x) + sd(x),
#'               fun.min = function(x) mean(x) - sd(x),
#'               width=.1,
#'               color="red", position=position_dodge(1)) +
#'  stat_summary(fun = mean,
#'               geom = "point",
#'               color="red", position=position_dodge(1))
#'
#'p
#'}
#' @export
ts_dyads <- function(sims, net1, simtype="notypespecified", forplot=TRUE) {

  if(is.list(sims) == TRUE){

    nsims <- length(sims)
  } else{

    nsims <- 1
  }

  #combine results of dyad.census
  df <- foreach::foreach(1:nsims, i=iterators::icount(), .combine="rbind") %dopar% {

    if(is.list(sims) == TRUE){
      net2 <- sims[[i]]
    } else{
      net2 <- sims
    }

    # t1 <- t1 #+ 1
    diag(net1) = NA #exclude the diagonal, not relevant data
    # t2 <- t2
    diag(net2) = NA #exclude the diagonal, not relevant data

    flips <- net1 + t(net2) - net2 #subtract s502 to ensure that Mutual ties are not included

    jumpst1 <- net1 + t(net1) #t1 plus its transpose lead to a value of 2 for mutual ties
    stablet1 <- jumpst1 #use this for stable assymetric ties
    jumpst1[lower.tri(jumpst1)] <- NA #remove duplicate ties

    jumpst2 <- net2 + t(net2)
    stablet2 <- jumpst2 #use this for stable assymetric ties
    jumpst2[lower.tri(jumpst2)] <- NA #remove duplicate ties

    #Dyad combinations
    stable00 <- as.data.frame(which(jumpst1 == 0 & jumpst2 == 0, arr.ind = TRUE)) #Null at t1 and Null at t2
    stable01 <- as.data.frame(which(flips == 0 & stablet1 == 1 & stablet2 == 1, arr.ind = TRUE)) #Assymetric at t1 and t2 and no flip
    stable11 <- as.data.frame(which(jumpst1 == 2 & jumpst2 == 2, arr.ind = TRUE)) #Mutual at t1 and Mutual at t2

    Null_Assym <- as.data.frame(which(jumpst1 == 0 & stablet2 == 1, arr.ind = TRUE))#Null -> assym
    Assym_Null <- as.data.frame(which(flips == 0 & stablet1 == 1 & stablet2 == 0, arr.ind = TRUE))#Assum -> Null
    Assym_Mut <- as.data.frame(which(stablet1 == 1 & jumpst2 == 2, arr.ind = TRUE))#Assym -> Mutual
    Mut_Assym <- as.data.frame(which(jumpst1 == 2 & jumpst2 == 1, arr.ind = TRUE))#Assym -> Mutual

    flip <- as.data.frame(which(flips == 2 & stablet1 != 2 & stablet2 != 2, arr.ind = TRUE)) #flips, so 01 at T1 and 10 at T2
    jump02 <- as.data.frame(which(jumpst1 == 0 & jumpst2 == 2, arr.ind = TRUE)) #Null jump, from Null to Mutual
    jump20 <- as.data.frame(which(jumpst1 == 2 & jumpst2 == 0, arr.ind = TRUE)) #Mutual jump, from Mutual to Null


    table <- cbind(nrow(stable00),
                   nrow(stable01),
                   nrow(stable11),
                   nrow(Null_Assym),
                   nrow(Assym_Null),
                   nrow(Assym_Mut),
                   nrow(Mut_Assym),
                   nrow(flip),
                   nrow(jump02),
                   nrow(jump20)) #bind the rownumbers to determine the census of the tie type

    colnames(table) = c("Null > Null",
                        "Assym > Assym",
                        "Mutual > Mutual",
                        "Null > Assym",
                        "Assym > Null",
                        "Assym > Mutual",
                        "Mutual > Assym",
                        "Tie flip",
                        "Null > Mutual",
                        "Mutual > Null") # Rename the columns to match the dyad change type

    return(data.frame(simnet = i, table))
  }
  df<- as.data.frame(df)
  df$type <- simtype

  if (forplot == TRUE) { #bit clumsy

    df <- df %>% tidyr::pivot_longer(cols = !c(simnet,type),
                                     names_to = "x",
                                     values_to = "y")
  }

  return(df)
}

#' @rdname ts_dyads
#' @export
ts_triads <- function(sims, net1, simtype="notypespecified", forplot=TRUE) {

  if(is.list(sims) == TRUE){
    nsims <- length(sims)
  } else{
    nsims <- 1
  }

  df <- foreach::foreach(1:nsims, x=iterators::icount(), .combine="rbind") %dopar% {

    if(is.list(sims) == TRUE){
      net2 <- sims[[x]]
    } else{
      net2 <- sims
    }

    #Davis & Leinhardt triad classification for identifying triads
    triads <- c("X003", "X012", "X102", "X021D",  "X021U",  "X021C",  "X111D",  "X111U",  "X030T",  "X030C",   "X201",  "X120D",  "X120U",  "X120C",   "X210",   "X300")

    # Create a dataframe(df) that determines a triad census for each possible triad in the matrix
    df <- foreach::foreach(a1=1:nrow(net1), i=iterators::icount() , .combine="rbind") %:%
      foreach::foreach(a2=1:nrow(net1), j=iterators::icount() , .combine="rbind") %:%
      foreach::foreach(a3=1:nrow(net1),  k=iterators::icount() , .combine="rbind") %do% {

        if (i>j & j>k ) data.frame(i=i, j=j, k=k, #(i>j & j>k ) determines with or without repititions [current = NO REP]
                                   t1_ij = as.character(net1[a1,a2]), t1_ji = as.character(net1[a2,a1]),
                                   t1_ik = as.character(net1[a1,a3]),t1_ki = as.character(net1[a3,a1]),
                                   t1_jk = as.character(net1[a2,a3]),t1_kj = as.character(net1[a3,a2]), #Configuration of T1 triad
                                   typeT1 = triads[which(sna::triad.census(net1[c(a1,a2,a3),c(a1,a2,a3) ]) == 1) ], # Census of specific triad, based on the values that are in a1,a2,a3 for T1
                                   t2_ij = as.character(net2[a1,a2]), t2_ji = as.character(net2[a2,a1]),
                                   t2_ik = as.character(net2[a1,a3]),t2_ki = as.character(net2[a3,a1]),
                                   t2_jk = as.character(net2[a2,a3]),t2_kj = as.character(net2[a3,a2]), # Configuration of T2 triad
                                   typeT2 = triads[which(sna::triad.census(net2[c(a1,a2,a3),c(a1,a2,a3) ]) == 1) ], # Census of specific triad, based on the values that are in a1,a2,a3 for T2
                                   tie_change = sum(c(net1[a1,a2],net1[a2,a1],net1[a1,a3],net1[a3,a1],net1[a2,a3],net1[a3,a2]) != c(net2[a1,a2],net2[a2,a1],net2[a1,a3],net2[a3,a1],net2[a2,a3],net2[a3,a2])))
      }

    df$name <- paste0(df$i,".",df$j,".",df$k) #name is each actor in triad [actor1.actor2.actor3]

    df <- subset(df, select=c(tie_change)) # Subset to remove redundant information

    df <- data.frame(simnet = x, table(df))
    df <- tidyr::pivot_wider(df, names_from = tie_change, values_from = Freq)

  }
  df <- as.data.frame(df)
  df$type <- simtype

  if (forplot == TRUE) {

    df <- df %>% tidyr::pivot_longer(cols = !c(simnet,type),
                                     names_to = "x",
                                     values_to = "y")
  }

  return(df)
}

#' @rdname ts_dyads
#' @export
ts_nacf <- function(sims, simtype="notypespecified", forplot=TRUE, cov) {
  nsims <- length(sims)
  df <- foreach::foreach(1:nsims, i=iterators::icount(), .combine="rbind") %dopar% {
    sna::nacf(sims[[i]], cov, type = "moran", neighborhood.type = "out", demean = TRUE)[2]
  }
  df <- as.data.frame(df)
  df$type <- simtype
  return(df)
}



#' @rdname ts_dyads
#' @export
ts_degreecount <- function(sims, mode="out", simtype="notypespecified", forplot=TRUE) {
  nsims <- length(sims)
  nactors <- dim(sims[[1]])[1]
  mattemp <- matrix(0, nrow=nsims, ncol=nactors)
  for (j in 1:nsims) {
    t1 <- table(igraph::degree(igraph::graph_from_adjacency_matrix(sims[[j]]), mode = mode))
    degrees <- as.numeric(names(t1))
    freqs <- as.numeric(t1)
    for (i in 1:length(t1)) {
      mattemp[j, degrees[i] + 1] <- freqs[i]
    }
  }
  df <- as.data.frame(mattemp)
  names(df) <- paste0("deg", 1:nactors)
  df <- df[,colSums(df)!=0]

  if (forplot) {
    degrees <- names(df)
    dflist <- list()
    for (i in 1:length(degrees)) {
      dflist[[i]] <- df
    }
    df <- do.call(rbind, dflist)
    df$x <- rep(degrees, each=nsims)
    df$y <- NA

    for (i in 1:length(degrees)) {
      df$y[df$x==degrees[i]] <- df[,degrees[i]][df$x==degrees[i]]
    }
    df <- df[,c("x", "y")]
  }

  df$x <- factor(df$x, levels=unique(df$x))
  df$type <- simtype
  return(df)
}

#' @rdname ts_dyads
#' @export
ts_rsienanets <- function(ans) {
  n <- dim(ans$f$Data1$depvars$mynet[, , 1])[1]
  sims <- foreach(i = 1:length(ans$sims)) %dopar% {
    edges <- ans$sims[[i]][[1]][[1]][[1]]
    # create empty adjacency matrix
    adj <- matrix(0, n, n)
    # put edge values in desired places
    adj[edges[, 1:2]] <- edges[, 3]
    adj
  }
  return(sims)
}


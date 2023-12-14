#' @title Network census
#'
#' @description These functions calculate characteristics of the simulated
#'   networks. For now, only a dyad census and a triad census are implemented.
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
ts_dyads <- function(sims, simtype="notypespecified", forplot=TRUE) {
  nsims <- length(sims)
  #combine results of dyad.census
  df <- foreach::foreach(1:nsims, i=iterators::icount(), .combine="rbind") %dopar% {
    sna::dyad.census(sims[[i]])
  }
  df <- as.data.frame(df)

  if (forplot) { #bit clumsy
    df <- rbind(df, df, df)
    df$x <- rep(c("mut", "asym", "null"), each=nsims)
    df$y <- NA
    df$y[df$x=="mut"] <- df$Mut[df$x=="mut"]
    df$y[df$x=="asym"] <- df$Asym[df$x=="asym"]
    df$y[df$x=="null"] <- df$Null[df$x=="null"]
  }
  df$type <- simtype
  return(df)
}

#' @rdname ts_dyads
#' @export
ts_triads <- function(sims, simtype="notypespecified", forplot=TRUE) {
  nsims <- length(sims)
  df <- foreach::foreach(1:nsims, i=iterators::icount(), .combine="rbind") %dopar% {
    sna::triad.census(sims[[i]])
  }
  df <- as.data.frame(df)

  if (forplot) {
    triads <- names(df)
    dflist <- list()
    for (i in 1:length(triads)) {
      dflist[[i]] <- df
    }
    df <- do.call(rbind, dflist)
    df$x <- rep(triads, each=nsims)
    df$y <- NA

    for (i in 1:length(triads)) {
      df$y[df$x==triads[i]] <- df[,triads[i]][df$x==triads[i]]
    }
  }
  df$type <- simtype
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


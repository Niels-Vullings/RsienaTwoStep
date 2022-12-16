#' @title Network census
#'
#' @description These functions calculate characteristics of the simulated networks. For now, only a dyad census and a triad census are implemented.
#'
#' @details For examples on how to use these statistics see: vignette("ABMministep-vs-ABMtwostep").
#'
#' @family networkcensus
#' @param sims list, the result of `f_sims():` the adjacency matrix representing the relations between actors. Valid values are 0 and 1.
#' @param simtype string, name of the simulation type used (e.g. *ministep*, *twostep*).
#' @param forplot logical, if set to `FALSE` a dataframe is returned with in the column the network characteristic and each row represents a simulation outcome.
#' If set to `TRUE` this dataframe is manipulated a bit, so that each row represents one specific network characteristic for each simulation outcome, this is useful for plotting.
#'
#' @return `data.frame`
#'
#' @examples
#' \dontrun{
#' results_ministep <- f_sims(net=net1, rate=5, statistics=list(f_degree, f_recip),
#' parameters=c(-3,1))
#' results_twostep <- f_sims(net=net1, rate=5, statistics=list(f_degree, f_recip),
#' parameters=c(-3,1), p2step=1)
#'
#' df_ms <- f_dyads(sims=results_ministep, simtype="ministep")
#' df_ts <- f_dyads(sims=results_twostep, simtype="twostep")
#'
#' df <- rbind(df_ms, df_ts)
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
f_dyads <- function(sims, simtype="notypespecified", forplot=TRUE) {
  nsims <- length(sims)
  #combine results of dyad.census
  df <- foreach(1:nsims, i=icount(), .combine="rbind") %dopar% {
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

#' @rdname f_dyads
#' @export
f_triads <- function(sims, simtype="notypespecified", forplot=TRUE) {
  nsims <- length(sims)
  df <- foreach(1:nsims, i=icount(), .combine="rbind") %dopar% {
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

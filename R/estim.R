#' Estimation of parameters via Robbins-Monro algorithm
#'
#' @description `ts_estim` is the workhorse function of the package
#'   `RsienaTwoStep`. For details on the simulation function see [`ts_sims()`].
#'   This function aims to estimate the model parameters according to the
#'   Robbins-Monro algorithm as described in
#'   \insertCite{snijders2001statistical}{RsienaTwoStep}.
#' @details For examples on how to use `ts_estim()` see:
#'   `vignette("1.Introduction_RsienaTwoStep", package="RsienaTwoStep")` or
#'   the [package website](https://jochemtolsma.github.io/RsienaTwoStep/).
#'   Before you set `parallel` to TRUE make sure to set-up a cluster with the
#'   package `doParallel` (see `Examples`).
#' - `p2step==c(1,0,0)`: ministep
#' - `p2step==c(0,1,0)` & dist1==NULL & dist2==NULL: twostep-simultaneity
#' - `p2step==c(0,1,0)` & dist1!=NULL & dist2==NULL: twostep-strict coordination
#' - `p2step==c(0,1,0)` & dist1!=NULL & dist2!=NULL: twostep-weak coordination
#' - `p2step==c(0,0,1)`: two-ministeps
#'
#'
#' @param ans Results of class sienaFit, produced by a call to
#'   [`RSiena::siena07()`]
#' @param mydata Siena data object created by a call to
#'   [`RSiena::sienaDataCreate()`]
#' @param myeff Siena effects object created by a call to
#'   [`RSiena::getEffects()`]
#' @param startvalues if not provided manually, taken from results of
#'   `ans$theta`, or if `ans` is not provided from
#'   `summary(myeff)$initialValue`.
#' @param net1 adjacency matrix, the adjacency matrix representing the relations
#'   between actors at Time=1. Valid values are 0 and 1. If not provided
#'   manually, retrieved from `ans` or `mydata`.
#' @param net2 adjacency matrix, the adjacency matrix representing the relations
#'   between actors at Time=2. Valid values are 0 and 1. If not provided
#'   manually, retrieved from `ans` or `mydata`.
#' @param ccovar data frame with named time-constant covariates. If not provided
#'   manually, retrieved from `ans` or `mydata`.
#' @param statistics list of statistics of `RsienaTwoStep`, see: [`ts_degree()`]
#'   and DETAILS. If not provided manually, retrieved from `ans` or `myeff`.
#' @param b numeric between 0.1 and 1 (default =0.5), used in Robbins-Monro
#'   algorithm.
#' @param conv numeric. Robbins-Monro algorithm stops if the mean deviation of
#'   parameters after each update steps become smaller than `conv`
#' @param nite number of iterations phase2 (actual estimation of parameters)
#' @param itef1 number of total iterations phase1 (To get to a crude Jacobian
#'   matrix, which can be used in phase2) )
#' @param itef3 number of total iterations phase3 (To estimate SE). This phase
#'   will take the longest.
#' @param p2step numeric vector of length 3, setting the ratio of ministep,
#'   twostep and twoministeps. This parameter is passed to [`ts_sims()`].
#' @param dist1 numeric, minimal path length between ego1 and ego2 at time1 in
#'   order to be allowed to start a cooperation. If `NULL` all dyads are allowed
#'   to start a cooperation. This parameter is passed to [`ts_sims()`].
#' @param dist2 numeric, minimal path length between ego1 and ego2 at time2 in
#'   order for twostep to be counted as cooperation. This parameter is passed to
#'   [`ts_sims()`].
#' @param modet1 string, indicating the type of ties being evaluated at time1.
#'   "*degree*" considers all ties as undirected. "*outdegree*" only allows
#'   directed paths starting from ego1 and ending at ego2. "*indegree*" only
#'   allows directed paths starting from ego2 and ending at ego2. This parameter
#'   is passed to [`ts_sims()`].
#' @param modet2 string, indicating the type of ties being evaluated at time2.
#'   "*degree*" considers all ties as undirected. "*outdegree*" only allows
#'   directed paths starting from ego1 and ending at ego2. "*indegree*" only
#'   allows directed paths starting from ego2 and ending at ego2. This parameter
#'   is passed to [`ts_sims()`].
#' @param verbose, TRUE/FALSE. If set to true it shows the iteration steps and
#' some results.
#' @param parallel Boolean. TRUE/FALSE.
#' @param returnDeps Boolean. If `TRUE` the simulated dependent variables
#'   (networks, behaviour) of phase3 will be returned.
#' @param phase1 TRUE/FALSE, If False no Jacobian matrix is calculated.
#' @param phase3 TRUE/FALSE, if FALSE no SE are calculated
#' @references \insertRef{ripley2022manual}{RsienaTwoStep}
#'   \insertRef{snijders2001statistical}{RsienaTwoStep}
#' @returns If `phase3 = FALSE ` A dataframe of estimated parameters. The last row are the final
#'   solutions of the Robbins Monro algorithm. If `Phase 3 = TRUE` a list with first list element the data frame of estimated parameters and the second element the results of phase 3.
#' @export
#' @seealso [`ts_sims()`],
#' [`ts_degree()`],
#'   [`RSiena::siena07()`]
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @importFrom foreach foreach
#' @importFrom Rdpack reprompt
#' @importFrom RSiena sienaDependent
#' @importFrom RSiena sienaDataCreate
#' @importFrom RSiena getEffects
#' @importFrom stats cov
#' @importFrom stats runif
#'
#' @examples
#' \dontrun{
#' # It is good practice to check if `RsienaTwoStep` uses the same
#' # target values as `RSiena`, which can be retrieved like: `ans$targets`.
#' ts_targets(net1 = s501, net2 = s502, statistics = list(ts_degree, ts_recip))
#' # In normal use case you do not want to estimate the Jacobian matrix yourself
#' # for phase 1. We simply use the matrix from `RSiena`. Except if you really
#' # expect different outcomes from the twostep model and not using a
#' # Dinv matrix at all does not work.
#' # Phase 1 only
#' jac <- ts_phase1(net1 = s501,
#'  net2 = s502,
#'  statistics = list(ts_degree, ts_recip),
#'  ccovar = NULL,
#'  itef1 = 30)
#' dinv <- solve(jac)
#' #estimate without the use of RSiena
#' ts_estim(net1 = s501,
#' net2 = s502,
#' statistics = list(ts_degree, ts_recip),
#' nite = 30,
#'  phase1 = FALSE)
#'  #estimate without the use of RSiena
#' ts_estim(net1 = s501,
#' net2 = s502,
#' statistics = list(ts_degree, ts_recip),
#' nite = 30,
#' itef1 = 10,
#' phase1 = TRUE)
#' #estimate without the use of RSiena
#' ts_estim(net1 = s501,
#'  net2 = s502,
#'  statistics = list(ts_degree, ts_recip),
#'  nite = 30,
#'  itef1 = 10,
#'  phase1 = TRUE,
#'  itef3 = 10,
#'  phase3 = TRUE)
#' # Phase 3 only
#' startvalues <- c(5.5, -2.2, 2.4)
#' stat <- list(ts_degree, ts_recip)
#' ts_phase3(startvalues = startvalues,
#' net1 = s501,
#' statistics = stat,
#' itef3 = 10,
#' verbose = TRUE)
#' library(RSiena)
#' mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
#' mydata <- sienaDataCreate(mynet)
#' myalgorithm <- sienaAlgorithmCreate(cond=FALSE)
#' #toggle set conditional to retrieve the rate parameter in theta!
#' myeff <- getEffects(mydata)
#' ts_estim(mydata = mydata, myeff = myeff)
#' #warning, this may take a while!
#' ts_estim(mydata = mydata, myeff = myeff, phase3 = TRUE)
#' ans1 <- siena07(myalgorithm, data=mydata, effects=myeff)
#' ts_estim(ans1)
#' }
#'
#' @export
ts_estim <- function(ans = NULL,
                     mydata = NULL,
                     myeff = NULL,
                     startvalues = NULL,
                     net1 = NULL,
                     net2 = NULL,
                     ccovar = NULL,
                     statistics = NULL,
                     b = 0.5,
                     conv = 0.001,
                     nite = 1000,
                     itef1 = 50,
                     itef3 = 1000,
                     p2step = c(1, 0, 0),
                     dist1 = NULL,
                     dist2 = NULL,
                     modet1 = "degree",
                     modet2 = "degree",
                     verbose = TRUE,
                     parallel = FALSE,
                     returnDeps = FALSE,
                     phase1 = TRUE,
                     phase3 = TRUE) {

  ### initializing function
{
  # starting networks
  nets <- ts_netprep(ans=ans, mydata=mydata, net1=net1, net2=net2)
  net1 <- nets$net1
  net2 <- nets$net2

  # ccovar (and deps)
  ccovar <- ts_dataprep(ans= ans, mydata= mydata, ccovar=ccovar)

  # included statistics
  stats <- ts_statprep(ans = ans, myeff = myeff, statistics = statistics)
  statistics <- stats$statistics
  namesstatistics <- stats$namesstatistics
  ratebeh <- stats$ratebeh
  netstats <- stats$netstats

  #startvalues (afster ts_statprep)
  startvalues <- ts_startprep(ans=ans, myeff=myeff, net1 = net1, net2 = net2,
                              ccovar = ccovar, startvalues= startvalues,
                              namesstatistics = namesstatistics)

}


  ### Phase 1
  {
  # INVERSE OF jacobi matrix for phase1
  if (!is.null(ans)) {
    dinv <-
      solve(ans$dfra) # or use ans$dinv ?? not sure what difference is.
    # If we have ans, we want dfra from phase 3 as input not dfra1 from phase 1!!
    print("phase 1 from ans")
  } else if (phase1) {
    jacob <- ts_phase1(
      mydata = mydata,
      myeff = myeff,
      net1 = net1,
      net2 = net2,
      ccovar = ccovar,
      statistics = statistics,
      itef1 = itef1,
      p2step = p2step,
      dist1 = dist1,
      dist2 = dist2,
      modet1 = modet1,
      modet2 = modet2,
      verbose = verbose,
      parallel = parallel
    )
    dinv <- solve(jacob)
    if (verbose) print("phase 1 done, start phase2")
  } else {
    dinv <- NULL
    if (verbose) print("skipped phase 1, no dinv in Robbins-Monro algorithm used.")
  }

  print("targets")
  # target Z statistics
  targets <- ts_targets(
    ans = ans,
    mydata = mydata,
    myeff = myeff,
    net1 = net1,
    net2 = net2,
    ccovar = ccovar,
    statistics = statistics
  )
  print(targets)




}

  ### Phase 2
  if (verbose) {
    print("startvalues: ")
    print(startvalues)
  }

    ##initialization of RM
  if (verbose) print("start of phase 2")
  x <- startvalues
  dif <- Inf
  ite <- 0
  ite2 <- 1 # try to keep an constant as long as the sequence sn has not crossed the observed values
  r <- x
  r_sub <- x # for the sub-phase
  update <- 0
  diagdinv <- diag(dinv)
  AC <- 0
    ##end initialization of RM

  while (max(dif) > conv & ite < nite) {
    ite <- ite + 1 # number of iterations
    a <- (1 / ite2) ^ b
    update_old <- update

    if (!parallel) {

    # estimate model
    sims1 <- ts_sims(
      startvalues = x,
      net1 = net1,
      ccovar = ccovar,
      statistics = statistics,
      nsims = 1,
      p2step = p2step,
      dist1 = dist1,
      dist2 = dist2,
      modet1 = modet1,
      modet2 = modet2,
      chain = FALSE,
      verbose = FALSE,
    )

    # calculate the observed statistics
    if (!length(ratebeh)>0) {
    Z <- ts_targets(
      net1 = net1,
      net2 = sims1[[1]],
      ccovar = ccovar,
      statistics = statistics
    )
    } else {
      ccovar_n <- ccovar
      ccovar_n[,2] <- sims1[[1]]$beh_n
      Z <- ts_targets(
        net1 = net1,
        net2 = sims1[[1]]$net_n,
        ccovar = ccovar_n,
        statistics = statistics
      )
    }

    }

    if (parallel) {
      cores <- parallel::detectCores() - 1
      Zs <- foreach(i = 1:cores, .combine="rbind") %dopar% {
        # estimate model
        sims1 <- ts_sims(
          startvalues = x,
          net1 = net1,
          ccovar = ccovar,
          statistics = statistics,
          nsims = 1,
          p2step = p2step,
          dist1 = dist1,
          dist2 = dist2,
          modet1 = modet1,
          modet2 = modet2,
          chain = FALSE,
          verbose = FALSE,
        )

        # calculate the observed statistics
        if (!length(ratebeh)>0) {
          Z <- ts_targets(
            net1 = net1,
            net2 = sims1[[1]]$net_n,
            ccovar = ccovar,
            statistics = statistics
          )
        } else {
          ccovar_n <- ccovar
          ccovar_n[,2] <- sims1[[1]]$beh_n
          Z <- ts_targets(
            net1 = net1,
            net2 = sims1[[1]]$net_n,
            ccovar = ccovar_n,
            statistics = statistics
          )
        }

      }
      Z <- colMeans(Zs)
    }

    #deviations
    update <- Z - targets

    if (!is.null(dinv)) {
      x <-
        x - pmax(pmin(diagdinv * a * update, 1), -1) #avoiding big jumps by pmin/pmax
    } else {
      x <- x - pmax(pmin((1 / nrow(net1)) * a * update, 1),-1) # if ans is not provided divide by nnodes?
    }
    # interfere with the updating. Is this allowed?? I guess avoiding big jumps is better.
     if (x[1] < 0.5)
       x[1] <- 1 # avoid negative rates
     if (x[2] < -5)
       x[2] <- -5 # avoid too negative degrees
     x[3:length(x)][x[3:length(x)] < -10] <- -9
     x[3:length(x)][x[3:length(x)] > 10] <- 9


    #check if we need to update ite2, that is make a smaller
    #if (sum(abs(sign(update) - sign(update_old)) == 2) > 1) {
    # the mean value of the product of old/new deviations is 0
    #if (mean(update * update_old) < 0) {
    # all products of deviations are negative
    #AC <- AC + (update * update_old)
    #if (max(AC) <= 0) {
      # check with TS if this is correct

     if (mean(update * update_old) < 0) {
      r_sub <- rbind(r_sub, x) # save results
      ite2 <- ite2 + 1
      #AC <- 0
      x <-
        colMeans(r_sub) # for new subphase take average of sequence as new starting value
      r_sub <- x
    } else {
      r_sub <- rbind(r_sub, x) # save results
    }

    r <- rbind(r, x) # save results



    if ((ite %% 3) == 0)
      dif <- (abs(x - r[ite - 2,])) # TO DO: check how is implemented
    if (verbose) {
      print(x)
      print(paste0("ite: ", ite))
      print(paste0("ite2: ", ite2))
    }
  }

  colnames(r) <- namesstatistics

  ### Phase 3 and return
  if (!phase3) {
    return(r)
  } else {
    estim <- r[nrow(r),]
    phase3 <-
      ts_phase3(
        net1 = net1,
        net2 = net2,
        ccovar = ccovar,
        statistics = statistics,
        startvalues = estim,
        itef3 = itef3,
        p2step = p2step,
        dist1 = dist1,
        dist2 = dist2,
        modet1 = modet1,
        modet2 = modet2,
        verbose = verbose,
        parallel = parallel,
        returnDeps = returnDeps
      )
    #SE <- sqrt(diag(phase3$covtheta))
    #df <- data.frame(estim = estim, SE = SE)
    #row.names(df) <- namesstatistics
    return(list(robmon=r, phase3=phase3))
  }
}

#' @rdname ts_estim
#' @export
ts_targets <-
  function(ans = NULL,
           mydata = NULL,
           myeff = NULL,
           net1 = NULL,
           net2 = NULL,
           ccovar = NULL,
           statistics = NULL) {

    if (!is.null(ans)) {
      targets <- ans$targets
      names(targets) <- ans$effects$shortName
      return(targets)
      stop()
    }

    ### initialize
    # starting networks
    nets <- ts_netprep(ans=ans, mydata=mydata, net1=net1, net2=net2)
    net1 <- nets$net1
    net2 <- nets$net2
    # included statistics
    stats <- ts_statprep(ans = ans, myeff = myeff, statistics = statistics)
    statistics <- stats$statistics
    namesstatistics <- stats$namesstatistics
    ratebeh <- stats$ratebeh
    netstats <- stats$netstats
    # ccovar (and deps)
    ccovar <- ts_dataprep(ans= ans, mydata= mydata, ccovar=ccovar)
    ### end

    # calculate the observed statistics
     Z <- rep(NA, length(namesstatistics) ) # empty vector


     for (j in seq_along(namesstatistics)) {
       if (namesstatistics[j] == "Rate") {
        Z[j] <- sum(abs(net2 - net1))
      } else if (namesstatistics[j] == "Rate beh") {
        Z[j] <- sum(abs(ccovar[, 2] - ccovar[, 1] ))
      } else {
        stat <- statistics[[which(names(statistics) == namesstatistics[j])]]

        if (length(stat) == 1) {
          #network statistic
          if (netstats[j]) {
          # Single-argument statistic function
          Z[j] <- sum(sapply(1:nrow(net2), function(i) stat(net = net2, ego = i)))
          } else {
            #behavior statistic
            Z[j] <- sum(sapply(1:nrow(net2), function(i) stat(beh = ccovar[,2], ego = i)))
          }
          # } else if (length(stat) == 2) {
        } else {
          #network statistic
          if (netstats[j]) {
          # Two-argument statistic function
          if (stat[[2]] %in% names(ccovar)) { #hence it is a ccovar
            Z[j] <- sum(sapply(1:nrow(net2), function(i) stat[[1]](net = net2, ego = i, cov = ccovar[, stat[[2]]])))
          } else { #it is a depvar
            depvar <- paste0(stat[[2]],".1") ##I really do not understand why in these statistics to calculate the target value the depvar at wave 1 is used.
            Z[j] <- sum(sapply(1:nrow(net2), function(i) stat[[1]](net = net2, ego = i, cov = ccovar[,depvar])))
          }
          } else {
            #behavior statistic
            cov <- NULL
            if (stat[[2]] %in% names(ccovar)) cov <- ccovar[, stat[[2]]]
            Z[j] <- sum(sapply(1:nrow(net2), function(i) stat[[1]](beh = ccovar[,2], net = net1, ego = i, cov = cov)))
            ##I really do not understand why in these statistics to calculate the target value the network at wave 1 is used.
          }
        }
      }
    }

    #correct target for cycle3
    names(Z) <- namesstatistics
    Z[namesstatistics=="cycle3"] <- Z[namesstatistics=="cycle3"]/3
    Z
  }

#' @rdname ts_estim
#' @export
ts_phase1 <- function(ans = NULL,
                      mydata = NULL,
                      myeff = NULL,
                      net1 = NULL,
                      net2 = NULL,
                      ccovar = NULL,
                      statistics = NULL,
                      itef1 = 100,
                      p2step = c(1, 0, 0),
                      dist1 = NULL,
                      dist2 = NULL,
                      modet1 = NULL,
                      modet2 = NULL,
                      verbose = TRUE,
                      parallel = FALSE) {

    # if siena07 is used to estimate a model we simply use the jacobi-matrix stored in that object
  if (!is.null(ans)) {
    dinv_f1 <- ans$dfra1
    return(dinv_f1)
    stop()
  }

  ### start initialization
  {
  #starting networks
  nets <- ts_netprep(ans=ans, mydata=mydata, net1=net1, net2=net2)
  net1 <- nets$net1
  net2 <- nets$net2

  # ccovar (and deps)
  ccovar <- ts_dataprep(ans= ans, mydata= mydata, ccovar=ccovar)

  # included statistics
  stats <- ts_statprep(ans = ans, myeff = myeff, statistics = statistics)
  statistics <- stats$statistics
  namesstatistics <- stats$namesstatistics
  ratebeh <- stats$ratebeh
  netstats <- stats$netstats

  #startvalues (afster ts_statprep)
  startvalues <- ts_startprep(ans=ans, myeff=myeff, net1 = net1, net2 = net2,
                              ccovar = ccovar, startvalues= startvalues,
                              namesstatistics = namesstatistics)



  }
  ### end of initialization


  # here we start with the actual phase 1.

  crn <- runif(itef1, 1 - 2^31, 2^31 - 1) # common random numbers
  pn <- length(startvalues)
  deviation <- rep(0.1, pn)
  deviation[names(startvalues) == "Rate"] <- startvalues[names(startvalues) == "Rate"] / 10
  deviation[names(startvalues) == "Rate beh"] <- startvalues[names(startvalues) == "Rate beh"] / 10

  if (verbose)
    print("start phase 1")

if (!parallel) {
  # to save the simulated z-scores
  res <-
    matrix(NA, nrow = itef1, ncol = pn + (pn * pn))


  for (i in 1:itef1) {
    if (verbose) print(i)

    # to save the simulated z-scores (for one itef1)
    rest <- rest2 <- rep(NA, pn)

    # Zbar
    set.seed(crn[i]) # is this just our crn
    sims1 <-
      ts_sim(
        net1 = net1,
        ccovar = ccovar,
        startvalues = startvalues,
        statistics = statistics,
        p2step = p2step,
        dist1 = dist1,
        dist2 = dist2,
        modet1 = modet1,
        modet2 = modet2
      )

    # calculate the observed statistics
    if (!length(ratebeh)>0) {
      rest <- ts_targets(
        net1 = net1,
        net2 = sims1$net_n,
        ccovar = ccovar,
        statistics = statistics
      )
    } else {
      ccovar_n <- ccovar
      ccovar_n[,2] <- sims1$beh_n
      rest <- ts_targets(
        net1 = net1,
        net2 = sims1$net_n,
        ccovar = ccovar_n,
        statistics = statistics
      )
    }

    # start chaning the startvalues for rate and statistics
    for (j in 1:length(startvalues)) {
      #print(j)
      #set new startvalues
      startvalues_n <- startvalues
      startvalues_n[j] <- startvalues_n[j] + deviation[j]

      #simulate network with new startvalues
      #set the crn
      set.seed(crn[i])
      sims1 <-
        ts_sim(
          net1 = net1,
          ccovar = ccovar,
          startvalues = startvalues_n,
          statistics = statistics,
          p2step = p2step,
          dist1 = dist1,
          dist2 = dist2,
          modet1 = modet1,
          modet2 = modet2
        )


      # calculate the observed statistics
      if (!length(ratebeh)>0) {
        rest2 <- ts_targets(
          net1 = net1,
          net2 = sims1$net_n,
          ccovar = ccovar,
          statistics = statistics
        )
      } else {
        ccovar_n <- ccovar
        ccovar_n[,2] <- sims1$beh_n
        rest2 <- ts_targets(
          net1 = net1,
          net2 = sims1$net_n,
          ccovar = ccovar_n,
          statistics = statistics
        )
      }

      rest <- c(rest, rest2)
    }
    res[i,] <- rest
  }
  # now we have all z-scores we can calculate the jacobi
}

if (parallel) {
    res <- foreach(i = 1: itef1, .combine="rbind", .verbose = verbose) %dopar% {

      # to save the simulated z-scores (for one itef1)
      rest <- rest2 <- rep(NA, pn)

      # Zbar
      set.seed(crn[i]) # is this just our crn
      sims1 <-
        ts_sim(
          net1 = net1,
          ccovar = ccovar,
          startvalues = startvalues,
          statistics = statistics,
          p2step = p2step,
          dist1 = dist1,
          dist2 = dist2,
          modet1 = modet1,
          modet2 = modet2
        )

      # calculate the observed statistics
      if (!length(ratebeh)>0) {
        rest <- ts_targets(
          net1 = net1,
          net2 = sims1$net_n,
          ccovar = ccovar,
          statistics = statistics
        )
      } else {
        ccovar_n <- ccovar
        ccovar_n[,2] <- sims1$beh_n
        rest <- ts_targets(
          net1 = net1,
          net2 = sims1$net_n,
          ccovar = ccovar_n,
          statistics = statistics
        )
      }

      # start chaning the startvalues for rate and statistics
      for (j in 1:length(startvalues)) {
        #print(j)
        #set new startvalues
        startvalues_n <- startvalues
        startvalues_n[j] <- startvalues_n[j] + deviation[j]

        #simulate network with new startvalues
        #set the crn
        set.seed(crn[i])
        sims1 <-
          ts_sim(
            net1 = net1,
            ccovar = ccovar,
            startvalues = startvalues_n,
            statistics = statistics,
            p2step = p2step,
            dist1 = dist1,
            dist2 = dist2,
            modet1 = modet1,
            modet2 = modet2
          )


        # calculate the observed statistics
        if (!length(ratebeh)>0) {
          rest2 <- ts_targets(
            net1 = net1,
            net2 = sims1$net_n,
            ccovar = ccovar,
            statistics = statistics
          )
        } else {
          ccovar_n <- ccovar
          ccovar_n[,2] <- sims1$beh_n
          rest2 <- ts_targets(
            net1 = net1,
            net2 = sims1$net_n,
            ccovar = ccovar_n,
            statistics = statistics
          )
        }
      rest <- c(rest, rest2)
      }
      rest
  }

  res_mat <- matrix(NA, nrow = pn, ncol = pn) # the jacobi matrix
  count <- pn + 1
  for (i in 1:pn) {
    for (j in 1:pn) {
      res_mat[i, j] <- mean((res[, count] - res[, j]) / deviation[i])
      count <- count + 1
    }
  }

  return(t(res_mat)) # perhaps return a list with more info, like the deviations and initial values and stuff
}

}



#' ts_phase1(ans=ans1)
#' ts_phase1(mydata=mydata, myeff = myeff, itef1=30)
#' ts_phase1(net1=s501, net2=s502, statistics=list(ts_degree, ts_recip), itef1=10)

#' @rdname ts_estim
#' @export
#' @returns list, containing: `estim`, `covtheta`, `tstat` `tconv.max`, `zdevs`, and if `returnDeps = TRUE` `simnets`
ts_phase3 <- function(ans = NULL,
                      mydata = NULL,
                      myeff = NULL,
                      net1 = NULL,
                      net2 = NULL,
                      ccovar = NULL,
                      statistics = NULL,
                      startvalues = NULL,
                      itef3 = 1000,
                      p2step = c(1, 0, 0),
                      dist1 = NULL,
                      dist2 = NULL,
                      modet1 = "degree",
                      modet2 = "degree",
                      verbose = TRUE,
                      parallel = FALSE,
                      returnDeps = FALSE) {
  if (verbose)
    print("start phase 3")

  ### start of initialization
  #starting networks
  nets <- ts_netprep(ans=ans, mydata=mydata, net1=net1, net2=net2)
  net1 <- nets$net1
  net2 <- nets$net2

  # ccovar (and deps)
  ccovar <- ts_dataprep(ans= ans, mydata= mydata, ccovar=ccovar)

  # included statistics
  stats <- ts_statprep(ans = ans, myeff = myeff, statistics = statistics)
  statistics <- stats$statistics
  namesstatistics <- stats$namesstatistics
  ratebeh <- stats$ratebeh
  netstats <- stats$netstats

  #startvalues (afster ts_statprep)
  startvalues <- ts_startprep(ans=ans, myeff=myeff, net1 = net1, net2 = net2,
                              ccovar = ccovar, startvalues= startvalues,
                              namesstatistics = namesstatistics)


  ### end of initialization



  # here we start with the actual phase 3. ###take start values from result of phase2


  crn <- runif(itef3, 1 - 2^31, 2^31 - 1) # common random numbers?
  pn <- length(startvalues)
  deviation <- rep(0.1, pn)
  deviation[names(startvalues) == "Rate"] <- startvalues[names(startvalues) == "Rate"] / 10
  deviation[names(startvalues) == "Rate beh"] <- startvalues[names(startvalues) == "Rate beh"] / 10

  if (returnDeps) resnets <-  vector(mode = "list", length = itef3) #list of matrix of simulated networks

  res_mat <- matrix(NA, nrow = pn, ncol = pn) # the jacobi matrix

  if (!parallel) {
    # to save the simulated z-scores
    res <-
      matrix(NA, nrow = itef3, ncol = pn + (pn * pn))


    for (i in 1:itef3) {
      if (verbose) print(i)

      # to save the simulated z-scores (for one itef3)
      rest <- rest2 <- rep(NA, pn)

      # Zbar
      set.seed(crn[i]) # is this just our crn
      sims1 <-
        ts_sim(
          net1 = net1,
          ccovar = ccovar,
          startvalues = startvalues,
          statistics = statistics,
          p2step = p2step,
          dist1 = dist1,
          dist2 = dist2,
          modet1 = modet1,
          modet2 = modet2
        )
      if (returnDeps) resnets[[i]] <- sims1$net_n

      # calculate the observed statistics
      if (!length(ratebeh)>0) {
        rest <- ts_targets(
          net1 = net1,
          net2 = sims1$net_n,
          ccovar = ccovar,
          statistics = statistics
        )
      } else {
        ccovar_n <- ccovar
        ccovar_n[,2] <- sims1$beh_n
        rest <- ts_targets(
          net1 = net1,
          net2 = sims1$net_n,
          ccovar = ccovar_n,
          statistics = statistics
        )
      }

      # start chaning the startvalues for rate and statistics
      for (j in 1:length(startvalues)) {
        #print(j)
        #set new startvalues
        startvalues_n <- startvalues
        startvalues_n[j] <- startvalues_n[j] + deviation[j]

        #simulate network with new startvalues
        #set the crn
        set.seed(crn[i])
        sims1 <-
          ts_sim(
            net1 = net1,
            ccovar = ccovar,
            startvalues = startvalues_n,
            statistics = statistics,
            p2step = p2step,
            dist1 = dist1,
            dist2 = dist2,
            modet1 = modet1,
            modet2 = modet2
          )


        # calculate the observed statistics
        if (!length(ratebeh)>0) {
          rest2 <- ts_targets(
            net1 = net1,
            net2 = sims1$net_n,
            ccovar = ccovar,
            statistics = statistics
          )
        } else {
          ccovar_n <- ccovar
          ccovar_n[,2] <- sims1$beh_n
          rest2 <- ts_targets(
            net1 = net1,
            net2 = sims1$net_n,
            ccovar = ccovar_n,
            statistics = statistics
          )
        }

        rest <- c(rest, rest2)
      }
      res[i,] <- rest
    }
    # now we have all z-scores we can calculate the jacobi
  }

  if (parallel & !returnDeps) {

   res <- foreach(i = 1: itef3, .combine="rbind", .verbose = verbose) %dopar% {

       # to save the simulated z-scores (for one itef1)
       rest <- rest2 <- rep(NA, pn)

       # Zbar
       set.seed(crn[i]) # is this just our crn
       sims1 <-
         ts_sim(
           net1 = net1,
           ccovar = ccovar,
           startvalues = startvalues,
           statistics = statistics,
           p2step = p2step,
           dist1 = dist1,
           dist2 = dist2,
           modet1 = modet1,
           modet2 = modet2
         )

       # calculate the observed statistics
       if (!length(ratebeh)>0) {
         rest <- ts_targets(
           net1 = net1,
           net2 = sims1$net_n,
           ccovar = ccovar,
           statistics = statistics
         )
       } else {
         ccovar_n <- ccovar
         ccovar_n[,2] <- sims1$beh_n
         rest <- ts_targets(
           net1 = net1,
           net2 = sims1$net_n,
           ccovar = ccovar_n,
           statistics = statistics
         )
       }

       # start chaning the startvalues for rate and statistics
       for (j in 1:length(startvalues)) {
         #print(j)
         #set new startvalues
         startvalues_n <- startvalues
         startvalues_n[j] <- startvalues_n[j] + deviation[j]

         #simulate network with new startvalues
         #set the crn
         set.seed(crn[i])
         sims1 <-
           ts_sim(
             net1 = net1,
             ccovar = ccovar,
             startvalues = startvalues_n,
             statistics = statistics,
             p2step = p2step,
             dist1 = dist1,
             dist2 = dist2,
             modet1 = modet1,
             modet2 = modet2
           )


         # calculate the observed statistics
         if (!length(ratebeh)>0) {
           rest2 <- ts_targets(
             net1 = net1,
             net2 = sims1$net_n,
             ccovar = ccovar,
             statistics = statistics
           )
         } else {
           ccovar_n <- ccovar
           ccovar_n[,2] <- sims1$beh_n
           rest2 <- ts_targets(
             net1 = net1,
             net2 = sims1$net_n,
             ccovar = ccovar_n,
             statistics = statistics
           )
         }
         rest <- c(rest, rest2)
       }
       rest
     }
  }

  if (parallel & returnDeps) {
    restemp <- foreach(i = 1: itef3) %dopar% {
      # to save the simulated z-scores (for one itef1)
      rest <- rest2 <- rep(NA, pn)

      # Zbar
      set.seed(crn[i]) # is this just our crn
      sims1 <-
        ts_sim(
          net1 = net1,
          ccovar = ccovar,
          startvalues = startvalues,
          statistics = statistics,
          p2step = p2step,
          dist1 = dist1,
          dist2 = dist2,
          modet1 = modet1,
          modet2 = modet2
        )

      # calculate the observed statistics
      if (!length(ratebeh)>0) {
        rest <- ts_targets(
          net1 = net1,
          net2 = sims1$net_n,
          ccovar = ccovar,
          statistics = statistics
        )
      } else {
        ccovar_n <- ccovar
        ccovar_n[,2] <- sims1$beh_n
        rest <- ts_targets(
          net1 = net1,
          net2 = sims1$net_n,
          ccovar = ccovar_n,
          statistics = statistics
        )
      }

      # start changing the startvalues for rate and statistics
      for (j in 1:length(startvalues)) {
        #print(j)
        #set new startvalues
        startvalues_n <- startvalues
        startvalues_n[j] <- startvalues_n[j] + deviation[j]

        #simulate network with new startvalues
        #set the crn
        set.seed(crn[i])
        sims1 <-
          ts_sim(
            net1 = net1,
            ccovar = ccovar,
            startvalues = startvalues_n,
            statistics = statistics,
            p2step = p2step,
            dist1 = dist1,
            dist2 = dist2,
            modet1 = modet1,
            modet2 = modet2
          )


        # calculate the observed statistics
        if (!length(ratebeh)>0) {
          rest2 <- ts_targets(
            net1 = net1,
            net2 = sims1$net_n,
            ccovar = ccovar,
            statistics = statistics
          )
        } else {
          ccovar_n <- ccovar
          ccovar_n[,2] <- sims1$beh_n
          rest2 <- ts_targets(
            net1 = net1,
            net2 = sims1$net_n,
            ccovar = ccovar_n,
            statistics = statistics
          )
        }
        rest <- c(rest, rest2)

      list(rest = rest, simnet = sims1$net_n)
      }
    }
    #convert nested list back to res and simnets
    resnets <- lapply(restemp, `[[`, 2)
    res <- do.call(rbind, lapply(restemp, `[[`, 1))
  }

  res_mat <- matrix(NA, nrow = pn, ncol = pn) # the jacobi matrix
  count <- pn + 1
  for (i in 1:pn) {
    for (j in 1:pn) {
      res_mat[i, j] <- mean((res[, count] - res[, j]) / deviation[i])
      count <- count + 1
    }
  }
  jacob <-
    t(res_mat)
  # perhaps return a list with more info, like the deviations and initial values
  # and stuff
  Sigma <- stats::cov(res[, 1:pn])
  covtheta <- solve(jacob) %*% Sigma %*% t(solve(jacob))
  sigmaj <- sqrt(diag(Sigma))
  zbar <- colMeans(res[, 1:pn])
  zobs <- ts_targets(net1 = net1,
                     net2 = net2,
                     ccovar = ccovar,
                     statistics = statistics)
  tstat <- (zbar - zobs)/sigmaj
  bstar <- (zbar - zobs)
  tconv.max <- sqrt(t(zbar - zobs) %*% solve(Sigma)  %*% (zbar - zobs))

  if (!returnDeps) {
    resultsphase3 <- list(estim = startvalues, covtheta=covtheta, tstat = tstat, tconv.max=tconv.max, zdevs = res, simnets = NULL)
  } else {
    resultsphase3 <- list(estim = startvalues, covtheta=covtheta, tstat = tstat, tconv.max=tconv.max, zdevs = res, simnets = resnets)
  }

  return(resultsphase3)
}

ts_dataprep <- function(ans = NULL, mydata = NULL, ccovar=NULL, ...) {

  #### be aware that it could be possible that we have no ccovars but do have depvars. need to correct for that.
  # check if there are ccovars stored in ans
  if (is.null(ccovar)) {
    if (!is.null(ans) & length(ans$f$Data1$cCovars) > 0) {
      data <-
        matrix(
          NA,
          nrow = length(ans$f$Data1$cCovars[[1]]),
          ncol = length(ans$f$Data1$cCovars)
        )
      for (i in 1:length(ans$f$Data1$cCovars)) {
        data[, i] <-
          as.numeric(ans$f$Data1$cCovars[[i]]) + attr(ans$f$Data1$cCovars[[i]], "mean")
      }
      ccovar <- as.data.frame(data)
      colnames(ccovar) <- names(ans$f$Data1$cCovars)
      if (length(ans$f$Data1$depvars) > 1 ) {
        dep <- data.frame(unclass(ans$f$Data1$depvars[[2]]))
        names(dep)[1] <-  paste0(attributes(ans$f$Data1$depvars)$names[2], ".1")
        names(dep)[2] <-  paste0(attributes(ans$f$Data1$depvars)$names[2], ".2")
        ccovar <- cbind(dep, ccovar)
      }
    } else if (length(mydata$cCovars) > 0) {
    data <-
      matrix(
        NA,
        nrow = length(mydata$cCovars[[1]]),
        ncol = length(mydata$cCovars)
      )
    for (i in 1:length(mydata$cCovars)) {
      data[, i] <-
        as.numeric(mydata$cCovars[[i]]) + attr(mydata$cCovars[[i]], "mean")
    }
    ccovar <- as.data.frame(data)
    colnames(ccovar) <- names(mydata$cCovars)
    if (length(mydata$depvars) > 1 ) {
      dep <- data.frame(mydata$depvars[[2]][])
      names(dep)[1] <-  paste0(attributes(mydata$depvars)$names[2], ".1")
      names(dep)[2] <-  paste0(attributes(mydata$depvars)$names[2], ".2")
      ccovar <- cbind(dep, ccovar)
    }
    }
  }
    # prepare dataset
  ccovar <- ts_prepdata(ccovar)

  return(ccovar)

}

ts_netprep <- function(ans = NULL, mydata = NULL, net1 = NULL, net2 = NULL, ...) {
  # starting networks
  if (is.null(net1)) {
    if (!is.null(ans)) {
      net1 <- ans$f$Data1$depvars[[1]][, , 1]
    } else if (!is.null(mydata)) {
      net1 <- mydata$depvars[[1]][, , 1]
    }
  }
  if (is.null(net2)) {
    if (!is.null(ans)) {
      net2 <- ans$f$Data1$depvars[[1]][, , 2]
    } else if (!is.null(mydata)) {
      net2 <- mydata$depvars[[1]][, , 2]
    }
  }
  return(list(net1 = net1, net2 = net2))
}

ts_statprep <- function(ans = NULL, myeff=NULL, statistics = NULL, ...) {
  # included statistics
  if (is.null(statistics)) {

    if (!is.null(ans)) {
      # namesstatistics <- ans$effects$shortName
      namesstatistics  <- ans$effects$shortName
      namesstatistics[namesstatistics == "density"] <- "degree"
      statistics <- ans$effects$shortName
      statistics[statistics == "density"] <- "degree"
      # remove "Rate" from the list
      statistics <- as.list((paste0("ts_", statistics))[-which(namesstatistics == "Rate")])
      names(statistics) <- namesstatistics[-which(namesstatistics == "Rate")]
      statistics <- lapply(statistics, get)
      # attach interaction1
      interaction1 <- ans$effects$interaction1[-which(namesstatistics == "Rate")]
      for (i in 1:length(interaction1)) {
        if (interaction1[i] != "") {
          statistics[[i]] <-
            list(statistics[[i]], interaction1[i])
        }
      }
      if (length(which(namesstatistics == "Rate") == 2)) namesstatistics[which(namesstatistics == "Rate")[2]] <- "Rate beh"


    } else if (!is.null(myeff)) {
      namesstatistics <- myeff$shortName[myeff$include]
      namesstatistics[namesstatistics == "density"] <- "degree"
      statistics <- myeff$shortName[myeff$include]
      statistics[statistics == "density"] <- "degree"
      # remove "Rate" from the list
      statistics <- as.list((paste0("ts_", statistics))[-which(namesstatistics == "Rate")])
      names(statistics) <- namesstatistics[-which(namesstatistics == "Rate")]
      statistics <- lapply(statistics, get)
      # attach interaction1
      interaction1 <- myeff$interaction1[myeff$include][-which(namesstatistics == "Rate")]
      for (i in 1:length(interaction1)) {
        if (interaction1[i] != "") {
          statistics[[i]] <-
            list(statistics[[i]], interaction1[i])
        }
      }
      if (length(which(namesstatistics == "Rate") == 2)) namesstatistics[which(namesstatistics == "Rate")[2]] <- "Rate beh"
    }
  } else {
    names(statistics) <- sapply(statistics, ts_names)
    namesstatistics <- c("Rate", sapply(statistics, ts_names))
    ratebeh <- which(namesstatistics == "linear")
    if (length(ratebeh) > 0) {
      namesstatistics <- c(namesstatistics[1:ratebeh - 1], "Rate beh", namesstatistics[ratebeh:length(namesstatistics)])
    }
  }


  ratebeh <- which(namesstatistics == "linear")
  netstats <- rep(TRUE, length(namesstatistics))
  if ("Rate beh" %in% namesstatistics) {
    netstats[which(namesstatistics == "Rate beh"): length(netstats)] <- FALSE
  }
  return(list(statistics = statistics, namesstatistics = namesstatistics, ratebeh=ratebeh, netstats = netstats))
}

ts_startprep <- function(ans = NULL, myeff=NULL, ccovar=NULL, net1 = NULL, net2 = NULL, startvalues = NULL, namesstatistics = NULL, ...) {


  if (is.null(startvalues)) {
    if (!is.null(ans)) {
      startvalues <- ans$theta
    } else if (!is.null(myeff)) {
      startvalues <- summary(myeff)$initialValue
    } else {
      # this is far from perfect. ASK TS how he gets initial values
      # startvalue_rate <- targets[1]/dim(net1)[1]
      # startvalue_param <- rep(0, length(statistics))
      # startvalue_param[1] <- qlogis(sum(net2)/(nrow(net2)*(nrow(net2) - 1)))
      # startvalues <- c(startvalue_rate, startvalue_param)
      mynet <-
        RSiena::sienaDependent(array(c(net1, net2), dim = c(dim(net1), 2)))
      if ("Rate beh" %in% namesstatistics) {
        dep_beh <- RSiena::sienaDependent(as.matrix(ccovar[,1:2]), type = "behavior")
        mydata <- RSiena::sienaDataCreate(mynet, dep_beh)
      } else {
        mydata <- RSiena::sienaDataCreate(mynet)
      }

      myeff <- RSiena::getEffects(mydata)

      startvalues <- rep(0, length(namesstatistics))
      startvalues[namesstatistics == "Rate"] <- summary(myeff)$initialValue[1]
      startvalues[namesstatistics == "degree"] <- summary(myeff)$initialValue[summary(myeff)$shortName == "density"]
      try(startvalues[namesstatistics == "Rate beh"] <- summary(myeff)$initialValue[summary(myeff)$shortName == "Rate" & summary(myeff)$name == "dep_beh"], silent = TRUE)
      try(startvalues[namesstatistics == "linear"] <- summary(myeff)$initialValue[summary(myeff)$shortName == "linear"], silent = TRUE)



    }
  }
  names(startvalues) <- namesstatistics
  return(startvalues)
}





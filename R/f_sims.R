#' Simulation of Network evolution
#'
#' @description `ts_sims` simulates the network evolution `nsims` times given
#' the existing network `net1`, the defined evaluation function [`ts_eval()`],
#' the included network statistics (taken from `ans`, `myeff`, or `statistics`)
#' and the corresponding parameter estimates (or starting values) taken from
#' `ans`, `myeff`, or `statistics`.
#'
#' @details For examples on how to use `ts_sims` see:
#' `vignette("1.Introduction_RsienaTwoStep", package="RsienaTwoStep")` or the
#' [package website](https://jochemtolsma.github.io/RsienaTwoStep/). Before you
#' set `parallel` to TRUE make sure to set-up a cluster with the package
#' `doParallel` (see `Examples`).
#' - `p2step==c(1,0,0)`: ministep
#' - `p2step==c(0,1,0)` & dist1==NULL & dist2==NULL: twostep-simultaneity
#' - `p2step==c(0,1,0)` & dist1!=NULL & dist2==NULL: twostep-strict coordination
#' - `p2step==c(0,1,0)` & dist1!=NULL & dist2!=NULL: twostep-weak coordination
#' - `p2step==c(0,0,1)`: two-ministeps
#' @inheritParams ts_estim
#' @param nsims numeric, number of simulations.
#' @param chain TRUE/FALSE, set to `TRUE` if you want to save all the subsequent
#'   networks (after the ministep or twostep) during the simulation. If `FALSE`
#'   only the end network is saved.
#' @param preparedata logical, if set to TRUE variables are demeaned and range
#'   is added. This is set to FALSE when running [`ts_estim()`] because here
#'   `ccovar` has already been prepared.
#' @return If `chain=FALSE` a `list` (of length `nsims`) of adjacency matrices
#'   representing the final network after the simulated evolution. If
#'   `chain=TRUE` a `list` of lists of adjacency matrices. Each inner list
#'   represents the complete network evolution of one simulation. The outer list
#'   refers to the simulation run (with length `nsims`).
#' @export
#' @seealso [`ts_estim()`], [`ts_alternatives_ministep()`],
#'   [`ts_alternatives_twostep()`], [`ts_alternatives_simstep()`], [`ts_eval()`]
#' @examples
#' ts_sims(
#'   net = ts_net2,
#'   nsims = 2,
#'   parallel = FALSE,
#'   statistics = list(ts_degree, ts_recip),
#'   startvalues = c(3, -2, 1),
#'   p2step = c(0, 1, 0)
#' )
#' #start with RSiena objects
#' \dontrun{
#' library(RSiena)
#' mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
#' mydata <- sienaDataCreate(mynet)
#' #toggle set conditional to retrieve the rate parameter in theta!
#' myalgorithm <- sienaAlgorithmCreate(cond=FALSE)
#' myeff <- getEffects(mydata)
#' ts_sims(
#'   mydata = mydata,
#'   myeff = myeff
#'   nsims = 2,
#'   parallel = FALSE,
#'   p2step = c(0, 1, 0))
#'  # or if you already used RSiena to estimate a model:
#'  ans1 <- siena07(myalgorithm, data=mydata, effects=myeff)
#'  ts_sims(
#'   ans = ans1
#'   nsims = 2,
#'   parallel = FALSE,
#'   p2step = c(0, 1, 0)
#' # Use a cluster
#' library(parallel)
#' n.cores <- parallel::detectCores() - 1  #save one core for other work
#' # create the cluster
#' my.cluster <- parallel::makeCluster(n.cores) #default PSOCK cluster
#' # register it to be used by %dopar%
#' doParallel::registerDoParallel(cl = my.cluster)
#' ts_sims(mydata = mydata,
#' myeff = myeff,
#' nsims = 20,
#' parallel = TRUE,
#' p2step = c(0, 1, 0))
#' stopCluster(my.cluster)
#' #Since the only official way to "unregister" a foreach backend is to register
#' the sequential backend: registerDoSEQ()
#'   }
#' @importFrom foreach %dopar%

ts_sims <- function(ans = NULL,
                    mydata = NULL,
                    myeff = NULL,
                    startvalues = NULL,
                    net1 = NULL,
                    ccovar = NULL,
                    preparedata = TRUE,
                    statistics = NULL,
                    nsims = 1000,
                    p2step = c(1, 0, 0),
                    dist1 = NULL,
                    dist2 = NULL,
                    modet1 = "degree",
                    modet2 = "degree",
                    chain = FALSE,
                    verbose = TRUE,
                    parallel = FALSE
                    ) {
  #### initialize function

  if (is.null(net1)) {
    if (!is.null(ans)) {
      net1 <- (ans$f$Data1$depvars$mynet)[, , 1]
    } else if (!is.null(mydata)) {
      net1 <- (mydata$depvars$mynet[, , 1])
    }
  }

  # included statistics
  if (is.null(statistics)) {
    if (!is.null(ans)) {
      # namesstatistics <- ans$effects$shortName
      statistics <- ans$effects$shortName
      statistics[statistics == "density"] <- "degree"
      statistics <- as.list((paste0("ts_", statistics))[-1])
      # names(statistics) <- ans$effects$shortName[-1]
      statistics <- lapply(statistics, get)
      for (i in 1:length(statistics)) {
        if (ans$effects$interaction1[i + 1] != "") {
          statistics[[i]] <-
            list(statistics[[i]], ans$effects$interaction1[i + 1])
        }
      }
    } else if (!is.null(myeff)) {
      # namesstatistics <- myeff$shortName[myeff$include]
      statistics <- myeff$shortName[myeff$include]
      statistics[statistics == "density"] <- "degree"
      statistics <- as.list((paste0("ts_", statistics))[-1])
      # names(statistics) <- myeff$shortName[myeff$include][-1]
      statistics <- lapply(statistics, get)
      for (i in 1:length(statistics)) {
        if (myeff$interaction1[myeff$include][i + 1] != "") {
          statistics[[i]] <-
            list(statistics[[i]], myeff$interaction1[myeff$include][i + 1])
        }
      }
    }
  }

  namesstatistics <- c("rate", sapply(statistics, ts_names))
  names(statistics) <- namesstatistics[-1]

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
          as.numeric(ans$f$Data1$cCovars[[i]]) +
          attr(ans$f$Data1$cCovars[[i]], "mean")
      }
      ccovar <- as.data.frame(data)
      colnames(ccovar) <- names(ans$f$Data1$cCovars)
    } else if (length(mydata$cCovars) > 0) {
      data <-
        matrix(
          NA,
          nrow = length(mydata$cCovars[[1]]),
          ncol = length(mydata$cCovars$cCovars)
        )
      for (i in 1:length(mydata$cCovars)) {
        data[, i] <-
          as.numeric(mydata$cCovars[[i]]) + attr(mydata$cCovars[[i]], "mean")
      }
      ccovar <- as.data.frame(data)
      colnames(ccovar) <- names(mydata$cCovars)
    }
  }


  if (is.null(startvalues)) {
    if (!is.null(ans)) {
      startvalues <- ans$theta
    } else if (!is.null(myeff)) {
      startvalues <- summary(myeff)$initialValue
    }
  }

  # prepare dataset
  if (!is.null(ccovar) & preparedata) {
    for (i in 1:ncol(ccovar)) {
      ccovar[, i] <- ts_centering(ccovar[, i])
      ccovar[, i] <- ts_simij(ccovar[, i])
    }
  }


  ### end of initialization



  if (parallel) {
    foreach::foreach(Nsim = 1:nsims) %dopar% {
      ts_sim(
        net1 = net1,
        ccovar = ccovar,
        startvalues = startvalues,
        statistics = statistics,
        p2step = p2step,
        chain = chain,
        dist1 = dist1,
        dist2 = dist2,
        modet1 = modet1,
        modet2 = modet2
      )
    }
  } else {
    sims <- list()
    for (i in 1:nsims) {
      if (verbose) print(paste0("nsim: ", i))
      sims[[i]] <- ts_sim(
        net1 = net1,
        ccovar = ccovar,
        startvalues = startvalues,
        statistics = statistics,
        p2step = p2step,
        chain = chain,
        dist1 = dist1,
        dist2 = dist2,
        modet1 = modet1,
        modet2 = modet2
      )
    }
    return(sims)
  }
}

#' @rdname ts_sims
#' @export
ts_sim <- function(
    ans = NULL,
    mydata = NULL,
    myeff = NULL,
    startvalues = NULL,
    net1 = NULL,
    ccovar = NULL,
    preparedata = FALSE,
    statistics = NULL,
    p2step = c(1, 0, 0),
    dist1 = NULL,
    dist2 = NULL,
    modet1 = "degree",
    modet2 = "degree",
    chain = FALSE,
    verbose = TRUE
  ) {

  ###I have experimented with doing the evaluations of the alternative networks in parallel.
  ###But the sapply function remains faster than parSapply or foreach.


  ###Initialize function

  if (is.null(net1)) {
    if (!is.null(ans)) {
      net1 <- (ans$f$Data1$depvars$mynet)[, , 1]
    } else if (!is.null(mydata)) {
      net1 <- (mydata$depvars$mynet[, , 1])
    }
  }

  # included statistics
  if (is.null(statistics)) {
    if (!is.null(ans)) {
      # namesstatistics <- ans$effects$shortName
      statistics <- ans$effects$shortName
      statistics[statistics == "density"] <- "degree"
      statistics <- as.list((paste0("ts_", statistics))[-1])
      # names(statistics) <- ans$effects$shortName[-1]
      statistics <- lapply(statistics, get)
      for (i in 1:length(statistics)) {
        if (ans$effects$interaction1[i + 1] != "") {
          statistics[[i]] <-
            list(statistics[[i]], ans$effects$interaction1[i + 1])
        }
      }
    } else if (!is.null(myeff)) {
      # namesstatistics <- myeff$shortName[myeff$include]
      statistics <- myeff$shortName[myeff$include]
      statistics[statistics == "density"] <- "degree"
      statistics <- as.list((paste0("ts_", statistics))[-1])
      # names(statistics) <- myeff$shortName[myeff$include][-1]
      statistics <- lapply(statistics, get)
      for (i in 1:length(statistics)) {
        if (myeff$interaction1[myeff$include][i + 1] != "") {
          statistics[[i]] <-
            list(statistics[[i]], myeff$interaction1[myeff$include][i + 1])
        }
      }
    }
  }

  namesstatistics <- c("rate", sapply(statistics, ts_names))
  names(statistics) <- namesstatistics[-1]

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
          as.numeric(ans$f$Data1$cCovars[[i]]) +
          attr(ans$f$Data1$cCovars[[i]], "mean")
      }
      ccovar <- as.data.frame(data)
      colnames(ccovar) <- names(ans$f$Data1$cCovars)
    } else if (length(mydata$cCovars) > 0) {
      data <-
        matrix(
          NA,
          nrow = length(mydata$cCovars[[1]]),
          ncol = length(mydata$cCovars$cCovars)
        )
      for (i in 1:length(mydata$cCovars)) {
        data[, i] <-
          as.numeric(mydata$cCovars[[i]]) + attr(mydata$cCovars[[i]], "mean")
      }
      ccovar <- as.data.frame(data)
      colnames(ccovar) <- names(mydata$cCovars)
    }
  }


  if (is.null(startvalues)) {
    if (!is.null(ans)) {
      startvalues <- ans$theta
    } else if (!is.null(myeff)) {
      startvalues <- summary(myeff)$initialValue
    }
  }

  # prepare dataset
  if (!is.null(ccovar) & preparedata) {
    for (i in 1:ncol(ccovar)) {
      ccovar[, i] <- ts_centering(ccovar[, i])
      ccovar[, i] <- ts_simij(ccovar[, i])
    }
  }


  ###End of Initialization

  rate <- startvalues[1]
  parameters <- startvalues[-1]
  nministep <- rate * nrow(net1) + 1
  net_n <- net1
  nets <- list()
  ministep <- 1
  iteration <- 1
  while (ministep < nministep) {
    # normal ministep or 2step?
    normal <- sample(c(1, 2, 3), 1, prob = p2step)

    # normal model
    if (normal == 1) {
      # sample agent
      ego <- ts_select(net1)
      # options
      options <- ts_alternatives_ministep(net = net_n, ego = ego)
      # evaluations
      eval <-
        sapply(
          options,
          FUN = ts_eval,
          ccovar = ccovar,
          ego = ego,
          statistics = statistics,
          parameters = parameters
        )
      eval <- eval - max(eval)
      # pick new network
      net_n <-
        options[[sample(1:length(eval),
          size = 1,
          prob = exp(eval) / sum(exp(eval))
        )]]
      if (chain) {
        nets[[iteration]] <- net_n
      }
      iteration <- iteration + 1
      ministep <- ministep + 1
    }

    # model with simultaneity
    if (normal == 2) {
      results <-
        ts_alternatives_twostep(
          net = net_n,
          dist1 = dist1,
          dist2 = dist2,
          modet1 = modet1,
          modet2 = modet2
        )
      egos <- results[[1]] # sampled dyad
      options <-
        results[[2]] # all possible future networks after the twostep
      if (is.null(egos) & is.null(options)) {
        iteration <- iteration + 1
      }

      if (!is.null(egos) &
        !is.null(options)) {
        # check if it was possible to sample agents
        # evaluations ego1
        eval1 <-
          sapply(
            options,
            FUN = ts_eval,
            ccovar = ccovar,
            ego = egos[1],
            statistics = statistics,
            parameters = parameters
          )
        # evaluations ego2
        eval2 <-
          sapply(
            options,
            FUN = ts_eval,
            ccovar = ccovar,
            ego = egos[2],
            statistics = statistics,
            parameters = parameters
          )
        # pick new network
        eval <- eval1 + eval2
        eval <- eval - max(eval)
        net_n <-
          options[[sample(1:length(eval),
            size = 1,
            prob = exp(eval) / sum(exp(eval))
          )]]
        if (chain) {
          nets[[iteration]] <- net_n
        }
        iteration <- iteration + 1
        ministep <- ministep + 2
      }
    }

    # model with two simultaneous ministeps of the same ego
    if (normal == 3) {
      # sample agent
      ego <- ts_select(net1)
      # options
      options <- ts_alternatives_simstep(net = net_n, ego = ego)
      # evaluations
      eval <-
        sapply(
          options,
          FUN = ts_eval,
          ccovar = ccovar,
          ego = ego,
          statistics = statistics,
          parameters = parameters
        )
      eval <- eval - max(eval)
      # pick new network
      net_n <-
        options[[sample(1:length(eval),
          size = 1,
          prob = exp(eval) / sum(exp(eval))
        )]]
      if (chain) {
        nets[[iteration]] <- net_n
      }
      iteration <- iteration + 1
      ministep <- ministep + 2
    }
  }
  if (chain) {
    return(nets)
  } else {
    return(net_n)
  }
}

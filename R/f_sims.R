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
#' @importFrom foreach %dopar%

ts_sims <- function(ans = NULL,
                    mydata = NULL,
                    myeff = NULL,
                    startvalues = NULL,
                    net1 = NULL,
                    ccovar = NULL,
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

  # starting networks
  nets <- ts_netprep(ans=ans, mydata=mydata, net1=net1)
  net1 <- nets$net1

  # included statistics
  stats <- ts_statprep(ans = ans, myeff = myeff, statistics = statistics)
  statistics <- stats$statistics
  namesstatistics <- stats$namesstatistics
  ratebeh <- stats$ratebeh
  netstats <- stats$netstats

  #startvalues (afster ts_statprep)
  startvalues <- ts_startprep(ans=ans, myeff=myeff, startvalues= startvalues, namesstatistics = namesstatistics)

  # ccovar (and deps)
  ccovar <- ts_dataprep(ans= ans, mydata= mydata, ccovar=ccovar)


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
  ###Need to tweak if behavioral only model.

  ###Initialize function

  # starting networks
  nets <- ts_netprep(ans=ans, mydata=mydata, net1=net1)
  net1 <- nets$net1

  # included statistics
  stats <- ts_statprep(ans = ans, myeff = myeff, statistics = statistics)
  statistics <- stats$statistics
  namesstatistics <- stats$namesstatistics
  ratebeh <- stats$ratebeh
  netstats <- stats$netstats

  #startvalues (afster ts_statprep)
  startvalues <- ts_startprep(ans=ans, myeff=myeff, startvalues= startvalues, namesstatistics = namesstatistics)

  # ccovar (and deps)
  ccovar <- ts_dataprep(ans= ans, mydata= mydata, ccovar=ccovar)



  ###End of Initialization

  rate <- startvalues[1]
  parameters <- startvalues[-1]
  nministep <- rate * nrow(net1) + 1
  nministepb <- Inf
  net_n <- net1
  ccovar_n <- ccovar
  nets <- list()
  #if only network part
  parametersn <- parameters
  statisticsn <- statistics

  # if we have a behavioral part
  if (length(ratebeh) > 0) {
    rateb <- startvalues[namesstatistics == "Rate beh"]
    nministepb <- rateb * nrow(net1) + 1
    beh_n <- ccovar[,1]
    behs <- list()
    #split statistics and parameters in network and behavior part
    parametersn <- startvalues[2:(which(namesstatistics == "Rate beh") -1)]
    parametersb <- startvalues[(which(namesstatistics == "Rate beh") +1):length(namesstatistics)]
    statisticsn <- statistics[1:(which(names(statistics) == "linear")-1)]
    statisticsb <- statistics[which(names(statistics) == "linear"):length(names(statistics))]
  }

  #it would be nice if I could preallocate the nets list, but since I can mix p2steps with ministeps this is not easy.
  # if (chain & normal == 1) nets <- vector("list", length = nministep -1 )
  # if (chain & normal != 1) nets <- vector("list", length = nministep/2 )
  ministep <- 1
  ministepb <- 1
  iteration <- 1
  typestep <- 1
  typesteps <- 0
  while ((ministep < nministep) & (ministepb < nministepb)) {

    # normal ministep or 2step?
    normal <- sample(c(1, 2, 3), 1, prob = p2step) #improve! we do not need to sample if no mixing.

    # net change or beh change
    if (length(ratebeh) > 0) {
    typestep <- sample(c(1,2), 1, prob = c(rate, rateb))
    typesteps <- c(typesteps, typestep)
    }

    #network change
    if (typestep == 1) {
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
          ccovar = ccovar_n,
          ego = ego,
          statistics = statisticsn,
          parameters = parametersn
        )

      # #eval <- eval - max(eval)
      # # pick new network
      # net_n <-
      #   options[[sample(1:length(eval),
      #     size = 1,
      #     prob = exp(eval) / sum(exp(eval))
      #   )]]


      # # pick new network
      net_n <-
        options[[sample(1:length(eval),
          size = 1,
          prob = exp(eval)
        )]]

      if (chain) {
        nets[[iteration]] <- net_n
      }
      iteration <- iteration + 1
      ministep <- ministep + 1
    } else if (normal == 2) {  # twostep
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
        if (chain) {
          nets[[iteration]] <- net_n
        }
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
            ccovar = ccovar_n,
            ego = egos[1],
            statistics = statisticsn,
            parameters = parametersn
          )
        # evaluations ego2
        eval2 <-
          sapply(
            options,
            FUN = ts_eval,
            ccovar = ccovar_n,
            ego = egos[2],
            statistics = statisticsn,
            parameters = parametersn
          )
        # pick new network
        eval <- eval1 + eval2


        #McFadden choice function.
        eval <- eval - max(eval)
        net_n <-
          options[[sample(1:length(eval),
            size = 1,
            prob = exp(eval)
          )]]

        if (chain) {
          nets[[iteration]] <- net_n
        }
        iteration <- iteration + 1
        ministep <- ministep + 2

      }
    } else if (normal == 3) { #simstep
      # sample agent
      ego <- ts_select(net1)
      # options
      options <- ts_alternatives_simstep(net = net_n, ego = ego)
      # evaluations
      eval <-
        sapply(
          options,
          FUN = ts_eval,
          ccovar = ccovar_n,
          ego = ego,
          statistics = statisticsn,
          parameters = parametersn
        )

      # pick new network
      eval <- eval - max(eval)
      net_n <-
        options[[sample(1:length(eval),
          size = 1,
          prob = exp(eval)
        )]]

      if (chain) {
        nets[[iteration]] <- net_n
      }
      iteration <- iteration + 1
      ministep <- ministep + 2
    }
    }

    #behavior change
    if (typestep == 2) {
      if (normal == 1) {
        # sample agent
        ego <- ts_select(net1)
        # options
        options <- ts_alternatives_ministep_beh(beh = beh_n, ego = ego)
        # evaluations
        eval <-
          sapply(
            options,
            FUN = ts_eval_beh,
            net = net_n,
            ccovar = ccovar_n,
            ego = ego,
            statistics = statisticsb,
            parameters = parametersb
          )

        # #eval <- eval - max(eval)
        # # pick new network
        # net_n <-
        #   options[[sample(1:length(eval),
        #     size = 1,
        #     prob = exp(eval) / sum(exp(eval))
        #   )]]


        # # pick new behavior
        beh_n <-
          options[[sample(1:length(eval),
                          size = 1,
                          prob = exp(eval)
          )]]
        ccovar_n[,1] <- beh_n
        if (chain) {
          behs[[iteration]] <- beh_n
        }
        iteration <- iteration + 1
        ministepb <- ministepb + 1
      } else if (normal == 2) {  # twostep
        results <-  ts_alternatives_twostep_beh(beh = beh_n, net=net_n, dist1 = dist1, modet1 = modet1)
        egos <- results[[1]] # sampled dyad
        options <-
          results[[2]] # all possible future networks after the twostep
        if (is.null(egos) & is.null(options)) {
          if (chain) {
            behs[[iteration]] <- beh_n
          }
          iteration <- iteration + 1
        }

        if (!is.null(egos) &
            !is.null(options)) {
          # check if it was possible to sample agents
          # evaluations ego1
          eval1 <-
            sapply(
              options,
              FUN = ts_eval_beh,
              net = net_n,
              ccovar = ccovar_n,
              ego = egos[1],
              statistics = statisticsb,
              parameters = parametersb
            )
          # evaluations ego2
          eval2 <-
            sapply(
              options,
              FUN = ts_eval_beh,
              net = net_n,
              ccovar = ccovar_n,
              ego = egos[2],
              statistics = statisticsb,
              parameters = parametersb
            )
          # pick new behavior
          eval <- eval1 + eval2


          #McFadden choice function.
          eval <- eval - max(eval)
          beh_n <-
            options[[sample(1:length(eval),
                            size = 1,
                            prob = exp(eval)
            )]]
          ccovar_n[,1] <- beh_n
          if (chain) {
            behs[[iteration]] <- beh_n
          }
          iteration <- iteration + 1
          ministepb <- ministepb + 2

        }
      } else if (normal == 3) {
        # model with two simultaneous ministeps of the same ego
        # sample agent
        ego <- ts_select(net1)
        # options
        options <- ts_alternatives_simstep_beh(beh = beh_n, ego = ego)
        # evaluations
        eval <-
          sapply(
            options,
            FUN = ts_eval_beh,
            net = net_n,
            ccovar = ccovar_n,
            ego = ego,
            statistics = statisticsb,
            parameters = parametersb
          )

        # pick new network
        eval <- eval - max(eval)
        beh_n <-
          options[[sample(1:length(eval),
                          size = 1,
                          prob = exp(eval)
          )]]
        ccovar_n[,1] <- beh_n
        if (chain) {
          behs[[iteration]] <- beh_n
        }
        iteration <- iteration + 1
        ministepb <- ministepb + 2
      }
    }

  }



  # return depends on chain and ratebeh
  if (!length(ratebeh) > 0) {
    if (chain) {
      return(list(final = list(net_n = net_n), chain = list(nets = nets)))
    } else {
      return(net_n)
    }
  } else {
    if (chain) {
      typesteps <- typesteps[-1]
      return(list(final = list(net_n = net_n, beh_n = beh_n), chain = list(nets = nets, behs = behs, typesteps = typesteps)))
    } else {
      return(list(net_n = net_n, beh_n = beh_n))
    }
  }

}

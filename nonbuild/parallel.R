#' If you have set-up a cluster and want to estimate a model with [`ts_estim()`] and set `parallel == TRUE`, in `ts_sims()`, `parallel2` is set to `TRUE`. This means the evaluation scores of the alternative networks are calculated in [`ts_eval()`] in parallel.
#' If you use `ts_sims()` to simulate networks, you may want to set either `parallel1` or `parallel2` to `TRUE`. Setting both to true will (probably) not work. `parallel1 == TRUE` means the possible outcome networks are simulated in parallel.
library(doParallel)
no_cores <- detectCores() - 1
mycl <- makeCluster(rep("localhost", no_cores))
clusterEvalQ(mycl, library(RsienaTwoStep)) #I think this is crucial
#clusterExport(mycl, varlist=c("ts_sim")) #I think this is crucial
registerDoParallel(mycl)
#stopCluster(mycl)


# via RSiena

library(RSiena)
mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet)
myalgorithm <- sienaAlgorithmCreate(seed=1293, cond=FALSE, findiff = TRUE) #toggle set conditional to retrieve the rate parameter in theta!
myeff <- getEffects(mydata)
#myeff <- includeEffects(myeff, transTrip, inAct)
ans1 <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE, returnDeps = TRUE)
ans1

system.time(jac <- ts_phase1(net1 = s501,net2 = s502, statistics = list(ts_degree, ts_recip),
                 ccovar = NULL,
                  itef1 = 30, p2step = c(1,0,0))
)
# user  system elapsed
# 23.70    6.08   30.03

system.time(jac <- ts_phase1(net1 = s501,net2 = s502, statistics = list(ts_degree, ts_recip),
                             ccovar = NULL,
                             itef1 = 30,
                             parallel = TRUE)
)
# user  system elapsed
# 0.10    0.05   24.64

system.time(ans <- ts_estim(net1 = s501, net2 = s502,
                 statistics = list(ts_degree, ts_recip),
                 p2step = c(1,0,0),
                 nite = 30,
         itef1 = 30,
         itef3 = 50,
         phase1 = TRUE,
         phase3 = TRUE,
         parallel = TRUE,
         verbose = TRUE)
)

system.time(ans2 <- ts_estim(ans = ans1,
                            p2step = c(0,1,0),
                            nite = 30,
                            itef1 = 10,
                            itef3 = 30,
                            dist1 = 2,
                            phase1 = TRUE,
                            phase3 = TRUE,
                            parallel = FALSE,
                            verbose = TRUE)
)

ts_estim(ans = ans1,
         p2step = c(0,1,0),
         nite = 30,
         itef1 = 10,
         itef3 = 10,
         phase1 = TRUE,
         phase3 = TRUE,
         parallel = TRUE,
         verbose = TRUE)
)


# user  system elapsed
# 1224.84  130.03 1383.06
# > ans2
# estim        SE
# rate    5.426663 1.1825682
# degree -2.138892 0.2599478
# recip   2.530833 0.3751216


ans
# user  system elapsed
# 70.36   15.64   86.45

system.time(ans <- ts_estim(net1 = s501, net2 = s502,
                            statistics = list(ts_degree, ts_recip),
                            p2step = c(1,0,0),
                            nite = 100,
                            itef1 = 30,
                            itef3 = 200,
                            phase1 = TRUE,
                            phase3 = TRUE,
                            parallel = TRUE,
                            verbose = TRUE)
)
ans
# user  system elapsed
# 12.38    2.45  202.09
# > ans
# estim        SE
# rate    5.255209 0.7208893
# degree -2.091314 0.1117334
# recip   2.157001 0.2106615

Estimate   Standard   Convergence
Error      t-ratio
1. rate basic rate parameter mynet  5.5096  ( 0.8475   )   -0.0810
2. eval outdegree (density)        -2.2331  ( 0.1284   )   -0.0473
3. eval reciprocity                 2.4378  ( 0.2420   )    0.0092

Overall maximum convergence ratio:    0.1093



ts_estim(mydata = mydata, myeff = myeff, phase1 = TRUE, phase3 = FALSE, verbose = TRUE, p2step=c(1,0,0))

system.time(ans <- ts_estim(ans=ans1,
                            p2step = c(1,0,0),
                            nite = 100,
                            itef1 = 10,
                            itef3 = 10,
                            phase1 = TRUE,
                            phase3 = TRUE,
                            parallel = FALSE,
                            verbose = TRUE)
)
# user  system elapsed
# 88.52    3.21   92.47

system.time(ans <- ts_estim(ans=ans1,
                            p2step = c(1,0,0),
                            nite = 100,
                            itef1 = 10,
                            itef3 = 10,
                            phase1 = TRUE,
                            phase3 = TRUE,
                            parallel = TRUE,
                            verbose = TRUE)
)

                 #' ts_estim(net1 = s501,
                 #' net2 = s502,
                 #' statistics = list(ts_degree, ts_recip),
                 #' nite = 30,
                 #' itef1 = 10,
                 #' phase1 = TRUE)



ego<- 3
statistics = list(ts_degree, ts_recip, ts_cycle3, ts_inAct, ts_inPop)
parameters = c(-2,2)
ccovar = NULL

options <- ts_alternatives_simstep(net = s501, ego = ego)

options <- ts_alternatives_twostep(net = s501)
options <- options[[2]]
str(options, 1)
# evaluations

system.time(
eval <-
  sapply(
    options,
    FUN = ts_eval,
    ccovar = ccovar,
    ego = ego,
    statistics = statistics,
    parameters = parameters
  )
)

system.time(
  eval <-
    parSapply(cl = mycl,
      options,
      FUN = ts_eval,
      ccovar = ccovar,
      ego = ego,
      statistics = statistics,
      parameters = parameters
    )
)


system.time(
  foreach(i = 1:length(options), .combine="c") %dopar% {
  ts_eval(options[[i]],
          ccovar = ccovar,
          ego = ego,
          statistics = statistics,
          parameters = parameters)
}
)




ts_phase1(net1=s501, net2=s502, statistics = list(ts_degree, ts_recip), itef1=100, parallel = TRUE )


itef1 <- 100
crn <- sample(12345:4567890, itef1)
net1 <- s501
startvalue_rate <- 5
statistics <- list(ts_degree, ts_recip)
startvalue_param <- c(-2,2)
p2step <- c(1,0,0)
deviation <- rep(0.1, length(statistics) + 1)
deviation[1] <- startvalue_rate / 10

res <- foreach(i = 1: itef1, .combine="rbind") %dopar% {
    set.seed(crn[i])

    # Zbar
    set.seed(crn[i])
    sim_net <-
      ts_sim(
        net = net1,
        rate = startvalue_rate,
        statistics = statistics,
        parameters = startvalue_param,
        p2step = p2step,
        dist1 = dist1,
        dist2 = dist2,
        modet1 = modet1,
        modet2 = modet2
      )
    rest <-
      ts_targets(net1 = net1,
                 net2 = sim_net,
                 statistics = statistics)

    # rate
    set.seed(crn[i])
    sim_net <-
      ts_sim(
        net = net1,
        rate = startvalue_rate + deviation[1],
        statistics = statistics,
        parameters = startvalue_param,
        p2step = p2step,
        dist1 = dist1,
        dist2 = dist2,
        modet1 = modet1,
        modet2 = modet2
      )
    rest2 <-
      ts_targets(net1 = net1,
                 net2 = sim_net,
                 statistics = statistics)
    rest <- c(rest, rest2)

    # statistics
    rest3 <-
    foreach(j = 1:length(statistics), .combine="cbind") %do% {
      set.seed(crn[i])
      startvalue_param2 <- startvalue_param
      startvalue_param2[j] <-
        startvalue_param2[j] + deviation[1 + j]
      sim_net <-
        ts_sim(
          net = net1,
          rate = startvalue_rate,
          statistics = statistics,
          parameters = startvalue_param2,
          p2step = p2step,
          dist1 = dist1,
          dist2 = dist2,
          modet1 = modet1,
          modet2 = modet2
        )
      rest3 <-
        ts_targets(net1 = net1,
                   net2 = sim_net,
                   statistics = statistics)
    }
    rest <- c(rest, rest3)
    rest
    }
  }

test <- NA
for (i in 1:nrow(s501)) {
  test[i] <- ts_transTrip(s501, i)
}
test
library(microbenchmark)
library(RsienaTwoStep)
ts_eval2(net= s502, ego=10, statistics=list(ts_degree, ts_recip, ts_transTrip), ccovar = NULL, parameters= c(-3,3, 4))
microbenchmark(times = 1000, unit = "ms", ts_eval(net= s502, ego=10, statistics=list(ts_degree, ts_recip, ts_transTrip), ccovar = NULL, parameters= c(-3,3, 4)))
microbenchmark(times = 1000, unit = "ms", ts_eval(net=s502, ego=2, ccovar=df_ccovar1, statistics=list(ts_degree, ts_recip, ts_transTrip, ts_transMedTrip, list(ts_egoX, "cov1")), parameters=c(-2,2,7,7,1)))
microbenchmark(times = 1000, unit = "ms", ts_eval2(net=s502, ego=2, ccovar=df_ccovar1, statistics=list(ts_degree, ts_recip, ts_transTrip, ts_transMedTrip, list(ts_egoX, "cov1")), parameters=c(-2,2,7,7,1)))



microbenchmark(times = 1000, unit = "ms",
        ts_eval(net=ts_net1, ego=10, ccovar=df_ccovar1, statistics=list(ts_degree, ts_recip, ts_transTrip,
        ts_transMedTrip, list(ts_egoX, "cov1")), parameters=c(-2,2,7,7,1))
)
library(microbenchmark)
microbenchmark(times = 1000, unit = "ms",
               ts_eval2(net=ts_net1, ego=10, ccovar=df_ccovar1, statistics=list(ts_degree, ts_recip, ts_transTrip,
                                                                               ts_transMedTrip, list(ts_egoX, "cov1")), parameters=c(-2,2,7,7,1))
)



ts_eval(net= s502, ego=5, statistics=list(ts_degree, ts_recip, ts_transTrip))

microbenchmark(times = 1000, unit = "ms", ts_transTrip2(s501, 5))
microbenchmark(times = 1000, unit = "ms", ts_3(s501, 5))
library(compiler)
ts_transTrip2 <- cmpfun(ts_transTrip)

ts_3 <- function (net, ego) {
  statistic <- 0
  alters <- which(net[ego, ] == 1)
  if (length(alters) > 1) {
    statistic <- statistic + sum(net[alters, alters])
  }
  return(statistic)
}



ts_net1[1:10, 1:10]

ts_transTrip

ts_targets(net1 = s501, net2 = s502, statistics=list(ts_degree, ts_recip, ts_transTrip,
                                                                  ts_transMedTrip))
j=4
for (j in seq_along(statistics)) {
  stat <- statistics[[j]]

  if (length(stat) == 1) {
    # Single-argument statistic function
    Z[j + 1] <- sum(stat(net2, 1:nrow(net2)))
  } else if (length(stat) == 2) {
    # Two-argument statistic function
    Z[j + 1] <- sum(stat[[1]](net2, 1:nrow(net2), ccovar[, stat[[2]]]))
  }
}

# Calculate the observed statistics
Z <- numeric(length(statistics) + 1)
Z[1] <- sum(abs(net2 - net1))  # Rate

for (j in seq_along(statistics)) {
  stat <- statistics[[j]]

  if (length(stat) == 1) {
    # Single-argument statistic function
    Z[j + 1] <- sum(sapply(1:nrow(net2), function(i) stat(net2, i)))
  #} else if (length(stat) == 2) {
  } else {
     # Two-argument statistic function
    Z[j + 1] <- sum(sapply(1:nrow(net2), function(i) stat[[1]](net2, i, ccovar[, stat[[2]]])))
  }
}

foreach::foreach(i = 1:nrow(net2), .combine = "c") %dopar%
         statistics[[j]](net = net2, ego = i) |>
         sum()


# Check if a function is precompiled with cmpfun
is_precompiled <- function(fun) {
  return("bytecode" %in% attributes(fun))
}

# Example usage:
my_function <- function(x) {
  x + 1
}

# Precompile the function
my_function <- compiler::cmpfun(my_function)

# Check if the function is precompiled
is_precompiled(ts_degree)  # Should return TRUE

library(utils)
library(Rprof)

pool <- which(ts_net1==1, arr.ind = TRUE)
dyad <- pool[sample(nrow(pool), 1),]

ts_select(net=ts_net1, steps=2, dist1=2, modet1="degree")
res <- ts_alternatives_twostep(net=ts_net1, dist1=1, modet1="degree")
res[[1]]
ts_net1


install.packages("doSNOW")
library(doSNOW)
library(tcltk)
nw <- 4  # number of workers
cl <- makeSOCKcluster(nw)
registerDoSNOW(cl)

x <- iris[which(iris[,5] != 'setosa'), c(1,5)]
niter <- 15e+4
chunksize <- 4000  # may require tuning for your machine
maxcomb <- nw + 1  # this count includes fobj argument
totaltasks <- ceiling(niter / chunksize)

comb <- function(fobj, ...) {
  for(r in list(...))
    writeBin(r, fobj)
  fobj
}

final <- function(fobj) {
  close(fobj)
  t(matrix(readBin('temp.bin', what='double', n=niter*2), nrow=2))
}

mkprogress <- function(total) {
  pb <- tkProgressBar(max=total,
                      label=sprintf('total tasks: %d', total))
  function(n, tag) {
    setTkProgressBar(pb, n,
                     label=sprintf('last completed task: %d of %d', tag, total))
  }
}
opts <- list(progress=mkprogress(totaltasks))
resultFile <- file('temp.bin', open='wb')

r <-
  foreach(n=idiv(niter, chunkSize=chunksize), .combine='comb',
          .maxcombine=maxcomb, .init=resultFile, .final=final,
          .inorder=FALSE, .options.snow=opts) %dopar% {
            do.call('c', lapply(seq_len(n), function(i) {
              ind <- sample(100, 100, replace=TRUE)
              result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
              coefficients(result1)
            }))
          }

r <-
  foreach(n=idiv(niter, chunkSize=chunksize), .combine='comb',
          .maxcombine=maxcomb,
          .inorder=FALSE) %dopar% {
            do.call('c', lapply(seq_len(n), function(i) {
              ind <- sample(100, 100, replace=TRUE)
              result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
              coefficients(result1)
            }))
          }


# I included a progress bar since this example takes several hours to execute.
#
# Note that this example also uses the idiv function from the iterators package to increase the amount of work in each of the tasks. This technique is called chunking, and often improves the parallel performance. However, using idiv messes up the task indices, since the variable i is now a per-task index rather than a global index. For a global index, you can write a custom iterator that wraps idiv:

  idivix <- function(n, chunkSize) {
    i <- 1
    it <- idiv(n, chunkSize=chunkSize)
    nextEl <- function() {
      m <- nextElem(it)  # may throw 'StopIterator'
      value <- list(i=i, m=m)
      i <<- i + m
      value
    }
    obj <- list(nextElem=nextEl)
    class(obj) <- c('abstractiter', 'iter')
    obj
  }

# The values emitted by this iterator are lists, each containing a starting index and a count. Here's a simple foreach loop that uses this custom iterator:

r <-
  foreach(a=idivix(10, chunkSize=3), .combine='c') %dopar% {
    do.call('c', lapply(seq(a$i, length.out=a$m), function(i) {
      i
    }))
  }

require("doParallel")
nw <- 8
registerDoParallel(nw)
x <- iris[which(iris[,5] != "setosa"), c(1,5)]
niter <- 10000000
eval <- 1:10000000
r <- foreach(n=idiv(niter, chunks=nw), .combine='rbind') %dopar% {
  do.call('rbind', lapply(seq_len(n), function(i) {
    ind <- sample(100, 100, replace=TRUE)
    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
    coefficients(result1)
  }))
}

r <- foreach(n=idiv(niter, chunks=nw), .combine='rbind', .inorder = FALSE) %dopar% {
  do.call('rbind', lapply(seq_len(n), function(i) {

    i^2

  }))
}
stopCluster()

#trying to do the evaluations in chunks
no_cores <- 2
no_cores <- detectCores() - 1
mycl <- makeCluster(rep("localhost", no_cores))
clusterEvalQ(mycl, library(RsienaTwoStep))
clusterEvalQ(mycl, library("network"))
clusterEvalQ(mycl, library("RSiena"))
clusterEvalQ(mycl, library("sna"))
registerDoParallel(mycl)
stopCluster(cl = mycl)

library(future.apply)

## The below is same as plan(multisession, workers=4)
plan(cluster, workers=mycl)

xs <- 1:100
results <- future_lapply(xs, FUN=function(x) {
  Sys.sleep(0.1)
  sqrt(x)
})

STATS <- list(ts_degree, ts_recip, ts_transTrip, ts_cycle3, list(ts_simX, "smoke"), list(ts_altX, "alcohol"), list(ts_egoX, "alcohol"), list(ts_egoXaltX, "alcohol"))
ccovar <- data.frame(alcohol = s50a[, 1], smoke = s50s[, 1])
preparedata <- TRUE
# prepare dataset
if (!is.null(ccovar) & preparedata) {
  for (i in 1:ncol(ccovar)) {
    ccovar[, i] <- ts_centering(ccovar[, i])
    ccovar[, i] <- ts_simij(ccovar[, i])
  }
}
STARTS <- ans3$theta

results <-
  ts_alternatives_twostep(
    net = s501,
    dist1 = 2,
    dist2 = NULL
  )
egos <- results[[1]] # sampled dyad
options <-
  results[[2]] # all possible future networks after the twostep

niter <- length(options)
nw <- detectCores() - 1

library(microbenchmark)
microbenchmark(times = 10, unit = "ms",
r <- foreach(n=idiv(niter, chunks=4), .combine='c', .inorder = FALSE) %dopar% {
  do.call('c', lapply(seq_len(n), function(i) {
    eval1 <- ts_eval(options[[i]], ccovar = ccovar,
                     ego = egos[1],
                     statistics = STATS,
                     parameters = STARTS )
    eval2 <- ts_eval(options[[i]], ccovar = ccovar,
                     ego = egos[1],
                     statistics = STATS,
                     parameters = STARTS )
    eval <- eval1 + eval2
    eval
  }))
}
)

library(microbenchmark)
microbenchmark(times = 10, unit = "ms",
{
eval1 <- future_sapply(options, future.chunk.size=Inf, FUN=function(option) {ts_eval(option,ccovar = ccovar,
                                              ego = egos[1],
                                              statistics = STATS,
                                              parameters = STARTS)})
eval2 <- future_sapply(options, future.chunk.size=Inf, FUN=function(option) {ts_eval(option,ccovar = ccovar,
                                                              ego = egos[2],
                                                              statistics = STATS,
                                                              parameters = STARTS)})
eval <- eval1 + eval2
}
)

library(microbenchmark)
microbenchmark(times = 10, unit = "ms",
{
 eval1 <-
    sapply(
      options,
      FUN = ts_eval,
      ccovar = ccovar,
      ego = egos[1],
      statistics = STATS,
      parameters = STARTS
    )
  eval2 <-
    sapply(
      options,
      FUN = ts_eval,
      ccovar = ccovar,
      ego = egos[2],
      statistics = STATS,
      parameters = STARTS
    )
  eval <- eval1 + eval2
}
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

  https://datascienceplus.com/strategies-to-speedup-r-code/




    library(pbapply)
  library(doParallel)
  registerDoParallel(cl <- makeCluster(3))

  my_fun <- function(repetitions = 10){
    foreach( i = 1:repetitions, .combine = 'rbind' ) %dopar% {
      #Sys.sleep(.5)
    }
  }

  tmp <- pbreplicate(3, my_fun(10))

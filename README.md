
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RsienaTwoStep

<!-- badges: start -->
<!-- badges: end -->

The goal of `RsienaTwoStep` is to provide a method to asses the extent
to which results obtained by `RSiena::siena07` depend on the validity of
the ministep assumption.

A crucial assumption of the SAOM as implemented in RSiena is the
‘ministep’: only one actor at the time is allowed to make one
tie-change. This package aims to provide a method to asses the extent to
which results obtained by RSiena hinge on this crucial ministep
assumption. It does so by simulating the evolution of networks by
allowing two actors to make decisions on tie-changes simultaneously.

Depending on the precise implementation of the two-step procedure,
allowing for twosteps can be interpreted simply as:  
- allowing for simultaneity (when two actors are picked at random to
simultaneously make a ministep),  
- or as coordination (when there are either rules on which two actors
are allowed to make a twostep and/or which options after the twostep
count as coordination and are included in the choice set of the two
actors).

The package can mix the ratio of ministeps and twosteps.

## Installation

You can install the development version of RsienaTwoStep from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JochemTolsma/RsienaTwoStep")
```

## Example

ABM-ministep and ABM-twostep simulations on toy data

``` r
library(RsienaTwoStep)
#> Loading required package: foreach
## basic example code
```

### our network

``` r
net1g <- igraph::graph_from_adjacency_matrix(net1, mode="directed")
plot(net1g)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r
unloadNamespace("igraph") # to avoid clashes with package `sna`. 
```

### setting up cluster

This is of course not necessary but will make the simulations a lot
faster.

``` r
library(doParallel)
#> Loading required package: iterators
#> Loading required package: parallel
#stopCluster(my.cluster) 
n.cores <- parallel::detectCores() - 1  #save one core for other work
# create the cluster
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
# register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
```

## Degree and reciprocity

Let us assume people really don’t like to have a non-reciprocal tie but
do like reciprocal ties a lot. (I only want to help you if you help me!)

Thus if you start with a network without many (reciprocal) ties it would
be very difficult to get more reciprocal ties in the normal ministep
model. However, with simultaneity this should be possible.

### simulate networks for three conditions

``` r
sims1 <- f_sims(nsims=1000, parallel=TRUE, net=net1, rate=10, statistics=list(f_degree, f_recip), parameters=c(-1,2), p2step=0, chain=FALSE) #ministep only

sims2 <- f_sims(nsims=1000, parallel=TRUE, net=net1, rate=10, statistics=list(f_degree, f_recip), parameters=c(-1,2), p2step=1, chain=FALSE) #twostep, random selection of dyads

sims3 <- f_sims(nsims=1000, parallel=TRUE, net=net1, rate=10, statistics=list(f_degree, f_recip), parameters=c(-1,2), p2step=1, dist1=2, dist2=2, chain=FALSE) #twostep, but not all options are seen as coordination 
```

### counting dyads

``` r
nsims <- 1000

#combine results of dyad.census
df <- foreach(1:nsims, i=icount(), .combine="rbind") %dopar% {
  sna::dyad.census(sims3[[i]])
}
df <- as.data.frame(df)

#bit clumsy
df <- rbind(df, df, df)
df$x <- rep(c("mut", "asym", "null"), each=nsims)
df$y <- NA
df$y[df$x=="mut"] <- df$Mut[df$x=="mut"]
df$y[df$x=="asym"] <- df$Asym[df$x=="asym"]
df$y[df$x=="null"] <- df$Null[df$x=="null"]
df$type <- "twostep-coordination"
dftwostepNR <- df

df <- foreach(1:nsims, i=icount(), .combine="rbind") %dopar% {
  sna::dyad.census(sims2[[i]])
}
df <- as.data.frame(df)

#bit clumsy
df <- rbind(df, df, df)
df$x <- rep(c("mut", "asym", "null"), each=nsims)
df$y <- NA
df$y[df$x=="mut"] <- df$Mut[df$x=="mut"]
df$y[df$x=="asym"] <- df$Asym[df$x=="asym"]
df$y[df$x=="null"] <- df$Null[df$x=="null"]
df$type <- "twostep-simultaneity"
dftwostepR <- df

df <- foreach(1:nsims, i=icount(), .combine="rbind") %dopar% {
  sna::dyad.census(sims1[[i]])
}
df <- as.data.frame(df)

#bit clumsy
df <- rbind(df, df, df)
df$x <- rep(c("mut", "asym", "null"), each=nsims)
df$y <- NA
df$y[df$x=="mut"] <- df$Mut[df$x=="mut"]
df$y[df$x=="asym"] <- df$Asym[df$x=="asym"]
df$y[df$x=="null"] <- df$Null[df$x=="null"]
df$type <- "ministep"
dfministep <- df

## combine datasets

df <- rbind(dfministep, dftwostepR, dftwostepNR)
```

### plot results of the three dyadcensus

``` r
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.2.2

p <- ggplot(df, aes(x=x, y=y, fill=type)) + 
  geom_violin(position=position_dodge(1)) + 
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x),
               width=.1,
               color="red", position=position_dodge(1)) + 
  stat_summary(fun = mean,
               geom = "point",
               color="red", position=position_dodge(1)) 
  
p
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

<!--- 

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:


```r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN.

-->

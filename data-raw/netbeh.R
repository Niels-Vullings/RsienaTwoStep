## code to prepare `netbeh` dataset goes here
library(RSiena)

STATS <- list(ts_degree,
              ts_recip,
              ts_transTrip,
              ts_cycle3,
              list(ts_simX, "smoke"),
              list(ts_absDiffX, "alcohol"),
              ts_linear,
              ts_quad,
              list(ts_avAlt, "friendship"),
              list(ts_effFrom, "smoke"))

DF <- data.frame(alcohol.1 = s50a[, 1], alcohol.2 = s50a[, 2], smoke = s50s[, 1])

STARTS <- c(4,   # rate
            -1.5,  # outdegree
            2,   # reciprocity
            .8,  # transitive triplets
            .2,   # 3-cycles
            1.3, # smoke similarity,
            -1.7,# alcohol abs. difference
            1.4, # Rate beh
            .32, # linear shape
            -.4, # quadratic shape
            2.3, # average alter effect
            1.3  # effect from smoke
)

#strict coordination but at a relatively large distance of 3
set.seed(56784357)
ts_simnet1 <- ts_sims(
  nsims = 1,
  startvalues = STARTS,
  net1 = s501,
  ccovar = DF,
  statistics = STATS,
  p2step= c(2,8,0),
  dist1 = 3
)


net_sim1 <- ts_simnet1[[1]]$net_n
alcohol_sim1 <- as.numeric(ts_simnet1[[1]]$beh_n) +  attributes(ts_simnet1[[1]]$beh_n)$mean

sum(DF$alcohol.1)
sum(DF$alcohol.2)
sum(alcohol_sim1)
table(DF$alcohol.1, alcohol_sim1)
#let us simulate some random increase in alcohol.
alcohol_sim2 <- alcohol_sim1 + sample(c(0,1,2), 50, replace=TRUE, prob=c(7,2,1))
alcohol_sim2 <- ifelse(alcohol_sim2 > 5, 5, alcohol_sim2)
sum(alcohol_sim2)

DF_new <- data.frame(alcohol.1 = alcohol_sim1, alcohol.2 = alcohol_sim2, smoke = s50s[, 1])


set.seed(35684357)
ts_simnet2 <- ts_sims(
  nsims = 1,
  startvalues = STARTS,
  net1 = net_sim1,
  ccovar = DF_new,
  statistics = STATS,
  p2step= c(2,8,0),
  dist1 = 3
)

net_sim2 <- ts_simnet2[[1]]$net_n
alcohol_sim2 <- as.numeric(ts_simnet2[[1]]$beh_n) +  attributes(ts_simnet2[[1]]$beh_n)$mean
DF_sim <- data.frame(alcohol.1 = alcohol_sim1, alcohol.2 = alcohol_sim2, smoke = s50s[, 1])


ts_simdata <- list(net_sim1 = net_sim1, net_sim2 = net_sim2, DF_sim = DF_sim, startvalues_sim = STARTS, statistics_sim = STATS, p2_step_sim = c(0,1,0), dist1_sim = 3)
usethis::use_data(ts_simdata, overwrite = TRUE)

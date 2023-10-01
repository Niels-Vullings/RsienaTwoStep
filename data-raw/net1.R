## code to prepare `net1` dataset goes here
set.seed(56787457)
ts_net1 <- matrix(sample(c(0,1), 100, replace = TRUE, prob = c(0.7, 0.1)), nrow=10, ncol=10)
diag(ts_net1) <- 0
#include an isolate
ts_net1[5,10] <- 0
#include some reciprocity
ts_net1[4,9] <- 1
ts_net1[8,10] <- 1
usethis::use_data(ts_net1, overwrite = TRUE)

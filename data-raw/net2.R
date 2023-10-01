## code to prepare `net2` dataset goes here
set.seed(56787457)
ts_net2 <- matrix(sample(c(0,1), 25, replace = TRUE, prob = c(0.7, 0.1)), nrow=5, ncol=5)
diag(ts_net2) <- 0
usethis::use_data(ts_net2, overwrite = TRUE)


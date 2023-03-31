## code to prepare `df_ccovar1` dataset goes here
set.seed(56787457)
cov1 <- sample(c(0,1), 10, replace=T)
cov2 <- round(rnorm(10,2,3))
df_ccovar1 <- data.frame(cov1=cov1, cov2=cov2)
usethis::use_data(df_ccovar1, overwrite = TRUE)

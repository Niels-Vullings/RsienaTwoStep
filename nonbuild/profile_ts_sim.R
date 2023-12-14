Rprof("./nonbuild/Rprof.txt", append=FALSE, memory.profiling = TRUE, line.profiling = TRUE, numfiles = 100L)

ts_estim(ans = ans1,
         p2step = c(0,1,0),
         nite = 30,
         itef1 = 3,
         itef3 = 30,
         dist1 = 2,
         phase1 = TRUE,
         phase3 = FALSE,
         parallel = FALSE,
         verbose = TRUE)

  Rprof(NULL)



  summaryRprof("./nonbuild/Rprof.txt", lines = "both")

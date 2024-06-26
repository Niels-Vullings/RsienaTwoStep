#----Dyad Function----
ts_dyad.evo <- function(t1 = "network",t2 = "network or simulated data" , sim_sample = c(1:1000), simN = 1000, seed = 345657343, save.dyads = FALSE){

  if((isTRUE("simnets" %in% names(t2)) == TRUE)){ # check the datatype of T2

    t2 <- t2$simnets # allow option to use entire simulated data as input

  }


  save <- save.dyads #savethe boolean for later

  #inner function for dyad census
  dyad.census.type <- function(t1, t2, save.dyads = FALSE){

    t1 <- t1 #+ 1
    diag(t1) = NA #exclude the diagonal, not relevant data
    t2 <- t2
    diag(t2) = NA #exclude the diagonal, not relevant data

    flips <- t1 + t(t2) - t2 #subtract s502 to ensure that Mutual ties are not included


    jumpst1 <- t1 + t(t1) #t1 plus its transpose lead to a value of 2 for mutual ties
    stablet1 <- jumpst1 #use this for stable assymetric ties
    jumpst1[lower.tri(jumpst1)] <- NA #remove duplicate ties

    jumpst2 <- t2 + t(t2)
    stablet2 <- jumpst2 #use this for stable assymetric ties
    jumpst2[lower.tri(jumpst2)] <- NA #remove duplicate ties

    #Dyad combinations
    stable00 <- as.data.frame(which(jumpst1 == 0 & jumpst2 == 0, arr.ind = TRUE)) #Null at t1 and Null at t2
    stable01 <- as.data.frame(which(flips == 0 & stablet1 == 1 & stablet2 == 1, arr.ind = TRUE)) #Assymetric at t1 and t2 and no flip
    stable11 <- as.data.frame(which(jumpst1 == 2 & jumpst2 == 2, arr.ind = TRUE)) #Mutual at t1 and Mutual at t2

    Null_Assym <- as.data.frame(which(jumpst1 == 0 & stablet2 == 1, arr.ind = TRUE))#Null -> assym
    Assym_Null <- as.data.frame(which(flips == 0 & stablet1 == 1 & stablet2 == 0, arr.ind = TRUE))#Assum -> Null
    Assym_Mut <- as.data.frame(which(stablet1 == 1 & jumpst2 == 2, arr.ind = TRUE))#Assym -> Mutual
    Mut_Assym <- as.data.frame(which(jumpst1 == 2 & jumpst2 == 1, arr.ind = TRUE))#Assym -> Mutual

    flip <- as.data.frame(which(flips == 2 & stablet1 != 2 & stablet2 != 2, arr.ind = TRUE)) #flips, so 01 at T1 and 10 at T2
    jump02 <- as.data.frame(which(jumpst1 == 0 & jumpst2 == 2, arr.ind = TRUE)) #Null jump, from Null to Mutual
    jump20 <- as.data.frame(which(jumpst1 == 2 & jumpst2 == 0, arr.ind = TRUE)) #Mutual jump, from Mutual to Null


    table <- cbind(nrow(stable00),
                   nrow(stable01),
                   nrow(stable11),
                   nrow(Null_Assym),
                   nrow(Assym_Null),
                   nrow(Assym_Mut),
                   nrow(Mut_Assym),
                   nrow(flip),
                   nrow(jump02),
                   nrow(jump20)) #bind the rownumbers to determine the census of the tie type

    colnames(table) = c("Null > Null",
                        "Assym > Assym",
                        "Mutual > Mutual",
                        "Null > Assym",
                        "Assym > Null",
                        "Assym > Mutual",
                        "Mutual > Assym",
                        "Tie flip",
                        "Null > Mutual",
                        "Mutual > Null")

    if(save.dyads == TRUE){
      #give each dyad a relationship identifier
      stable00$tie <- ifelse(is.na(stable00$row),NA, "Null > Null")
      stable01$tie <- ifelse(is.na(stable01$row),NA, "Assym > Assym")
      stable11$tie <- ifelse(is.na(stable11$row),NA, "Mutual > Mutual")
      Null_Assym$tie <- ifelse(is.na(Null_Assym$row),NA, "Null > Assym")
      Assym_Null$tie <- ifelse(is.na(Assym_Null$row),NA, "Assym > Null")
      Assym_Mut$tie <- ifelse(is.na(Assym_Mut$row),NA, "Assym > Mutual")
      Mut_Assym$tie <- ifelse(is.na(Mut_Assym$row),NA, "Mutual > Assym")
      flip$tie <- ifelse(is.na(flip$row),NA, "Tie flip")
      jump02$tie <- ifelse(is.na(jump02$row),NA, "Null > Mutual")
      jump20$tie <- ifelse(is.na(jump20$row),NA, "Mutual > Null")

      dyads <<- rbind(stable00,stable01,stable11, Null_Assym, Assym_Null, Assym_Mut, Mut_Assym, flip,jump02,jump20) #save a dataframe which contains the relationships of each dyad

      # print("Dyad transitions from T1 to T2")
      return(as.data.frame(table))

    } else{

      # print("Dyad transitions from T1 to T2")
      return(as.data.frame(table))


    }
  }

  if(is.list(t2)){

    set.seed(seed) #set seed for similar outcomes with the simulations
    if(simN > length(t2)){
      print("Warning, simN higher than the number of simulations")
    } else{

      sample <- as.vector(sample(sim_sample, simN, replace = TRUE)) #create a vector that chooses simulations at random
    }


    df <- foreach::foreach(sim_sample, i=iterators::icount(), .combine="rbind") %dopar% {

      df <- as.data.frame(dyad.census.type(t1, t2= t2[[i]])) #loop the dyad.census.type function over all relevant simulations and save in df

    }

    assign("Simulated ties", df, envir = .GlobalEnv) #save the census df to the global environment

    output <- round(colMeans(df)) #Take the mean of the dyad occurrences to get the census mean

    return(output)

  } else{
    dyad.census.type(t1,t2, save.dyads = save)
  }
}


#----Triad Function----
ts_triads.evo <- function(t1,t2, filter = FALSE) {

  #-- Check required packages --
  fpackage.check <- function(packages) {
    lapply(packages, FUN = function(x) {
      if (x %in% .packages() == FALSE) {
        warning(paste0("please load package:", x))

      }
    })
  }
  packages <- c("foreach", "sna", "iterators", "doParallel")
  fpackage.check(packages)

  #-- Check how large the network is --
  if(length(t1) >= 1225){

    print("Large networks can cause higher runtimes")
  }

  #-- Check if Parallel Cluster is loaded in --
  if(foreach::getDoParRegistered() == FALSE){

    warning("Function uses doParallel: No parallel backend registered")  #Check if a parallel cluster has been registered and return warning if FALSE

  } else{

    if(length(t1) == length(t2)){ #Check if matrices at T1 and T2 have equal dimmensions

      #Davis & Leinhardt triad classification for identifying triads
      triads <- c("X003", "X012", "X102", "X021D",  "X021U",  "X021C",  "X111D",  "X111U",  "X030T",  "X030C",   "X201",  "X120D",  "X120U",  "X120C",   "X210",   "X300")

      # Create a dataframe(df) that determines a triad census for each possible triad in the matrix
      df <- foreach(a1=1:nrow(t1), i=icount() , .combine="rbind") %:%
        foreach(a2=1:nrow(t1), j=icount() , .combine="rbind") %:%
        foreach(a3=1:nrow(t1),  k=icount() , .combine="rbind") %dopar% {

          if (i>j & j>k ) data.frame(i=i, j=j, k=k, #(i>j & j>k ) determines with or without repititions [current = NO REP]
                                     t1_ij = as.character(t1[a1,a2]), t1_ji = as.character(t1[a2,a1]),
                                     t1_ik = as.character(t1[a1,a3]),t1_ki = as.character(t1[a3,a1]),
                                     t1_jk = as.character(t1[a2,a3]),t1_kj = as.character(t1[a3,a2]), #Configuration of T1 triad
                                     typeT1 = triads[which(sna::triad.census(t1[c(a1,a2,a3),c(a1,a2,a3) ]) == 1) ], # Census of specific triad, based on the values that are in a1,a2,a3 for T1
                                     t2_ij = as.character(t2[a1,a2]), t2_ji = as.character(t2[a2,a1]),
                                     t2_ik = as.character(t2[a1,a3]),t2_ki = as.character(t2[a3,a1]),
                                     t2_jk = as.character(t2[a2,a3]),t2_kj = as.character(t2[a3,a2]), # Configuration of T2 triad
                                     typeT2 = triads[which(sna::triad.census(t2[c(a1,a2,a3),c(a1,a2,a3) ]) == 1) ], # Census of specific triad, based on the values that are in a1,a2,a3 for T2
                                     equal = sum(c(t1[a1,a2],t1[a2,a1],t1[a1,a3],t1[a3,a1],t1[a2,a3],t1[a3,a2]) != c(t2[a1,a2],t2[a2,a1],t2[a1,a3],t2[a3,a1],t2[a2,a3],t2[a3,a2])))
        }

      df$name <- paste0(df$i,".",df$j,".",df$k) #name is each actor in triad [actor1.actor2.actor3]

      df <- subset(df, select=c(name, t1_ij:t1_kj, typeT1, t2_ij:t2_kj, typeT2, equal)) # Subset to remove redundant information

      print(as.data.frame(df[which(df$equal > 2),c("name", "typeT1", "typeT2", "equal")])) #prints the "invalid" triad changes

      triads <<- df # Write out all triads for visual inspection

      #Start printing output
      #Based on the established parameter, return the table for filtered (filter == TRUE) or unfiltered (filter == FALSE) data
      if(filter == TRUE){

        filt <- df[df$equal <=2,] #Filters data to only include rows that classify as "valid"
        df <- reshape2::melt(table(filt$typeT1, filt$typeT2, dnn = c("Timepoint.1", "Timepoint.2")))
        df$percentage <- round(df$value/sum(df$value), digits = 3)*100
        print(as.data.frame(df[order(df$Timepoint.1, decreasing = FALSE),]))

      } else{

        df <- reshape2::melt(table(df$typeT1, df$typeT2, dnn = c("Timepoint.1", "Timepoint.2")))
        df$percentage <- round(df$value/sum(df$value), digits = 3)*100
        print(as.data.frame(df[order(df$Timepoint.1, decreasing = FALSE),]))

      }


    } else{

      warning("Execution halted: Matrices are not equal")

    }
  }
}

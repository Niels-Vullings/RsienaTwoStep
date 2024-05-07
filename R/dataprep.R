#' Data preparation
#'
#' @description
#' `ts_prepdata` performs centering and similarity score and set attribute prepared to TRUE
#' `ts_centering` centers the variables before use.
#' `ts_simij` calculates the similarity scores before use.
#'
#' @details
#' I really need to update the dataprep part, so to have behavioral dependents
#' ccovars and time varying covars.
#' @param cov numeric, behavioral scores of actors
#' @param min numeric, minimum value of behavioral scores of actors. If `NULL` the empirically observed minimum is used.
#' @param max numeric, maximum value of behavioral scores of actors. If `NULL` the empirically observed maximum is used.
#' @param ccovar data frame with named time-constant covariates.
#'
#' @examples
#' ts_centering(cov=df_ccovar1[,"cov1"])
#' ts_simij(cov=df_ccovar1[,"cov2"])
#' ts_simij(cov=df_ccovar1[,"cov2"], min=-5, max=5)
#' @importFrom Rdpack reprompt
#' @importFrom stringr str_sub
#' @export
ts_centering <- function(cov) {
  attr(cov, "mean") <- round(mean(cov), 4)
  centered <- cov - attributes(cov)$mean
  return(centered)
}

#' @rdname ts_centering
#' @export
ts_prep_dep <- function(cov) {
  #mean and range calculated on the two waves
  attr(cov[,1], "mean") <- round(mean(c(cov[,1],cov[,2])), 4)
  attr(cov[,2], "mean") <- round(mean(c(cov[,1],cov[,2])), 4)
  max <- max(c(cov[,1],cov[,2]), na.rm=TRUE)
  min <- min(c(cov[,1],cov[,2]), na.rm=TRUE)
  rv <- max - min

  #attach range/range2
  attr(cov[,1], "range") <- rv
  attr(cov[,1], "range2") <- c(min, max)
  attr(cov[,2], "range") <- rv
  attr(cov[,2], "range2") <- c(min, max)

  #center depvars on total mean
  cov[,1] <- cov[,1] - attributes(cov[,1])$mean
  cov[,2] <- cov[,2] - attributes(cov[,2])$mean

  #calculate similarity scores
  mat1 <- matrix(cov[,1], nrow = length(cov[,1]), ncol = length(cov[,1]), byrow = TRUE)
  mat2 <- t(mat1)
  simij <- 1 - (abs(mat1 - mat2)/rv)
  diag(simij) <- NA
  attr(cov[,1], "simMean") <- mean(simij, na.rm=TRUE)
  attr(cov[,1], "simij") <- simij

  mat1 <- matrix(cov[,2], nrow = length(cov[,2]), ncol = length(cov[,2]), byrow = TRUE)
  mat2 <- t(mat1)
  simij <- 1 - (abs(mat1 - mat2)/rv)
  diag(simij) <- NA
  attr(cov[,2], "simMean") <- mean(simij, na.rm=TRUE)
  attr(cov[,2], "simij") <- simij


  return(cov)
}

#' @rdname ts_centering
#' @export
ts_simij <- function(cov, min=NULL, max=NULL) {
  if (is.null(min) & is.null(max)) {
    max <- max(cov, na.rm=TRUE)
    min <- min(cov, na.rm=TRUE)
    rv <- max - min
  }
  if (!is.null(min) & !is.null(max)) rv <- max - min
  mat1 <- matrix(cov, nrow = length(cov), ncol = length(cov), byrow = TRUE)
  mat2 <- t(mat1)
  simij <- 1 - (abs(mat1 - mat2)/rv)
  diag(simij) <- NA
  attr(cov, "simMean") <- mean(simij, na.rm=TRUE)
  attr(cov, "range") <- rv
  attr(cov, "range2") <- c(min, max) + attributes(cov)$mean
  attr(cov, "simij") <- simij
  return(cov)
}

#' @rdname ts_centering
#' @export
ts_prepdata <- function(ccovar) {
  if (!is.null(ccovar) & is.null(attributes(ccovar)$prepared)) {
    # check if there are dependent variables stored in ccovar
    if (stringr::str_sub(names(ccovar)[1], -2, -1) == ".1") {
      nodeps <- FALSE
    } else {
      nodeps <- TRUE
    }

    if (nodeps) {
      for (i in 1:ncol(ccovar)) {
        ccovar[, i] <- ts_centering(ccovar[, i])
        ccovar[, i] <- ts_simij(ccovar[, i])
      }
    } else {
      ccovar[, 1:2] <- ts_prep_dep(ccovar[, 1:2])
      for (i in 3:ncol(ccovar)) {
        ccovar[, i] <- ts_centering(ccovar[, i])
        ccovar[, i] <- ts_simij(ccovar[, i])
      }
    }
    attr(ccovar, "prepared") <- TRUE
  }
  return(ccovar)
}

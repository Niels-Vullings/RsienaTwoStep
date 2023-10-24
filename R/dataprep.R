#' Data preparation
#'
#' @description
#' `ts_dataprep` performs centering and similarity score and set attribute prepared to TRUE
#' `ts_centering` centers the variables before use.
#' `ts_simij` calculates the similarity scores before use.
#'
#' @details
#' to do
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
#' @export
ts_centering <- function(cov) {
  attr(cov, "mean") <- mean(cov)
  centered <- cov - attributes(cov)$mean
  return(centered)
}

#' @rdname ts_centering
#' @export
ts_simij <- function(cov, min=NULL, max=NULL) {
  if (is.null(min) & is.null(max)) rv <- max(cov, na.rm=TRUE) - min(cov, na.rm=TRUE)
  if (!is.null(min) & !is.null(max)) rv <- max - min
  mat1 <- matrix(cov, nrow = length(cov), ncol = length(cov), byrow = TRUE)
  mat2 <- t(mat1)
  simij <- 1 - (abs(mat1 - mat2)/rv)
  diag(simij) <- NA
  attr(cov, "simMean") <- mean(simij, na.rm=TRUE)
  attr(cov, "range") <- rv
  attr(cov, "simij") <- simij
  return(cov)
}

#' @rdname ts_centering
#' @export
ts_prepdata <- function(ccovar) {
  if (!is.null(ccovar) & is.null(attributes(ccovar)$prepared)) {
    for (i in 1:ncol(ccovar)) {
      ccovar[, i] <- ts_centering(ccovar[, i])
      ccovar[, i] <- ts_simij(ccovar[, i])
    }
    attr(ccovar, "prepared") <- TRUE
  }
  return(ccovar)
}

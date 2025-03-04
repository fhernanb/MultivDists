#' Bivariate Poisson Distribution - Geoffroy
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function obtains the probability for the Bivariate Poisson distribution
#' under the parameterization of Geoffroy et. al (2021).
#'
#' @param x vector or matrix of quantiles. When \code{x} is a matrix,
#' each row is taken to be a quantile and columns correspond to the number of dimensions \code{p}.
#' @param n number of random observations.
#' @param l1 mean for \eqn{Z_1} variable with Poisson distribution.
#' @param l2 mean for \eqn{Z_2} variable with Poisson distribution.
#' @param l0 mean for \eqn{U} variable with Poisson distribution.
#' @param log logical; if \code{TRUE}, densities d are given as log(d).
#'
#' @returns
#' Returns the density for a given data \code{x}.
#'
#' @references
#' Kouakou, K. J. G., Hili, O., & Dupuy, J. F. (2021). Estimation in the zero-inflated bivariate Poisson model with an application to health-care utilization data. Afrika Statistika, 16(2), 2767-2788.
#'
#' @example examples/examples_dBP_Geoffroy.R
#'
#' @export
#' @useDynLib MultivDists
#' @importFrom Rcpp sourceCpp
dBP_Geoffroy <- function(x, l1, l2, l0, log=FALSE) {
  # Initial checks
  if(any(l1 <= 0)) stop("lambda_1 must be positive.")
  if(any(l2 <= 0)) stop("lambda_2 must be positive.")
  if(any(l0 <= 0)) stop("lambda_0 must be positive.")

  if(is.vector(x)) {
    if (length(x) !=2) stop("The vector length must be 2.")
    x <- matrix(x, nrow=1)
  }

  ## Begin Auxiliar function
  phi <- function(x, l1, l2, l0) {
    s <- 0:min(x)
    num <- l0^s*l1^(x[1]-s)*l2^(x[2]-s)
    den <- factorial(s) * factorial(x[1]-s) * factorial(x[2]-s)
    res <- sum(num/den)
    res
  }
  ## End Auxiliar function

  res <- apply(X=x, MARGIN=1, FUN=phi, l1=l1, l2=l2, l0=l0)
  res <- exp(-l1-l2-l0) * res
  if (log)
    log(res)
  else
    res
}
#'
#' @rdname dBP_Geoffroy
#' @importFrom stats rpois cov
#' @export
#' @useDynLib MultivDists
#' @importFrom Rcpp sourceCpp
rBP_Geoffroy <- function(n, l1, l2, l0) {

  z1 <- rpois(n=n, lambda=l1)
  z2 <- rpois(n=n, lambda=l2)
  u  <- rpois(n=n, lambda=l0)
  x1 <- z1 + u
  x2 <- z2 + u
  x <- cbind(x1, x2)
  colnames(x) <- c("X1", "X2")
  res <- as.matrix(x)
  return(x)
}
#' Moment estimations for Bivariate Poisson Distribution - Geoffroy
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function obtains moment estimators for the Bivariate Poisson distribution
#' under the parameterization of Geoffroy et. al (2021).
#'
#' @param x vector or matrix of quantiles. When \code{x} is a matrix,
#' each row is taken to be a quantile and columns correspond to the number of dimensions \code{p}.
#'
#' @returns
#' Returns a vector with \eqn{\hat{\lambda_1}}, \eqn{\hat{\lambda_2}} and \eqn{\hat{\alpha}}.
#'
#' @references
#' Kouakou, K. J. G., Hili, O., & Dupuy, J. F. (2021). Estimation in the zero-inflated bivariate Poisson model with an application to health-care utilization data. Afrika Statistika, 16(2), 2767-2788.
#'
#' @example examples/examples_dBP_Geoffroy.R
#'
#' @importFrom stats cov
#' @export
moments_estim_BP_Geoffroy <- function(x) {

  # For l0
  l0_hat <- abs(cov(x)[2])

  # For l1
  l1_hat <- mean(x[,1]) - l0_hat

  # For l1
  l2_hat <- mean(x[,2]) - l0_hat

  theta_hat <- c(l1_hat=l1_hat, l2_hat=l2_hat, l0_hat=l0_hat)
  return(theta_hat)
}


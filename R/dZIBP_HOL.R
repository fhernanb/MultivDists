#' Zero Inflated Bivariate Poisson Distribution - Holgate
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function obtains the probability for the Zero Inflated
#' Bivariate Poisson distribution
#' under the parameterization of Kouakou et. al (2021).
#'
#' @param x vector or matrix of quantiles. When \code{x} is a matrix,
#' each row is taken to be a quantile and columns correspond to the number of dimensions \code{p}.
#' @param n number of random observations.
#' @param l1 mean for \eqn{Z_1} variable with Poisson distribution.
#' @param l2 mean for \eqn{Z_2} variable with Poisson distribution.
#' @param l0 mean for \eqn{U} variable with Poisson distribution.
#' @param psi parameter with the contamination proportion, \eqn{0 \leq \psi \leq 1}.
#' @param log logical; if \code{TRUE}, densities d are given as log(d).
#'
#' @returns
#' Returns the density for a given data \code{x}.
#'
#' @references
#' P. Holgate. Estimation for the bivariate poisson distribution.
#' Biometrika, 51(1/2):241–245, 1964. ISSN 00063444, 14643510.
#' URL http://www.jstor.org/stable/2334210.
#'
#' Konan Jean Geoffroy Kouakou, Ouagnina Hili, and Jean-Francois
#' Dupuy. Estimation in the zero-inflated bivariate poisson model
#' with an application to health-care utilization data.
#' Afrika Statistika, 16(2):2767–2788, 2021.
#'
#' @example examples/examples_dZIBP_HOL.R
#'
#' @export
#' @useDynLib MultivDists
#' @importFrom Rcpp sourceCpp
dZIBP_HOL <- function(x, l1, l2, l0, psi, log=FALSE) {
  # Initial checks
  if(any(l1 <= 0)) stop("lambda_1 must be positive.")
  if(any(l2 <= 0)) stop("lambda_2 must be positive.")
  if(any(l0 <= 0)) stop("lambda_0 must be positive.")
  if(any(psi<0 | psi>1)) stop("psi must be in the (0, 1) interval.")

  if(is.vector(x)) {
    if (length(x) != 2) stop("The vector length must be 2.")
    x <- matrix(x, nrow=1)
  }

  ## Begin Auxiliar function
  aux_pmf <- function(input) {
    y1 <- input[1]
    y2 <- input[2]
    l1 <- input[3]
    l2 <- input[4]
    l0 <- input[5]
    psi <- input[6]
    if (y1==0 & y2==0) {
      res <- psi + (1-psi) * exp(-l0-l1-l2)
    }
    else {
      res <- (1-psi)*exp(-l0-l1-l2)*phi(c(y1, y2), l1=l1, l2=l2, l0=l0)
    }
    return(res)
  }

  phi <- function(x, l1, l2, l0) {
    s <- 0:min(x)
    num <- l0^s*l1^(x[1]-s)*l2^(x[2]-s)
    den <- factorial(s) * factorial(x[1]-s) * factorial(x[2]-s)
    res <- sum(num/den)
    res
  }
  ## End Auxiliar function

  res <- apply(X=cbind(x, l1, l2, l0, psi), MARGIN=1, FUN=aux_pmf)
  if (log)
    log(res)
  else
    res
}
#'
#' @rdname dZIBP_HOL
#' @importFrom stats runif
#' @export
#' @useDynLib MultivDists
#' @importFrom Rcpp sourceCpp
rZIBP_HOL <- function(n, l1, l2, l0, psi) {
  if(any(l1 <= 0)) stop("lambda_1 must be positive.")
  if(any(l2 <= 0)) stop("lambda_2 must be positive.")
  if(any(l0 <= 0)) stop("lambda_0 must be positive.")
  if(any(psi<0 | psi>1)) stop("psi must be in the (0, 1) interval.")
  if(n <= 0) stop("n must be a positive integer")

  i <- sample(x=c(TRUE, FALSE), size=n, prob=c(psi, 1-psi), replace=TRUE)
  z1 <- rpois(n=n, lambda=l1)
  z2 <- rpois(n=n, lambda=l2)
  u  <- rpois(n=n, lambda=l0)
  X1 <- z1 + u
  X2 <- z2 + u
  x <- cbind(X1, X2)
  x[i, ] <- 0
  return(x)
}
#'
#'
#'
#' Moment estimations for Zero Inflated Bivariate Poisson
#' Distribution - Holgate
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function obtains moment estimators for the Zero Inflated
#' Bivariate Poisson distribution
#' under the parameterization of Holgate (1964).
#'
#' @param x vector or matrix of quantiles. When \code{x} is a matrix,
#' each row is taken to be a quantile and columns correspond to the number of dimensions \code{p}.
#'
#' @returns
#' Returns a vector with \eqn{\hat{\lambda_1}}, \eqn{\hat{\lambda_2}}, \eqn{\hat{\alpha}} and \eqn{\hat{\psi}}.
#'
#' @references
#' P. Holgate. Estimation for the bivariate poisson distribution.
#' Biometrika, 51(1/2):241–245, 1964. ISSN 00063444, 14643510.
#' URL http://www.jstor.org/stable/2334210.
#'
#' @example examples/examples_dZIBP_HOL.R
#'
#' @importFrom stats cor var
#' @importFrom nleqslv nleqslv
#' @export
moments_estim_ZIBP_HOL <- function(x) {

  # This is an auxiliar function to create a equation system
  # with the samples statistics
  aux_fun <- function(theta) {
    l1 <- theta[1]
    l2 <- theta[2]
    l0 <- theta[3]
    psi <- theta[4]
    EX1 <- colMeans(x)[1]
    EX2 <- colMeans(x)[2]
    z <- numeric(4) # contiene LD de las ecuaciones
    z[1] <- EX1 - (1-psi) * (l1+l0)
    z[2] <- var(x)[1] - EX1 * (1+psi*(l1+l0))
    z[3] <- var(x)[4] - EX2 * (1+psi*(l2+l0))
    z[4] <- var(x)[2] - (1-psi) * (l0+psi*(l1+l0)*(l2+l0))
    return(z)
  }

  res <- nleqslv::nleqslv(x=c(1, 1, 1, 0.5), fn=aux_fun,
                          method="Newton",
                          control=list(btol=0.01))
  theta_hat <- res$x
  names(theta_hat) <- c("l1_hat", "l2_hat", "l0_hat", "psi_hat")
  return(theta_hat)
}


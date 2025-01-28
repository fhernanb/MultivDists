#' Bivariate Poisson Distribution - Laksh
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function obtains the probability for the Bivariate Poisson distribution
#' under the parameterization of Lakshminarayana et. al (1993).
#'
#' @param x vector or matrix of quantiles. When \code{x} is a matrix,
#' each row is taken to be a quantile and columns correspond to the number of dimensions \code{p}.
#' @param l1 mean for the marginal \eqn{X_1} variable with Poisson distribution.
#' @param l2 mean for the marginal \eqn{X_2} variable with Poisson distribution.
#' @param alpha third parameter.
#' @param log logical; if \code{TRUE}, densities d are given as log(d).
#'
#' @returns
#' Returns the density for a given data \code{x}.
#'
#' @example examples/examples_dBP_Laksh.R
#'
#' @export
#' @useDynLib MultivDists
#' @importFrom Rcpp sourceCpp
dBP_Laksh <- function(x, l1, l2, alpha, log=FALSE) {
  # Initial checks
  if(any(l1 <= 0)) stop("lambda_1 must be positive.")
  if(any(l2 <= 0)) stop("lambda_2 must be positive.")

  if(is.vector(x)) {
    if (length(x) !=2) stop("The vector length must be 2.")
    x <- matrix(x, nrow=1)
  }

  augmented_x <- cbind(x, l1, l2, alpha)

  res <- apply(X=augmented_x, MARGIN=1, FUN=aux_BP_Laksh)

  if (log)
    res
  else
    exp(res)
}
#'
#'
#'
#' Correct value for alpha parameter - Lakshminarayana
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function obtains the interval in which the parameter \eqn{\alpha} must be.
#'
#' @param l1 mean for the marginal \eqn{X_1} variable with Poisson distribution.
#' @param l2 mean for the marginal \eqn{X_2} variable with Poisson distribution.
#'
#' @returns
#' Returns the density for a given data \code{x}.
#'
#' @examples
#' correct_alpha_BP_Laksh(l1=3, l2=4)
#' correct_alpha_BP_Laksh(l1=4, l2=3)
#' correct_alpha_BP_Laksh(l1=4, l2=4)
#'
#' @export
correct_alpha_BP_Laksh <- function(l1, l2) {
  k <- 1 - exp(-1)
  A <- exp(-l1*k)
  B <- exp(-l2*k)
  max_alpha <-  1 / ((1-A)*(1-B))
  min_alpha <- -1 / ((1-A)*(1-B))
  return(list(min_alpha=min_alpha, max_alpha=max_alpha))
}

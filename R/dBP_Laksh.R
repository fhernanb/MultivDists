#' Bivariate Poisson Distribution - Lakshminarayana
#'
#' @param x vector or matrix of quantiles. When \code{x} is a matrix, each row is taken to be a quantile and columns correspond to the number of dimensions \code{p}.
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

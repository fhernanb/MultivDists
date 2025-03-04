#' Zero Inflated Bivariate Poisson Distribution - Laksh
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function obtains the probability for the Zero Inflated
#' Bivariate Poisson distribution
#' under the parameterization of Lakshminarayana et. al (1993).
#'
#' @param x vector or matrix of quantiles. When \code{x} is a matrix,
#' each row is taken to be a quantile and columns correspond to the number of dimensions \code{p}.
#' @param n number of random observations.
#' @param l1 mean for the marginal \eqn{X_1} variable with Poisson distribution.
#' @param l2 mean for the marginal \eqn{X_2} variable with Poisson distribution.
#' @param alpha third parameter.
#' @param psi parameter with the contamination proportion, \eqn{0 \leq \psi \leq 1}.
#' @param max_val_x1 maximum value for \eqn{X_1} that is expected, by default it is 100.
#' @param max_val_x2 maximum value for \eqn{X_2} that is expected, by default it is 100.
#' @param log logical; if \code{TRUE}, densities d are given as log(d).
#'
#' @returns
#' Returns the density for a given data \code{x}.
#'
#' @references
#' Lakshminarayana, J., Pandit, S. N., & Srinivasa Rao, K. (1999). On a bivariate Poisson distribution. Communications in Statistics-Theory and Methods, 28(2), 267-276.
#'
#' @example examples/examples_dZIBP_Laksh.R
#'
#' @export
#' @useDynLib MultivDists
#' @importFrom Rcpp sourceCpp
dZIBP_Laksh <- function(x, l1, l2, alpha, psi, log=FALSE) {
  # Initial checks
  if(any(l1 <= 0)) stop("lambda_1 must be positive.")
  if(any(l2 <= 0)) stop("lambda_2 must be positive.")
  if(any(psi<0 | psi>1)) stop("psi must be in [0, 1].")

  if(is.vector(x)) {
    if (length(x) != 2) stop("The vector length must be 2.")
    x <- matrix(x, nrow=1)
  }

  augmented_x <- cbind(x, l1, l2, alpha, psi)
  res <- apply(X=augmented_x, MARGIN=1, FUN=aux_ZIBP_Laksh)

  if (log)
    res
  else
    exp(res)
}
#'
#' @rdname dZIBP_Laksh
#' @importFrom stats runif
#' @export
#' @useDynLib MultivDists
#' @importFrom Rcpp sourceCpp
rZIBP_Laksh <- function(n, l1, l2, alpha, psi,
                        max_val_x1=NULL, max_val_x2=NULL) {

  if(any(l1 <= 0)) stop("lambda_1 must be positive.")
  if(any(l2 <= 0)) stop("lambda_2 must be positive.")
  if(any(psi<0 | psi>1)) stop("psi must be in [0, 1].")
  if(n <= 0) stop("n must be a positive integer")

  if(is.null(max_val_x1)) max_val_x1 <- 100
  if(is.null(max_val_x2)) max_val_x2 <- 100

  val_x1 <- 0:max_val_x1
  val_x2 <- 0:max_val_x2
  grid_x <- expand.grid(val_x1, val_x2)
  grid_x <- as.matrix(grid_x)
  probs <- dZIBP_Laksh(x=grid_x, l1=l1, l2=l2, alpha=alpha, psi=psi)
  distri <- data.frame(grid_x, probs)
  distri <- distri[order(-distri$probs), ] # Organizando
  distri$cumul_probs <- cumsum(distri$probs)

  # Function to extract positions in a probability table
  pseudo_r_BP_Laksh <- function(u, cumul_probs) {
    position <- NULL
    for (i in 1:length(u)) {
      position[i] <- which(u[i] <= cumul_probs)[1]
    }
    return(position)
  }

  # Generating n values from a BP
  positions <- pseudo_r_BP_Laksh(u=runif(n),
                                 cumul_probs=distri$cumul_probs)
  res <- distri[positions, 1:2]

  colnames(res) <- c("X1", "X2")

  rownames(res) <- NULL
  res <- as.matrix(res)
  return(res)
}
#'
#'
#'
#' Moment estimations for Zero Inflated Bivariate Poisson Distribution - Laksh
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function obtains moment estimators for the Zero Inflated
#' Bivariate Poisson distribution
#' under the parameterization of Lakshminarayana et. al (1993).
#'
#' @param x vector or matrix of quantiles. When \code{x} is a matrix,
#' each row is taken to be a quantile and columns correspond to the number of dimensions \code{p}.
#'
#' @returns
#' Returns a vector with \eqn{\hat{\lambda_1}}, \eqn{\hat{\lambda_2}}, \eqn{\hat{\alpha}} and \eqn{\hat{\psi}}.
#'
#' @references
#' Lakshminarayana, J., Pandit, S. N., & Srinivasa Rao, K. (1999). On a bivariate Poisson distribution. Communications in Statistics-Theory and Methods, 28(2), 267-276.
#'
#' @example examples/examples_moments_estim_ZIBP_Laksh.R
#'
#' @importFrom stats cor
#' @export
moments_estim_ZIBP_Laksh <- function(x) {
  n <- nrow(x)

  # zero positions
  positions_zeros <- x[, 1] == 0 & x[, 2] == 0
  n_zero <- sum(positions_zeros)

  # subset with no zeros
  x_without_zeros <- x[!positions_zeros, ]

  # Correlation value without zeros
  correlation <- cor(x_without_zeros)[1, 2]

  # Estimating l1
  l1_hat <- mean(x_without_zeros[, 1])

  # Estimating l2
  l2_hat <- mean(x_without_zeros[, 2])

  psi_hat <- 1 - mean(x[, 1]) / l1_hat

  # subset discounting excess zeros
  small_zeros <- matrix(0, ncol=2, nrow=n_zero * psi_hat)
  small_zeros <- as.data.frame(small_zeros)
  colnames(small_zeros) <- colnames(x)
  x_discounting_zeros <- rbind(x_without_zeros,
                               small_zeros)

  # Auxiliar functio to check if alpha_hat is inside the interval
  aux_check <- function(alpha_hat, l1_hat, l2_hat) {
    min_alpha <- correct_alpha_BP_Laksh(l1=l1_hat, l2=l2_hat)$min_alpha
    max_alpha <- correct_alpha_BP_Laksh(l1=l1_hat, l2=l2_hat)$max_alpha
    alpha_hat <- ifelse(alpha_hat > max_alpha, max_alpha,
                        ifelse(alpha_hat < min_alpha, min_alpha, alpha_hat))
    return(alpha_hat)
  }

  # For alpha
  aux <- moments_estim_BP_Laksh(x_without_zeros)[3]
  alpha_hat <- as.numeric(aux)

  # Adjusting the alpha_hat's to ensure that it is inside in a correct interval
  alpha_hat <- aux_check(alpha_hat=alpha_hat, l1_hat=l1_hat, l2_hat=l2_hat)

  res <- c(l1_hat=l1_hat,
           l2_hat=l2_hat,
           alpha_hat=alpha_hat,
           psi_hat=psi_hat)

  return(round(res, digits=4))
}


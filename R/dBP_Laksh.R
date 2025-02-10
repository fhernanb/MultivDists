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
#' @param n number of random observations.
#' @param l1 mean for the marginal \eqn{X_1} variable with Poisson distribution.
#' @param l2 mean for the marginal \eqn{X_2} variable with Poisson distribution.
#' @param alpha third parameter.
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
#' @rdname dBP_Laksh
#' @importFrom stats runif
#' @export
#' @useDynLib MultivDists
#' @importFrom Rcpp sourceCpp
rBP_Laksh <- function(n, l1, l2, alpha,
                      max_val_x1=NULL, max_val_x2=NULL) {

  if(is.null(max_val_x1)) max_val_x1 <- 100
  if(is.null(max_val_x2)) max_val_x2 <- 100

  val_x1 <- 0:max_val_x1
  val_x2 <- 0:max_val_x2
  grid_x <- expand.grid(val_x1, val_x2)
  grid_x <- as.matrix(grid_x)
  probs <- dBP_Laksh(x=grid_x, l1=l1, l2=l2, alpha=alpha)
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

  positions <- pseudo_r_BP_Laksh(u=runif(n),
                                 cumul_probs=distri$cumul_probs)
  res <- distri[positions, 1:2]
  rownames(res) <- NULL
  colnames(res) <- c("X1", "X2")
  res <- as.matrix(res)
  return(res)
}
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
#' A interval with possible values for \eqn{\alpha}.
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
#'
#'
#'
#' Moment estimations for Bivariate Poisson Distribution - Laksh
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function obtains moment estimators for the Bivariate Poisson distribution
#' under the parameterization of Lakshminarayana et. al (1993).
#'
#' @param x vector or matrix of quantiles. When \code{x} is a matrix,
#' each row is taken to be a quantile and columns correspond to the number of dimensions \code{p}.
#'
#' @returns
#' Returns a vector with \eqn{\hat{\lambda_1}}, \eqn{\hat{\lambda_2}} and \eqn{\hat{\alpha}}.
#'
#' @references
#' Lakshminarayana, J., Pandit, S. N., & Srinivasa Rao, K. (1999). On a bivariate Poisson distribution. Communications in Statistics-Theory and Methods, 28(2), 267-276.
#'
#' @example examples/examples_moments_estim_BP_Laksh.R
#'
#' @importFrom stats cor
#' @export
moments_estim_BP_Laksh <- function(x) {

  # For l1
  l1_hat <- mean(x[,1])

  # For l1
  l2_hat <- mean(x[,2])

  n <- nrow(x)
  k <- 1 - exp(-1)
  A <- exp(-l1_hat*k)
  B <- exp(-l2_hat*k)

  # For alpha
  c_x_y <- sum((x[,1]-l1_hat) * (x[,2]-l2_hat)) / n
  alpha_hat <- c_x_y / (l1_hat * l2_hat * A * B * k^2)

  # Auxiliar functio to check if alpha_hat is inside the interval
  aux_check <- function(alpha_hat, l1_hat, l2_hat) {
    min_alpha <- correct_alpha_BP_Laksh(l1=l1_hat, l2=l2_hat)$min_alpha
    max_alpha <- correct_alpha_BP_Laksh(l1=l1_hat, l2=l2_hat)$max_alpha
    alpha_hat <- ifelse(alpha_hat > max_alpha, max_alpha,
                        ifelse(alpha_hat < min_alpha, min_alpha, alpha_hat))
    return(alpha_hat)
  }

  # Using the correlation to estimate alpha
  alpha_hat_cor <- cor(x)[2] / (sqrt(l1_hat*l2_hat)*k^2*exp(-(l1_hat+l2_hat)*k))

  # Adjusting the alpha_hat's to ensure that it is inside in a correct interval
  alpha_hat     <- aux_check(alpha_hat=alpha_hat,     l1_hat=l1_hat, l2_hat=l2_hat)
  alpha_hat_cor <- aux_check(alpha_hat=alpha_hat_cor, l1_hat=l1_hat, l2_hat=l2_hat)

  res <- c(l1_hat=l1_hat,
           l2_hat=l2_hat,
           #alpha_hat=alpha_hat,
           alpha_hat_cor=alpha_hat_cor
           )

  return(round(res, digits=4))
}
#'
#'
#'
llBP_Laksh <- function(param, x) {
  l1    <- param[1]  # param is the vector parameters
  l2    <- param[2]
  alpha <- param[3]
  sum(dBP_Laksh(x=x, l1=l1, l2=l2, alpha=alpha, log=TRUE))
}

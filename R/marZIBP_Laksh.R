#' Marginal Zero Inflated Bivariate Poisson Distribution - Laksh
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function fit a marginal Zero Inflated Bivariate Poisson
#' distribution under the of Lakshminarayana et. al (1993).
#'
#' @param mu1.fo a formula object to explain the marginal mean \eqn{\mu_1}, the first response \eqn{X_1} on the left of an ~ operator, and the terms, separated by operators, on the right.
#' @param mu2.fo a formula object to explain the marginal mean \eqn{\mu_2}, the second response \eqn{X_2} on the left of an ~ operator, and the terms, separated by operators, on the right.
#' @param psi.fo a formula object to explain the \eqn{\psi} parameter. An example could be: \code{~ age + salary}.
#' @param data a data frame containing the variables occurring in the formula2.
#' @param initial.values a vector with possible initial values for the parameters.
#'
#' @returns
#' Returns a object with class "marZIBPLaksh".
#'
#' @references
#' Lakshminarayana, J., Pandit, S. N., & Srinivasa Rao, K. (1999). On a bivariate Poisson distribution. Communications in Statistics-Theory and Methods, 28(2), 267-276.
#'
#' @example examples/examples_marZIBP_Laksh.R
#'
#' @export
marZIBP_Laksh <- function(mu1.fo, mu2.fo, psi.fo, data,
                          initial.values=NULL) {
  stopifnot (class(mu1.fo) == "formula")
  stopifnot (class(mu2.fo) == "formula")
  stopifnot (class(psi.fo) == "formula")
  matri <- model_matrix_marZIBP_Laksh(mu1.fo, mu2.fo, psi.fo, data)
  res <- fit.marZIBP_Laksh(matri, initial.values)
  class(res) <- "marZIBPLaksh"
  res
}
#' @importFrom stats model.matrix model.frame
model_matrix_marZIBP_Laksh <- function(mu1.fo, mu2.fo, psi.fo, data=NULL) {
  stopifnot (class(mu1.fo) == "formula")
  stopifnot (class(mu2.fo) == "formula")
  stopifnot (class(psi.fo) == "formula")
  W  <- model.matrix(object=psi.fo, data=data)
  X1 <- model.matrix(object=mu1.fo, data=data)
  X2 <- model.matrix(object=mu2.fo, data=data)
  y1 <- model.frame(mu1.fo, data=data)[, 1]
  y2 <- model.frame(mu2.fo, data=data)[, 1]
  list(X1=X1, X2=X2, W=W, y1=y1, y2=y2)
}
#' @importFrom stats optim
#' @importFrom numDeriv hessian
fit.marZIBP_Laksh <- function(matri, initial.values) {
  np_mu1 <- ncol(matri$X1)  # Number of parameters in mu1
  np_mu2 <- ncol(matri$X2)  # Number of parameters in mu2
  np_psi <- ncol(matri$W)   # Number of parameters in psi
  X1 <- matri$X1  # Model matrix for mu1
  X2 <- matri$X2  # Model matrix for mu2
  W  <- matri$W   # Model matrix for psi
  y1 <- matri$y1  # Response variable 1
  y2 <- matri$y2  # Response variable 2
  names_X1 <- colnames(matri$X1)
  names_X2 <- colnames(matri$X2)
  names_W  <- colnames(matri$W)

  # Initial value for alpha
  theta <- moments_estim_ZIBP_Laksh(cbind(y1, y2))
  min_alpha <- correct_alpha_BP_Laksh(l1=theta[1], l2=theta[2])$min_alpha
  max_alpha <- correct_alpha_BP_Laksh(l1=theta[1], l2=theta[2])$max_alpha

  # Interval to explore parameters
  lower_values <- c(rep(-Inf, ncol(X1)+ncol(X2)+ncol(W)), min_alpha)
  upper_values <- c(rep( Inf, ncol(X1)+ncol(X2)+ncol(W)), max_alpha)

  if (is.null(initial.values) | length(initial.values) != (np_mu1+np_mu2+np_psi+1))
    initial.values <- rep(0, np_mu1+np_mu2+np_psi+1)

  fit <- optim(par=initial.values,
               method="L-BFGS-B",
               fn=ll_marZIBP_Laksh,
               lower=lower_values,
               upper=upper_values,
               y1=y1, y2=y2,
               X1=X1, X2=X2,
               W=W,
               control=list(trace=3, maxit=10000))
  fit$objective <- -fit$value

  # Obtaining fitted mu1, mu2, alpha and mu
  betas_X1 <- matrix(fit$par[1:ncol(X1)], ncol=1)
  betas_X2 <- matrix(fit$par[(ncol(X1)+1) : (ncol(X1)+ncol(X2))], ncol=1)
  betas_W  <- matrix(fit$par[(ncol(X1)+ncol(X2)+1) : (ncol(X1)+ncol(X2)+ncol(W))], ncol=1)

  # Auxiliar function
  logit_inv <- function(x) exp(x) / (1+exp(x))

  fit$fitted.mu1   <- exp(X1 %*% betas_X1)
  fit$fitted.mu2   <- exp(X2 %*% betas_X2)
  fit$fitted.psi   <- logit_inv(W %*% betas_W)
  fit$fitted.alpha <- fit$par[ncol(X1)+ncol(X2)+ncol(W)+1]

  fit$fitted.l1    <- fit$fitted.mu1 / (1-fit$fitted.psi)
  fit$fitted.l2    <- fit$fitted.mu2 / (1-fit$fitted.psi)

  # Obtaining the hessian
  fit$Hessian <- numDeriv::hessian(func=ll_marZIBP_Laksh,
                                   x=fit$par,
                                   method='Richardson',
                                   y1=y1, y2=y2,
                                   X1=X1, X2=X2,
                                   W=W)

  inputs <- list(y1=y1, y2=y2, X1=X1, X2=X2, W=W,
                 np_mu1=np_mu1, np_mu2=np_mu2, np_psi=np_psi)
  fit <- c(fit, inputs)
}
#'
ll_marZIBP_Laksh <- function(theta, y1, y2, X1, X2, W) {
  betas_X1 <- matrix(theta[1:ncol(X1)], ncol=1)
  betas_X2 <- matrix(theta[(ncol(X1)+1) : (ncol(X1)+ncol(X2))], ncol=1)
  betas_W  <- matrix(theta[(ncol(X1)+ncol(X2)+1) : (ncol(X1)+ncol(X2)+ncol(W))], ncol=1)
  alpha    <- theta[ncol(X1)+ncol(X2)+ncol(W)+1]
  # Now we obtain the mean for Y1 and Y2
  mu1 <- exp(X1 %*% betas_X1)
  mu2 <- exp(X2 %*% betas_X2)
  # To obtain the parameter psi
  logit_inv <- function(x) exp(x) / (1+exp(x))
  psi <- logit_inv(W %*% betas_W)
  # Now we obtain lambda1 and lambda2
  l1 <- mu1 / (1-psi)
  l2 <- mu2 / (1-psi)
  y <- as.matrix(cbind(y1, y2))
  ll <- sum(dZIBP_Laksh(x=y, l1=l1, l2=l2, alpha=alpha, psi=psi, log=TRUE))
  return(-ll)  # minus to use with optim/nlminb function
}
#'
#' Summary table for marZIP - Laksh
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function obtains the summary table for objects of class marZIBPLaksh.
#'
#' @param object of class marZIBPLaksh.
#' @param ... aditional arguments.
#'
#' @returns
#' Returns the summary table.
#'
#' @example examples/examples_marZIBP_Laksh.R
#'
#' @export
#' @importFrom stats pnorm printCoefmat
summary.marZIBPLaksh <- function(object, ...) {
  .myenv <- environment()
  var.list <- as.list(object)
  list2env(var.list , envir = .myenv)
  estimate <- object$par
  elements <- sqrt(diag(solve(object$Hessian))) # diagonal of Hessian^-1
  se <- elements
  zvalue   <- estimate / se
  pvalue   <- 2 * pnorm(abs(zvalue), lower.tail=F)
  res      <- cbind(estimate=estimate, se=se, zvalue=zvalue, pvalue=pvalue)
  res      <- data.frame(res)
  colnames(res) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')
  np_mu1 <- object$np_mu1
  np_mu2 <- object$np_mu2
  np_psi <- object$np_psi
  res.mu1 <- res[1:np_mu1,]
  res.mu2 <- res[(np_mu1+1):(np_mu1+np_mu2),]
  res.psi  <- res[(np_mu1+np_mu2+1):(np_mu1+np_mu2+np_psi),]
  res.alpha <- res[np_mu1+np_mu2+np_psi+1,]
  rownames(res.mu1)   <- colnames(object$X1)
  rownames(res.mu2)   <- colnames(object$X2)
  rownames(res.psi)   <- colnames(object$W)
  rownames(res.alpha) <- "(Intercept)"
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for log(mu1) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res.mu1, P.values=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for log(mu2) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res.mu2, P.values=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for logit(psi) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res.psi, P.values=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Estimation for alpha \n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res.alpha, P.values=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
}

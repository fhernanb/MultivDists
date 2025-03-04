#' Zero Inflated Bivariate Poisson model - Geoffroy
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function fit a Zero Inflated Bivariate Poisson model
#' under Geoffroy et. al (2021).
#'
#' @param l1.fo a formula object to explain the parameter \eqn{\lambda_1}, the first response \eqn{X_1} on the left of an ~ operator, and the terms, separated by operators, on the right.
#' @param l2.fo a formula object to explain the parameter \eqn{\lambda_2}, the second response \eqn{X_2} on the left of an ~ operator, and the terms, separated by operators, on the right.
#' @param psi.fo a formula object to explain the \eqn{\psi} parameter. An example could be: \code{~ age + salary}.
#' @param data a data frame containing the variables occurring in the formula.
#' @param initial.values a vector with possible initial values for the parameters.
#' @param optimizer the optimizer to be used: nlminb, optim.
#'
#' @returns
#' Returns a object with class "ZIBPGeoffroy".
#'
#' @references
#' Kouakou, K. J. G., Hili, O., & Dupuy, J. F. (2021). Estimation in the zero-inflated bivariate Poisson model with an application to health-care utilization data. Afrika Statistika, 16(2), 2767-2788.
#'
#' @example examples/examples_ZIBP_Geoffroy.R
#'
#' @export
ZIBP_Geoffroy <- function(l1.fo, l2.fo, psi.fo, data,
                          initial.values=NULL,
                          optimizer="nlminb") {
  stopifnot (class(l1.fo) == "formula")
  stopifnot (class(l2.fo) == "formula")
  stopifnot (class(psi.fo) == "formula")
  matri <- model.matrix.ZIBP_Geoffroy(l1.fo, l2.fo, psi.fo, data)
  res <- fit.ZIBP_Geoffroy(matri, initial.values, optimizer)
  class(res) <- "ZIBPGeoffroy"
  res
}
#' @importFrom stats model.matrix model.frame
model.matrix.ZIBP_Geoffroy <- function(l1.fo, l2.fo, psi.fo, data=NULL) {
  stopifnot (class(l1.fo) == "formula")
  stopifnot (class(l2.fo) == "formula")
  stopifnot (class(psi.fo) == "formula")
  W  <- model.matrix(object=psi.fo, data=data)
  X1 <- model.matrix(object=l1.fo, data=data)
  X2 <- model.matrix(object=l2.fo, data=data)
  y1 <- model.frame(l1.fo, data=data)[, 1]
  y2 <- model.frame(l2.fo, data=data)[, 1]
  list(X1=X1, X2=X2, W=W, y1=y1, y2=y2)
}
#' @importFrom stats optim nlminb
#' @importFrom numDeriv hessian
fit.ZIBP_Geoffroy <- function(matri, initial.values, optimizer) {
  np_l1  <- ncol(matri$X1)  # Number of parameters in l1
  np_l2  <- ncol(matri$X2)  # Number of parameters in l2
  np_psi <- ncol(matri$W)   # Number of parameters in psi
  X1 <- matri$X1  # Model matrix for l1
  X2 <- matri$X2  # Model matrix for l2
  W  <- matri$W   # Model matrix for psi
  y1 <- matri$y1  # Response variable 1
  y2 <- matri$y2  # Response variable 2
  names_X1 <- colnames(matri$X1)
  names_X2 <- colnames(matri$X2)
  names_W  <- colnames(matri$W)

  # Initial value for l0
  theta <- moments_estim_ZIBP_Geoffroy(cbind(y1, y2))
  min_l0 <- 0.001
  max_l0 <- Inf

  # Interval to explore parameters
  lower_values <- c(rep(-Inf, ncol(X1)+ncol(X2)+ncol(W)), min_l0)
  upper_values <- c(rep( Inf, ncol(X1)+ncol(X2)+ncol(W)), max_l0)

  if (is.null(initial.values) |
      length(initial.values) != (np_l1+np_l2+np_psi+1))
    initial.values <- c(rep(0, np_l1+np_l2+np_psi), 1)

  if (np_l1+np_l2+np_psi+1 == 4)
    initial.values <- theta

  if (optimizer=="nlminb") {
    fit <- nlminb(start=initial.values,
                  objective=ll_ZIBP_Geoffroy,
                  lower=lower_values,
                  upper=upper_values,
                  y1=y1, y2=y2,
                  X1=X1, X2=X2,
                  W=W)
    fit$objective <- -fit$objective
  }

  else {
    fit <- optim(par=initial.values,
                 method="Nelder-Mead",
                 fn=ll_ZIBP_Geoffroy,
                 lower=lower_values,
                 upper=upper_values,
                 y1=y1, y2=y2,
                 X1=X1, X2=X2,
                 W=W)
    fit$objective <- -fit$value
  }

  # Obtaining fitted l1, l2, l0 and psi
  betas_X1 <- matrix(fit$par[1:ncol(X1)], ncol=1)
  betas_X2 <- matrix(fit$par[(ncol(X1)+1) : (ncol(X1)+ncol(X2))], ncol=1)
  betas_W  <- matrix(fit$par[(ncol(X1)+ncol(X2)+1) : (ncol(X1)+ncol(X2)+ncol(W))], ncol=1)

  # Auxiliar function
  logit_inv <- function(x) exp(x) / (1+exp(x))

  fit$fitted.l1  <- exp(X1 %*% betas_X1)
  fit$fitted.l2  <- exp(X2 %*% betas_X2)
  fit$fitted.psi <- logit_inv(W %*% betas_W)
  fit$fitted.l0  <- fit$par[ncol(X1)+ncol(X2)+ncol(W)+1]

  # Obtaining the hessian
  fit$Hessian <- numDeriv::hessian(func=ll_ZIBP_Geoffroy,
                                   x=fit$par,
                                   method='Richardson',
                                   y1=y1, y2=y2,
                                   X1=X1, X2=X2,
                                   W=W)

  inputs <- list(y1=y1, y2=y2, X1=X1, X2=X2, W=W,
                 np_l1=np_l1, np_l2=np_l2, np_psi=np_psi)
  fit <- c(fit, inputs)
}
#'
ll_ZIBP_Geoffroy <- function(theta, y1, y2, X1, X2, W) {
  betas_X1 <- matrix(theta[1:ncol(X1)], ncol=1)
  betas_X2 <- matrix(theta[(ncol(X1)+1) : (ncol(X1)+ncol(X2))], ncol=1)
  betas_W  <- matrix(theta[(ncol(X1)+ncol(X2)+1) : (ncol(X1)+ncol(X2)+ncol(W))], ncol=1)
  l0       <- theta[ncol(X1)+ncol(X2)+ncol(W)+1]
  # Now we obtain the parameter l1 and l2 for Y1 and Y2
  l1 <- exp(X1 %*% betas_X1)
  l2 <- exp(X2 %*% betas_X2)
  # To obtain the parameter psi
  logit_inv <- function(x) exp(x) / (1+exp(x))
  psi <- logit_inv(W %*% betas_W)
  y <- as.matrix(cbind(y1, y2))
  ll <- sum(dZIBP_Geoffroy(x=y, l1=l1, l2=l2, l0=l0, psi=psi, log=TRUE))

  # resul <- c(betas_X1, betas_X2, betas_W, l0)
  # resul <- round(resul, digits=2)
  # names(resul) <- NULL
  # print(resul)

  return(-ll)  # minus to use with optim/nlminb function
}
#'
#' Summary table for marZIP - Geoffroy
#'
#' @author Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
#'
#' @description
#' This function obtains the summary table for objects of class marZIBPGeoffroy.
#'
#' @param object of class marZIBPGeoffroy.
#' @param ... aditional arguments.
#'
#' @returns
#' Returns the summary table.
#'
#' @example examples/examples_marZIBP_Geoffroy.R
#'
#' @export
#' @importFrom stats pnorm printCoefmat
summary.ZIBPGeoffroy <- function(object, ...) {
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
  np_l1 <- object$np_l1
  np_l2 <- object$np_l2
  np_psi <- object$np_psi
  res.l1 <- res[1:np_l1,]
  res.l2 <- res[(np_l1+1):(np_l1+np_l2),]
  res.psi  <- res[(np_l1+np_l2+1):(np_l1+np_l2+np_psi),]
  res.l0 <- res[np_l1+np_l2+np_psi+1,]
  rownames(res.l1)   <- colnames(object$X1)
  rownames(res.l2)   <- colnames(object$X2)
  rownames(res.psi)   <- colnames(object$W)
  rownames(res.l0) <- "(Intercept)"
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for log(l1) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res.l1, P.values=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for log(l2) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res.l2, P.values=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for logit(psi) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res.psi, P.values=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Estimation for l0 \n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res.l0, P.values=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
}
#' @rdname summary.marZIBPGeoffroy
print.ZIBPGeoffroy <- function(object, ...) {
  print("Please use the summary method to obtain complete information.")
}

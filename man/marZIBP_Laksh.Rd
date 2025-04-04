% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marZIBP_Laksh.R
\name{marZIBP_Laksh}
\alias{marZIBP_Laksh}
\title{Marginal Zero Inflated Bivariate Poisson Distribution - Laksh}
\usage{
marZIBP_Laksh(mu1.fo, mu2.fo, psi.fo, data, initial.values = NULL)
}
\arguments{
\item{mu1.fo}{a formula object to explain the marginal mean \eqn{\mu_1}, the first response \eqn{X_1} on the left of an ~ operator, and the terms, separated by operators, on the right.}

\item{mu2.fo}{a formula object to explain the marginal mean \eqn{\mu_2}, the second response \eqn{X_2} on the left of an ~ operator, and the terms, separated by operators, on the right.}

\item{psi.fo}{a formula object to explain the \eqn{\psi} parameter. An example could be: \code{~ age + salary}.}

\item{data}{a data frame containing the variables occurring in the formula2.}

\item{initial.values}{a vector with possible initial values for the parameters.}
}
\value{
Returns a object with class "marZIBPLaksh".
}
\description{
This function fit a marginal Zero Inflated Bivariate Poisson
distribution under the of Lakshminarayana et. al (1993).
}
\examples{
# Example 1 ---------------------------------------------------------------
l1 <- 3
l2 <- 4
alpha <- -0.90
psi <- 0.2

set.seed(12345678)
data1 <- rZIBP_Laksh(n=100, l1=l1, l2=l2, alpha=alpha, psi=psi)
data1 <- as.data.frame(data1)

# To fit the model
mod1 <- NULL
mod1 <- marZIBP_Laksh(mu1.fo=X1~1,
                      mu2.fo=X2~1,
                      psi.fo=~1,
                      data=data1)

# To obtain the usual summary table
summary(mod1)

# To explore the estimations of l1, l2, mu and p

# To obtain E(Y1)=v1 and E(Y2)=v2
mod1$fitted.mu1[1]
mod1$fitted.mu2[1]

# To compare sample means with v1 and v2
colMeans(data1)

# To obtain alpha and psi
mod1$fitted.alpha
mod1$fitted.psi[1]

# To obtain l1 and l2
mod1$fitted.l1[1]
mod1$fitted.l2[1]

# Example 2 ---------------------------------------------------------------

gen_data_ZIBP_Laksh <- function(n=100) {
  # To generate the covariates
  x1 <- runif(n=n)
  x2 <- runif(n=n)

  # To generate the means
  mu1 <- exp(-2 + 3.5 * x1 + 2.7 * x2)
  mu2 <- exp(-1 + 1.3 * x1 + 2.1 * x2)

  # To generate the psi
  logit_inv <- function(x) exp(x) / (1+exp(x))
  psi  <- logit_inv(-2.4 + 1.2 * x2)

  alpha <- -1

  # To obtain lambdas
  l1 <- mu1 / (1-psi)
  l2 <- mu2 / (1-psi)

  # To generate Y1 and Y2
  y <- NULL
  for (i in 1:n)
    y <- rbind(y, rZIBP_Laksh(n=1, l1=l1[i], l2=l2[i],
                              alpha=alpha, psi=psi[i]))

  # To create the dataset
  dataset <- data.frame(y1=y[,1], y2=y[,2],
                        x1=x1, x2=x2,
                        mu1=mu1, mu2=mu2,
                        alpha=alpha, psi=psi,
                        l1=l1, l2=l2)

  return(dataset)
}

set.seed(123456)
data2 <- gen_data_ZIBP_Laksh(n=100)
head(data2, n=8)

mod2 <- NULL
mod2 <- marZIBP_Laksh(mu1.fo=y1~x1+x2,
                      mu2.fo=y2~x1+x2,
                      psi.fo=~x2,
                      data=data2)

summary(mod2)

mod2$fitted.alpha

}
\references{
Lakshminarayana, J., Pandit, S. N., & Srinivasa Rao, K. (1999). On a bivariate Poisson distribution. Communications in Statistics-Theory and Methods, 28(2), 267-276.
}
\author{
Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ZIBP_Geoffroy.R
\name{ZIBP_Geoffroy}
\alias{ZIBP_Geoffroy}
\title{Zero Inflated Bivariate Poisson model - Geoffroy}
\usage{
ZIBP_Geoffroy(
  l1.fo,
  l2.fo,
  psi.fo,
  data,
  initial.values = NULL,
  optimizer = "nlminb"
)
}
\arguments{
\item{l1.fo}{a formula object to explain the parameter \eqn{\lambda_1}, the first response \eqn{X_1} on the left of an ~ operator, and the terms, separated by operators, on the right.}

\item{l2.fo}{a formula object to explain the parameter \eqn{\lambda_2}, the second response \eqn{X_2} on the left of an ~ operator, and the terms, separated by operators, on the right.}

\item{psi.fo}{a formula object to explain the \eqn{\psi} parameter. An example could be: \code{~ age + salary}.}

\item{data}{a data frame containing the variables occurring in the formula.}

\item{initial.values}{a vector with possible initial values for the parameters.}

\item{optimizer}{the optimizer to be used: nlminb, optim.}
}
\value{
Returns a object with class "ZIBPGeoffroy".
}
\description{
This function fit a Zero Inflated Bivariate Poisson model
under Geoffroy et. al (2021).
}
\examples{
# Example 1 ---------------------------------------------------------------
# Estimating parameters as a regression model

l1 <- 5
l2 <- 3
psi <- 0.20
l0 <- 1.20

set.seed(12345)
data1 <- rZIBP_Geoffroy(n=500, l1=l1, l2=l2, l0=l0, psi=psi)
data1 <- as.data.frame(data1)
colnames(data1) <- c("y1", "y2")

mod <- ZIBP_Geoffroy(l1.fo=y1~1,
                     l2.fo=y2~1,
                     psi.fo=~1,
                     data=data1)

summary(mod)

c(exp(mod$par[1:3]), mod$par[4])

moments_estim_ZIBP_Geoffroy(data1)


# Example 2 ---------------------------------------------------------------
# Replicating application from section 5.1 from Geoffroy et. al (2021).
\dontrun{
  library(AER)
  data("NMES1988")

  # To obtain some statistics as in the paper
  y1 <- NMES1988$nvisits
  y2 <- NMES1988$novisits

  cor(y1, y2)
  cor.test(y1, y2)
  cov(y1, y2)
  mean(y1==0 & y2==0)

  # Model as in section 5.1 from Geoffroy et. al (2021).
  mod3 <- NULL
  mod3 <- ZIBP_Geoffroy(l1.fo=nvisits ~health+chronic+age+gender+
                          married+school+income+medicaid,
                        l2.fo=novisits~health+chronic+age+gender+
                          married+school+income+medicaid,
                        psi.fo=~chronic+gender+school+medicaid,
                        data=NMES1988)

  summary(mod3)
}
}
\references{
Kouakou, K. J. G., Hili, O., & Dupuy, J. F. (2021). Estimation in the zero-inflated bivariate Poisson model with an application to health-care utilization data. Afrika Statistika, 16(2), 2767-2788.
}
\author{
Freddy Hernandez-Barajas, \email{fhernanb@unal.edu.co}
}

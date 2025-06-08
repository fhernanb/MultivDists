
# Example 1 ---------------------------------------------------------------
# Probability for single values of X1 and X2
dZIBP_HOL(c(0, 0), l1=3, l2=4, l0=1, psi=0.15)
dZIBP_HOL(c(1, 0), l1=3, l2=4, l0=1, psi=0.15)
dZIBP_HOL(c(0, 1), l1=3, l2=4, l0=1, psi=0.15)

# Probability for a matrix the values of X1 and X2
x <- matrix(c(0, 0,
              1, 0,
              0, 1), ncol=2, byrow=TRUE)
x
dZIBP_HOL(x=x, l1=3, l2=4, l0=1, psi=0.15)

# Checking if the probabilities sum 1
val_x1 <- val_x2 <- 0:50
space <- expand.grid(val_x1, val_x2)
space <- as.matrix(space)

l1 <- 3
l2 <- 4
l0 <- 1.27
psi <- 0.15

probs <- dZIBP_HOL(x=space, l1=l1, l2=l2, l0=l0, psi=psi)
sum(probs)

# Example 2 ---------------------------------------------------------------
# Heat map for a ZIBP_HOL

l1 <- 3
l2 <- 4
l0 <- 1.2
psi <- 0.03

X1 <- 0:10
X2 <- 0:10
data <- expand.grid(X1=X1, X2=X2)
data$Prob <- dZIBP_HOL(x=data, l1=l1, l2=l2, l0=l0, psi=psi)
data$X1 <- factor(data$X1)
data$X2 <- factor(data$X2)

library(ggplot2)
ggplot(data, aes(X1, X2, fill=Prob)) +
  geom_tile() +
  scale_fill_gradient(low="darkgreen", high="yellow")


# Example 3 ---------------------------------------------------------------
# Generating random values and moment estimations

l1 <- 13
l2 <- 8
l0 <- 5
psi <- 0.15

x <- rZIBP_HOL(n=500, l1, l2, l0, psi)
moments_estim_ZIBP_HOL(x)

# Example 4 ---------------------------------------------------------------
# Estimating the parameters using the loglik function

# Loglik function
llZIBP_HOL <- function(param, x) {
  l1 <- param[1]  # param: is the parameter vector
  l2 <- param[2]
  l0 <- param[3]
  psi <- param[4]
  sum(dZIBP_HOL(x=x, l1=l1, l2=l2,
                     l0=l0, psi=psi, log=TRUE))
}

# The known parameters
l1 <- 5
l2 <- 3
psi <- 0.20
l0 <- 1.20

set.seed(12345)
x <- rZIBP_HOL(n=500, l1=l1, l2=l2, l0=l0, psi=psi)

# To obtain reasonable values for l0
start_param <- moments_estim_ZIBP_HOL(x)
start_param

# Estimating parameters
res1 <- optim(fn = llZIBP_HOL,
              par = start_param,
              lower = c(0.001, 0.001, 0.001, 0.0001),
              upper = c(  Inf,   Inf,   Inf, 0.9999),
              method = "L-BFGS-B",
              control = list(maxit=100000, fnscale=-1),
              x=x)

res1


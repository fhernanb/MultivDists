
# Example 1 ---------------------------------------------------------------
# Probability for single values of X1 and X2
dZIBP_Laksh(c(0, 0), l1=3, l2=4, alpha=1, psi=0.15)
dZIBP_Laksh(c(1, 0), l1=3, l2=4, alpha=1, psi=0.15)
dZIBP_Laksh(c(0, 1), l1=3, l2=4, alpha=1, psi=0.15)

# Probability for a matrix the values of X1 and X2
x <- matrix(c(0, 0,
              1, 0,
              0, 1), ncol=2, byrow=TRUE)
x
dZIBP_Laksh(x=x, l1=3, l2=4, alpha=1, psi=0.15)

# Checking if the probabilities sum 1
val_x1 <- val_x2 <- 0:50
space <- expand.grid(val_x1, val_x2)
space <- as.matrix(space)

l1 <- 3
l2 <- 4
alpha <- -1.27
psi <- 0.15

probs <- dZIBP_Laksh(x=space, l1=l1, l2=l2, alpha=alpha, psi=psi)
sum(probs)

# Example 2 ---------------------------------------------------------------
# Heat map for a ZIBP_Laksh

l1 <- 3
l2 <- 4
alpha <- -1.2
psi <- 0.15

X1 <- 0:10
X2 <- 0:10
data <- expand.grid(X1=X1, X2=X2)
data$Prob <- dZIBP_Laksh(x=data, l1=l1, l2=l2, alpha=alpha, psi=psi)
data$X1 <- factor(data$X1)
data$X2 <- factor(data$X2)

library(ggplot2)
ggplot(data, aes(X1, X2, fill=Prob)) +
  geom_tile() +
  scale_fill_gradient(low="darkgreen", high="yellow")


# Example 3 ---------------------------------------------------------------
# Generating random values and moment estimations

l1 <- 15
l2 <- 13
correct_alpha_BP_Laksh(l1, l2)
alpha <- 0.9
psi <- 0.20

x <- rZIBP_Laksh(n=50000, l1, l2, alpha, psi)
moments_estim_ZIBP_Laksh(x)

# Example 4 ---------------------------------------------------------------
# Estimating the parameters using the loglik function

# Loglik function
llZIBP_Laksh <- function(param, x) {
  l1    <- param[1]  # param: is the parameter vector
  l2    <- param[2]
  alpha <- param[3]
  psi   <- param[4]
  sum(dZIBP_Laksh(x=x, l1=l1, l2=l2,
                  alpha=alpha, psi=psi, log=TRUE))
}

# The known parameters
l1 <- 5
l2 <- 3
correct_alpha_BP_Laksh(l1=l1, l2=l2)
alpha <- -1.20
psi <- 0.20

set.seed(12345)
x <- rZIBP_Laksh(n=500, l1=l1, l2=l2, alpha=alpha, psi=psi)

# To obtain reasonable values for alpha
theta <- as.numeric(moments_estim_ZIBP_Laksh(x))
theta

# To create start parameters
min_alpha <- correct_alpha_BP_Laksh(l1=theta[1],
                                    l2=theta[2])$min_alpha
max_alpha <- correct_alpha_BP_Laksh(l1=theta[1],
                                    l2=theta[2])$max_alpha

start_param <- theta
names(start_param) <- c("l1_hat", "l2_hat", "alpha_hat", "psi_hat")
start_param

# Estimating parameters
res1 <- optim(fn = llZIBP_Laksh,
              par = start_param,
              lower = c(0.001, 0.001, min_alpha, 0.0001),
              upper = c(  Inf,   Inf, max_alpha, 0.9999),
              method = "L-BFGS-B",
              control = list(maxit=100000, fnscale=-1),
              x=x)

res1




# Example 1 ---------------------------------------------------------------
# Probability for single values of X1 and X2
dBP_Geoffroy(c(0, 0), l1=3, l2=4, mu=1)
dBP_Geoffroy(c(1, 0), l1=3, l2=4, mu=1)
dBP_Geoffroy(c(0, 1), l1=3, l2=4, mu=1)

# Probability for a matrix the values of X1 and X2
x <- matrix(c(0, 0,
              1, 0,
              0, 1), ncol=2, byrow=TRUE)
x
dBP_Geoffroy(x=x, l1=3, l2=4, mu=1)

# Checking if the probabilities sum 1
val_x1 <- val_x2 <- 0:50
space <- expand.grid(val_x1, val_x2)
space <- as.matrix(space)

l1 <- 3
l2 <- 4
mu <- 5

probs <- dBP_Geoffroy(x=space, l1=l1, l2=l2, mu=mu)
sum(probs)

# Example 2 ---------------------------------------------------------------
# Heat map for a BP_Geoffroy

l1 <- 1
l2 <- 2
mu <- 4

X1 <- 0:10
X2 <- 0:10
data <- expand.grid(X1=X1, X2=X2)
data$Prob <- dBP_Geoffroy(x=data, l1=l1, l2=l2, mu=mu)
data$X1 <- factor(data$X1)
data$X2 <- factor(data$X2)

library(ggplot2)
ggplot(data, aes(X1, X2, fill=Prob)) +
  geom_tile() +
  scale_fill_gradient(low="darkgreen", high="pink")

# Example 3 ---------------------------------------------------------------
# Generating random values and moment estimations

l1 <- 1
l2 <- 2
mu <- 4

x <- rBP_Geoffroy(n=500, l1, l2, mu)
moments_estim_BP_Geoffroy(x)

# Example 4 ---------------------------------------------------------------
# Estimating the parameters using the loglik function

# Loglik function
llBP_Geoffroy <- function(param, x) {
  l1 <- param[1]  # param: is the parameter vector
  l2 <- param[2]
  mu <- param[3]
  sum(dBP_Geoffroy(x=x, l1=l1, l2=l2, mu=mu, log=TRUE))
}

# The known parameters
l1 <- 1
l2 <- 2
mu <- 4

set.seed(12345)
x <- rBP_Geoffroy(n=500, l1=l1, l2=l2, mu=mu)

# To obtain reasonable values for mu
theta <- as.numeric(moments_estim_BP_Geoffroy(x))
theta

start_param <- theta
names(start_param) <- c("l1_hat", "l2_hat", "mu_hat")

# Estimating parameters
res1 <- optim(fn = llBP_Geoffroy,
              par = start_param,
              lower = c(0.001, 0.001, 0.001),
              upper = c(  Inf,   Inf,   Inf),
              method = "L-BFGS-B",
              control = list(maxit=100000, fnscale=-1),
              x=x)

res1


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
  scale_fill_gradient(low="darkgreen", high="white")

# Generating random values and moment estimations
l1 <- 1
l2 <- 2
alpha <- -2.7
psi <- 0.30

x <- rZIBP_Laksh(n=500, l1, l2, alpha, psi)
moments_estim_ZIBP_Laksh(x)

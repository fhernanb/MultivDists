# Probability for single values of X1 and X2
dBP_Laksh(c(0, 0), l1=3, l2=4, alpha=1)
dBP_Laksh(c(1, 0), l1=3, l2=4, alpha=1)
dBP_Laksh(c(0, 1), l1=3, l2=4, alpha=1)

# Probability for a matrix the values of X1 and X2
x <- matrix(c(0, 0,
              1, 0,
              0, 1), ncol=2, byrow=TRUE)
x
dBP_Laksh(x=x, l1=3, l2=4, alpha=1)

# Checking if the probabilites sum 1
val_x1 <- val_x2 <- 0:50
space <- expand.grid(val_x1, val_x2)
space <- as.matrix(space)

l1 <- 3
l2 <- 4
alpha <- -1.27

probs <- dBP_Laksh(x=space, l1=l1, l2=l2, alpha=alpha)
sum(probs)

# Heatmap for a BP_Laksh

l1 <- 1
l2 <- 2
correct_alpha_BP_Laksh(l1=l1, l2=l2)
alpha <- -2.9

X1 <- 0:10
X2 <- 0:10
data <- expand.grid(X1=X1, X2=X2)
data$Prob <- dBP_Laksh(x=data, l1=l1, l2=l2, alpha=alpha)
data$X1 <- factor(data$X1)
data$X2 <- factor(data$X2)

library(ggplot2)
ggplot(data, aes(X1, X2, fill=Prob)) +
  geom_tile() +
  scale_fill_gradient(low="darkgreen", high="white")

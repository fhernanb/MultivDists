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


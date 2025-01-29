# Generating random values and moment estimations
l1 <- 1
l2 <- 2
alpha <- -2.7
x <- rBP_Laksh(n=500, l1, l2, alpha)

moments_estim_BP_Laksh(x)

# Analizing example from Famoye (2010)

freq <- c(34, 20, 4, 6, 4,
          17, 7, 0, 0, 0,
          6, 4, 1, 0, 0,
          0, 4, 0, 0, 0,
          0, 0, 0, 0, 0,
          2, 0, 0, 0, 0)

data_table <- matrix(freq, ncol=5, byrow=TRUE)
rownames(data_table) <- 0:5
colnames(data_table) <- 0:4

data_table

long <- as.data.frame.table(data_table)
x <- long[rep(1:nrow(long), long$Freq), -3]
x <- data.matrix(x)

colMeans(x)
var(x)
cor(x)
moments_estim_BP_Laksh(x)

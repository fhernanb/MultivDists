# Generating random values and moment estimations
l1 <- 1
l2 <- 2
correct_alpha_BP_Laksh(l1=l1, l2=l2)
alpha <- -2.7

x <- rBP_Laksh(n=500, l1, l2, alpha)
moments_estim_BP_Laksh(x)

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

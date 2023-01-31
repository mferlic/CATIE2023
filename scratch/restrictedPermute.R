N <- 1e4

A1 <- 2*rbinom(N,1,0.5) - 1

U <- rnorm(N)

p.R <- 0.5 + 0.2*A1 + 0.7*U
p.R[p.R > 1] <- 1
p.R[p.R < 0 ] <- 0
R <- rbinom(N, 1, p.R)
R.model <- lm(R ~ A1)
R.resid <- R.model$residuals
R.c <- R - mean(R)

# Mean Centered; TRUE DGP
Y <- 0.2*A1 + (0.5 - 0.1*A1)*R.resid + 0.8*U
lm(Y ~ A1) %>% summary()

lm(Y ~ A1*R.resid) %>% summary()

# Naive
lm(Y ~ A1*R.c) %>% summary()

X <- data.frame(A1, R.c, Y)

M <- 1e3
coefs <- matrix(0, M, 4)
permR <- R.c
for(i in 1:M){
  u <- runif(1)
  strata <- Y > quantile(Y, u)
  s1 <- sum(strata)
  s0 <- N - s1
  ps1 <- sample(s1)
  ps0 <- sample(s0)
  permR[strata] <- R.c[strata][ps1]
  permR[!strata] <- R.c[!strata][ps0]
  fit <- lm(Y ~ A1*permR)
  coefs[i, ] <- coefficients(fit)
}
colMeans(coefs)


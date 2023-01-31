# Test generation of Y_t data
# Y_t is Markov
library(dplyr)

# Baseline covariates
N <- 1000
X <- MASS::mvrnorm(N, c(0,0), matrix(c(2, 1, 1, 3), 2, 2))
X <- cbind(X, rbinom(N, 1, 0.5))

# Baseline
eta0 <- c(2, -1, 3) # random
Y0 <- X %*% eta0 + rnorm(N)
hist(Y0)

# First-stage
A1 <- rbinom(N, 1, 0.5)

# Response
eta1 <- c(0.2, -0.2, 0.3, 0.1, 0.2) # log-odds scale
pR <- plogis(cbind(X, Y0, A1) %*% eta1)
hist(pR)
mean(pR)
R <- rbinom(N, 1, pR)

# Adherence
eta2 <- c(-0.3, 0.1, -0.3, -0.3, -0.2, -0.4) # logit scale
lS1 <- cbind(X, Y0, A1, R) %*% eta2 # log-odds
pS1 <- plogis(lS1) #convert to probabilities
S1 <- rbinom(N, 1, pS1) # sample
hist(pS1)
mean(pS1)

# Time to Non-response
times = c(2,4,8)
eta3 <- c(intercept = -1, -0.2, 0.1, 0.1, -0.3, -0.2, t = -0.2) # logit scale
# Make person period dataframe
M <- cbind(1, X, Y0, A1)

lhaz <- cbind(1, X, Y0, A1) %*% eta3 #log-hazards, time-invariant
haz <- plogis(lhaz)
hist(haz)
et <- eventTime(c(2,4,8), haz)
hist(et)


eventTime <- function(times, hazard) {
  # Generate conditional draws for probability of failure in time periods
  # This assumes constant hazard!
  # Would need to loop person-period dataset
  n <- length(hazard)
  ET <- rep_len(0, n) # zero indicates no event-time
  for (i in times) {
    d <- runif(n) < hazard
    ET[ET == 0 & d == 1] = i # event in interval and survived until interval
  }
  return(ET)
}

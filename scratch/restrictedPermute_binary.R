N <- 1e3

A1 <- 2*rbinom(N,1,0.5) - 1
U <- rbinom(N, 1, 0.5)

Rbetas <- list(a1 = 0.2 )

# LPM
p.R <- 0.5 + Rbetas$a1*A1 - 0.1*U
R <- rbinom(N, 1, p.R)

R.model <- lm(R ~ A1)
ER_A1 <- 0.5 + Rbetas$a1*A1
R.resid <- R - ER_A1
R.c <- R - mean(R)

# Mean Centered; TRUE DGP
p.Y <- 0.5 + 0.2*A1 + (0.4 + 0.1*A1)*R.resid - 0.2*U
Y <- rbinom(N, 1, p.Y)

lm(Y ~ A1) %>% summary()

lm(Y ~ A1*R.resid) %>% summary()

# Naive
lm(Y ~ A1*R.c) %>% summary()


## Stratify on A1, permute R -----------------------------------------------------
# cov(R*, Y) = 0; betaA1 correct
M <- 1e3
coefs <- matrix(0, M, 4)
mcov <- array(0, c(3,3,M))

Rperm <- R.c
A1perm <- A1
Yperm <- Y

strata <- A1 == 1 # binary

s1 <- sum(strata)
s0 <- N - s1
for(i in 1:M){
  Rperm[strata] <- sample(R.c[strata])
  Rperm[!strata] <- sample(R.c[!strata])

  A1perm[strata] <- sample(A1[strata])
  A1perm[!strata] <- sample(A1[!strata])

  Yperm[strata] <- sample(Y[strata])
  Yperm[!strata] <- sample(Y[!strata])

  fit <- lm(Y ~ A1 * Rperm)
  coefs[i, ] <- coefficients(fit)

  mcov[, ,i] <- cov(data.frame(A1, Rperm, Y))
}
names(fit$coefficients)
colMeans(coefs)
apply(coefs, 2, sd)
apply(mcov, c(1,2), mean)

## Stratify on A1, permute Y -----------------------------------------------------
# cov(R*, Y) = 0; betaA1 correct
M <- 1e3
coefs <- matrix(0, M, 4)
mcov <- array(0, c(3,3,M))

Rperm <- R.c
A1perm <- A1
Yperm <- Y

strata <- A1 == 1 # binary

s1 <- sum(strata)
s0 <- N - s1
for(i in 1:M){
  Rperm[strata] <- sample(R.c[strata])
  Rperm[!strata] <- sample(R.c[!strata])

  A1perm[strata] <- sample(A1[strata])
  A1perm[!strata] <- sample(A1[!strata])

  Yperm[strata] <- sample(Y[strata])
  Yperm[!strata] <- sample(Y[!strata])

  fit <- lm(Yperm ~ A1*R.c)
  coefs[i, ] <- coefficients(fit)

  mcov[, ,i] <- cov(data.frame(A1, R, Yperm))
}
names(fit$coefficients)
colMeans(coefs)
apply(coefs, 2, sd)
apply(mcov, c(1,2), mean)


## Stratify on R, permute A1 -----------------------------------------------------
# cov(A1*, Y) != 0 ??, but coeff = 0
M <- 1e3
coefs <- matrix(0, M, 4)
mcov <- array(0, c(3,3,M))

Rperm <- R.c
A1perm <- A1
Yperm <- Y

strata <- R == 1 # binary

s1 <- sum(strata)
s0 <- N - s1
for(i in 1:M){
  Rperm[strata] <- sample(R.c[strata])
  Rperm[!strata] <- sample(R.c[!strata])

  A1perm[strata] <- sample(A1[strata])
  A1perm[!strata] <- sample(A1[!strata])

  Yperm[strata] <- sample(Y[strata])
  Yperm[!strata] <- sample(Y[!strata])

  fit <- lm(Y ~ A1perm*R.c)
  coefs[i, ] <- coefficients(fit)

  mcov[, ,i] <- cov(data.frame(A1perm, R, Y))
}
names(fit$coefficients)
colMeans(coefs)
apply(coefs, 2, sd)
apply(mcov, c(1,2), mean)


## Stratify on R, permute Y -----------------------------------------------------
# cov(A1*, Y) != 0 ??, but coeff = 0
M <- 1e3
coefs <- matrix(0, M, 4)
mcov <- array(0, c(3,3,M))

Rperm <- R.c
A1perm <- A1
Yperm <- Y

strata <- R == 1 # binary

s1 <- sum(strata)
s0 <- N - s1
for(i in 1:M){
  Rperm[strata] <- sample(R.c[strata])
  Rperm[!strata] <- sample(R.c[!strata])

  A1perm[strata] <- sample(A1[strata])
  A1perm[!strata] <- sample(A1[!strata])

  Yperm[strata] <- sample(Y[strata])
  Yperm[!strata] <- sample(Y[!strata])

  fit <- lm(Yperm ~ A1*R.c)
  coefs[i, ] <- coefficients(fit)

  mcov[, ,i] <- cov(data.frame(A1, R, Yperm))
}
names(fit$coefficients)
colMeans(coefs)
apply(coefs, 2, sd)
apply(mcov, c(1,2), mean)

## Stratify on Y, permute A1 -----------------------------------------------------
# cov(A1*, Y) != 0 ??, but coeff = 0
M <- 1e3
coefs <- matrix(0, M, 4)
mcov <- array(0, c(3,3,M))

Rperm <- R.c
A1perm <- A1
Yperm <- Y

strata <- Y == 1 # binary

s1 <- sum(strata)
s0 <- N - s1
for(i in 1:M){
  Rperm[strata] <- sample(R.c[strata])
  Rperm[!strata] <- sample(R.c[!strata])

  A1perm[strata] <- sample(A1[strata])
  A1perm[!strata] <- sample(A1[!strata])

  Yperm[strata] <- sample(Y[strata])
  Yperm[!strata] <- sample(Y[!strata])

  fit <- lm(Y ~ A1perm*R.c)
  coefs[i, ] <- coefficients(fit)

  mcov[, ,i] <- cov(data.frame(A1perm, R, Y))
}
names(fit$coefficients)
colMeans(coefs)
apply(coefs, 2, sd)
apply(mcov, c(1,2), mean)

## Stratify on Y, permute R -----------------------------------------------------
# cov(A1*, Y) != 0 ??, but coeff = 0
M <- 1e3
coefs <- matrix(0, M, 4)
mcov <- array(0, c(3,3,M))

Rperm <- R.c
A1perm <- A1
Yperm <- Y

strata <- Y == 1 # binary

s1 <- sum(strata)
s0 <- N - s1
for(i in 1:M){
  Rperm[strata] <- sample(R.c[strata])
  Rperm[!strata] <- sample(R.c[!strata])

  A1perm[strata] <- sample(A1[strata])
  A1perm[!strata] <- sample(A1[!strata])

  Yperm[strata] <- sample(Y[strata])
  Yperm[!strata] <- sample(Y[!strata])

  fit <- lm(Y ~ A1*Rperm)
  coefs[i, ] <- coefficients(fit)

  mcov[, ,i] <- cov(data.frame(A1, Rperm, Y))
}
names(fit$coefficients)
colMeans(coefs)
apply(coefs, 2, sd)
apply(mcov, c(1,2), mean)

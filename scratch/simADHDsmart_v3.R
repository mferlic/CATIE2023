
#' simADHDsmart
#'
#' @description
#' This function simulates synthetic data to match the Pelham ADHD SMART study.
#' Author: Mason Ferlic, December 2022, University of Michigan, D3C
#'
#' @details
#' @section Operating Characteristics:
#' TO DO
#'
#' @section Baseline Covariates:
#' @section First-Stage:
#' @section Second-Stage:
#'
#'
#' @param N number of observations to generate
#'
#' @returns data.frame
#' @export
#'



N = 1000
baseline.params = list(p.odd = 0.4, m.severity = 0, p.priormed = 0.3, p.race = 0.8)
Y0.coef = c(2, odd = -0.2, severity = -0.3)
U.params = list(mu = 0, sd = 1) # do not change
R.coef = c(-0.3, A1 = -0.2, "priormed:A1" = -0.3, U = 0.5) # probit model
adherence.coef = c(0.2, "priormed:A1" = -0.2, U = 0.5) # probit model

event_time.coef = NULL
Y1.coef = c(2.5, A1 = -0.3, U = 0.4)
Y2.baseline = c(3, odd = -0.3, severity = -0.4, priormed = 0, race = 0.5)
Y2.tx1 = c(A1 = 0.3, "priormed:A1" = -0.6)
Y2.n1 = c(R = 0.4, "A1:R" = -1.2, adherence = 0.1, U = 0.4)
Y2.tx2 = c(A2 = -0.3, "A1:A2" = -0.2, "adherence:A2" = 1.2)

# Baseline Covariates -----------------------------------------------------
# generator functions; centered
odd_f <- function(N){rbinom(N, 1, baseline.params$p.odd) -  baseline.params$p.odd}
severity_f <- function(N){rnorm(N, baseline.params$m.severity, 1) - baseline.params$m.severity}
priormed_f <- function(N){rbinom(N, 1, baseline.params$p.priormed) - baseline.params$p.priormed}
race_f <- function(N){rbinom(N, 1, baseline.params$p.race) - baseline.params$p.race}

# Internals
linearMult <- function(named.coefs, H) {
  nm <- names(named.coefs)
  if (nm[1] == "") {
    Int <- TRUE
    nm[1] <- "1"
  } else { Int <- FALSE}
  ff <- reformulate(nm, intercept = Int) # generate formula; look in caller environment. exclude first name which is intercept
  #env = rlang::caller_env()
  mt <- terms.formula(ff, keep.order = TRUE) # model terms
  mf <- model.frame.default(mt, data = H) # model frame; where variables are created. looks first in data then environment
  mm <- model.matrix(mt, data = H) # model matrix
  xb <- mm %*% named.coefs
  xb %>% as.vector()
}

sampleLogisticMean <- function(x, named.coefs, N = 1e6) {
  A1 <- rep(x, N)
  odd <- odd_f(N)
  severity <- severity_f(N)
  priormed <- priormed_f(N)
  race <- race_f(N)
  U <- U_f(N)
  H <- data.frame(A1, odd, severity, priormed, race, U)
  lo <- linearMult(named.coefs, H)
  p <- plogis(lo) # link function
  mean(p)
}

sampleProbitMean <- function(by, named.coefs, N = 1e6) {
  A1 <- A1_f(N)
  odd <- odd_f(N)
  severity <- severity_f(N)
  priormed <- priormed_f(N)
  race <- race_f(N)
  U <- U_f(N)
  H <- data.frame(A1, odd, severity, priormed, race, U)
  lat <- linearMult(named.coefs, H) + rnorm(N)
  D <- lat > 0
  data.frame(H, D) %>% group_by_at(vars(!!by)) %>% summarise(mean = mean(D))
}

### Integrate logistic distribution ??
f_u <- Vectorize(function(U, priormed, A1, named.coefs) {
  H <- data.frame(A1, priormed, U)
  plogis(linearMult(named.coefs, H)) * dnorm(U)
}, vectorize.args = "U")

dIntegrate <- function(f, z, p, ...) {
  f(z[1], ...)*p + f(z[2], ...)*(1-p)
}

g_pm_A1 <- function(priormed, A1, named.coefs) {integrate(f_u, -Inf, Inf, priormed=priormed, A1=A1, named.coefs)$value}

f_A1 <- function(A1, named.coefs) {dIntegrate(g_pm_A1, z=c(1, 0)-baseline.params$p.priormed, p=baseline.params$p.priormed, A1, named.coefs)}
ER_A1 <- function(A1) {f_A1(A1, R.coef)}
ER_A1(1)

### Integrate Probit distribution ??
f_u <- Vectorize(function(U, priormed, A1, named.coefs) {
  H <- data.frame(A1, priormed, U)
  pnorm(linearMult(named.coefs, H)) * dnorm(U)
}, vectorize.args = "U")

dIntegrate <- function(f, z, p, ...) {
  f(z[1], ...)*p + f(z[2], ...)*(1-p)
}

g_pm_A1 <- function(priormed, A1, named.coefs) {integrate(f_u, -Inf, Inf, priormed=priormed, A1=A1, named.coefs)$value}

f_A1 <- function(A1, named.coefs) {dIntegrate(g_pm_A1, z=c(1, 0)-baseline.params$p.priormed, p=baseline.params$p.priormed, A1, named.coefs)}
ER_A1 <- function(A1) {f_A1(A1, R.coef)}
ER_A1(1)

# centered
odd <- odd_f(N)
odd.nc <- odd +  baseline.params$p.odd
severity <- severity_f(N)
severity.nc <- severity + baseline.params$m.severity
priormed <- priormed_f(N)
priormed.nc <- priormed + baseline.params$p.priormed
race <- race_f(N)
race.nc <- race + baseline.params$p.race

H0.c <- data.frame(odd, severity, priormed, race)

Y0 <- linearMult(Y0.coef, H0.c) + rnorm(N)

#Y02 <- cbind(1, odd, severity) %*% c(2, -0.2, -0.3)

# First-stage -------------------------------------------------------------
#Med: -1, BMOD: 1
A1_f <- function(N){2 * rbinom(N, 1, 0.5) - 1}
A1 <- A1_f(N)

## Unknown
U_f <- function(N){rnorm(N, U.params$mu, U.params$sd)} # common cause of R, Adh, Y1, Y2
U <- U_f(N)

H1.c <- mutate(H0.c, A1, U)

## Response OC: Med (-1) positively associated with Response, priorMed + association
lo.R <- linearMult(R.coef, H1.c)
lo.R2 <- cbind(1, A1, priormed*A1, U) %*% c(Int = -0.7, A1 = -0.2, "pm:A1" = -0.4, U = 0.5)  # log-odds scale
p.R <- plogis(lo.R)
R <- rbinom(N, 1, p.R)
R.resid <- R - p.R

# E[R | A1]
a1 <- c(-1, 1)
ER.A1 <- sapply(a1, sampleLogisticMean, named.coefs = R.coef)

## R Probit model
R.lat <- linearMult(R.coef, H1.c) + rnorm(N)
R <- as.numeric(R.lat > 0)
p.R <- pnorm(R.lat)
R.resid <- R - p.R

# Integrate out U
sd <- sqrt(1 + R.coef['U']^2 * U.params$sd^2)
EU <- U.params$mu
Z <- linearMult(R.coef, H1.c %>% mutate(U = EU)) / sd
ER.A1.X <- pnorm(Z)

data.frame(ER.A1.X, A1) %>% group_by(A1) %>% summarise(mean(ER.A1.X)) # yay agrees!

ER.A1 <- sampleProbitMean("A1", R.coef)
ER.A1

## Adherence OC: For those on PriorMed, Med(-1) positively associated with Adherence
lo.adherence <- linearMult(adherence.coef, H1.c)
#lo.adherence2 <- cbind(1, priormed * A1, U) %*% c(Int = -0.2, "pm:A1" = -0.3, U = 0.2)  # log-odds scale
p.adherence <- plogis(lo.adherence)
adherence <- rbinom(N, 1, p.adherence)
adherence.resid <- adherence - p.adherence

# E[Adh | A1]
EAdh.A1 <- sapply(a1, sampleLogisticMean, named.coefs = adherence.coef)

## Adherence Probit Model
adherence.lat <- linearMult(adherence.coef, H1.c) + rnorm(N)
adherence <- as.numeric(adherence.lat > 0)
p.adherence <- pnorm(adherence.lat)
adherence.resid <- adherence - p.adherence
mean(adherence)

# Integrate out U
sd <- sqrt(1 + adherence.coef['U']^2 * U.params$sd^2)
Z <- linearMult(adherence.coef, H1.c %>% mutate(U = EU)) / sd
ES.A1.X <- pnorm(Z)
data.frame(ES.A1.X, A1) %>% group_by(A1) %>% summarise(mean(ES.A1.X)) # yay agrees!
ES.A1 <- sampleProbitMean("A1", adherence.coef)
ES.A1

## Time to Non-response (Event-time) OC: No association for now
event_time <- sample.int(8, N, replace = TRUE)  # uniform; for now
event_time[R == 1] <- NA

H1.c <- mutate(H1.c, R = R.resid, adherence = adherence.resid, event_time) # use centered R and Adherence!

## First-stage Outcome OC: Small change from baseline, Med (-1) is better initially
Y1 <- linearMult(Y1.coef, H1.c) + rnorm(N)

H1.c <- mutate(H1.c, Y1)


# Second-stage ------------------------------------------------------------

#AUG: -1, INT: 1

A2_f <- function(N){2 * rbinom(N, 1, 0.5) - 1}
A2 <- A2_f(N)
A2[R == 1] <- 0  # Fill NA at end

H2.c <- mutate(H1.c, A2)

# EOS Outcome -----------------------------------------------------------------

## Baseline
Y2 <- linearMult(Y2.baseline, H2.c)
#Y22 <- cbind(1, odd, severity, priormed, race) %*% c(3, -0.3, -0.4, 0, 0.5)

## First-stage OC: BMOD is marginally better in the long run. For those on PriorMed, Med(-1) is better. For those on PriorMed and
# Responded, Med(-1) is MUCH better.
tx1 <- linearMult(Y2.tx1, H2.c)  # small main effect!
Y2 <- Y2 + tx1

## Intermediate Associations
n1 <- linearMult(Y2.n1, H2.c)
Y2 <- Y2 + n1

## Second-stage: For non-responders OC: AUG(-1) is marginally better. The effect of ADD(-1) is neg for those that
## started with Med(-1). For those Adherent, INT is MUCH better.
#tx2 <- (1 - R) * cbind(A2, A1 * A2, adherence * A2) %*% c(-0.3, -0.2, 1.2)
tx2 <- (1 - R) * linearMult(Y2.tx2, H2.c)
Y2 <- Y2 + tx2

## Error
Y2 <- Y2 + rnorm(N, 0, 0.5) # determine what error is needed for sig effects

H3.c <- mutate(H2.c, Y2)

A2[R == 1] <- NA
event_time[R == 1] <- NA

## Marginal Model
marginalMeanModel <- function(a1, a2) {
  A1 <-  a1
  A2 <-  a2
  odd <-  0
  severity <-  0
  priormed <- 0
  race <- 0
  adherence <-  0
  R <- 0
  U <- 0
  EH2.c <- data.frame(odd, severity, priormed, race, A1, U, R, adherence, A2)
  linearMult(Y2.baseline, EH2.c) + linearMult(Y2.tx1, EH2.c) + (1 - ER.A1[-a1]) * linearMult(Y2.tx2, EH2.c)
  Y2.baseline[1] + Y2.tx1['A1']*a1 + (1 - ER.A1[-a1])*(Y2.tx2['A2'] + Y2.tx2['A1:A2']*a1)*a2
}

DTRmeans <- expand.grid(a1 = c(-1,1), a2 = c(-1,1))
DTRmeans %>% rowwise() %>% mutate(means = marginalMeanModel(a1, a2)[1])

data <- data.frame(ID = 1:N, odd = odd.nc, severity = severity.nc,
                   priormed = priormed.nc, race = race.nc, Y0,
                   A1, R, event_time, adherence, Y1, A2, Y2)
str(data)





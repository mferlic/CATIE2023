
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


simADHDsmart <- function(N = 150,
                         baseline.params = list(p.odd = 0.4, m.severity = 0, p.priormed = 0.3, p.race = 0.8),
                         Y0.coef = c(2, odd = -0.2, severity = -0.3),
                         U.params = list(mu = 0, sd = 1), # do not change
                         R.coef = c(-0.3, A1 = -0.1, "priormed:A1" = -0.2, U = 0.5), # probit model
                         adherence.coef = c(-0.1, "priormed:A1" = -0.2, U = 0.4), # probit model
                         event_time.coef = NULL,
                         Y1.coef = c(2.5, A1 = -0.3, U = 0.4),
                         Y2.baseline = c(3, odd = -0.3, severity = -0.4, priormed = 0, race = 0.5),
                         Y2.tx1 = c(A1 = 0.2, "priormed:A1" = -0.6),
                         Y2.n1 = c(R.resid = 0.4, adherence.resid = 0.2, U = 0.3),
                         Y2.tx2 = c(A2 = -0.3, "A1:A2" = -0.1, "adherence:A2" = 1.2)) {


# Generator functions -----------------------------------------------------
odd_f <- function(N){rbinom(N, 1, baseline.params$p.odd) -  baseline.params$p.odd}
severity_f <- function(N){rnorm(N, baseline.params$m.severity, 1) - baseline.params$m.severity}
priormed_f <- function(N){rbinom(N, 1, baseline.params$p.priormed) - baseline.params$p.priormed}
race_f <- function(N){rbinom(N, 1, baseline.params$p.race) - baseline.params$p.race}

U_f <- function(N){rnorm(N, 0, U.params$sd)} # common cause of R, Adh, Y1, Y2
A1_f <- function(N){2 * rbinom(N, 1, 0.5) - 1}
A2_f <- function(N){2 * rbinom(N, 1, 0.5) - 1}

## Internals
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

sampleProbitMean <- function(named.coefs, N = 1e6) {
  A1 <- A1_f(N)
  odd <- odd_f(N)
  severity <- severity_f(N)
  priormed <- priormed_f(N)
  race <- race_f(N)
  U <- U_f(N)
  H <- data.frame(A1, odd, severity, priormed, race, U)
  lat <- linearMult(named.coefs, H) + rnorm(N)
  D <- as.numeric(lat > 0)
  data.frame(H, D)
}

# Baseline Covariates -----------------------------------------------------
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

# First-stage -------------------------------------------------------------
#Med: -1, BMOD: 1
A1 <- A1_f(N)

## Unknown
U <- U_f(N)

H1.c <- mutate(H0.c, A1, U)

## Response OC: Med (-1) positively associated with Response, priorMed + association
## R Probit model
R.lat <- linearMult(R.coef, H1.c) + rnorm(N)
R <- as.numeric(R.lat > 0)
p.R <- pnorm(R.lat)
R.resid <- R - p.R

# Integrate out U
sd <- sqrt(1 + R.coef['U']^2 * U.params$sd^2)
EU <- 0
Z <- linearMult(R.coef, H1.c %>% mutate(U = EU)) / sd
ER.A1.X <- pnorm(Z) # analytic expression

ER.H <- sampleProbitMean(R.coef) # sample from the posterior, tibble
ER.A1 <- ER.H %>% group_by(A1) %>% summarize(mean = mean(D))
ER_a1 <- Vectorize(function(a1) {
  ER.A1$mean[ER.A1$A1 == a1]
}) # function if we want to center by first stage a1
ER <- mean(ER.H$D)
R.c = R - ER # grand mean centered

## Adherence (S) OC: For those on PriorMed, Med(-1) positively associated with Adherence
## Adherence Probit Model
adherence.lat <- linearMult(adherence.coef, H1.c) + rnorm(N)
adherence <- as.numeric(adherence.lat > 0)
p.adherence <- pnorm(adherence.lat)
adherence.resid <- adherence - p.adherence

# Integrate out U
sd <- sqrt(1 + adherence.coef['U']^2 * U.params$sd^2)
Z <- linearMult(adherence.coef, H1.c %>% mutate(U = EU)) / sd
ES.A1.X <- pnorm(Z)

ES.H <- sampleProbitMean(adherence.coef) # sample from the posterior, tibble
ES.A1 <- ES.H %>% group_by(A1) %>% summarize(mean = mean(D))
ES_a1 <- Vectorize(function(a1) {
  ES.A1$mean[ES.A1$A1 == a1]
}) # function if we want to center by first stage a1
ES <- mean(ES.H$D)
adherence.c = adherence - ES # grand mean centered

## Time to Non-response (Event-time) OC: No association for now
event_time <- sample.int(8, N, replace = TRUE)  # uniform; for now
event_time[R == 1] <- NA

H1.c <- mutate(H1.c, R, R.resid, adherence = adherence, adherence.resid, event_time) # use grand mean centered Adherence!

## First-stage Outcome OC: Small change from baseline, Med (-1) is better initially
Y1 <- linearMult(Y1.coef, H1.c) + rnorm(N)

H1.c <- mutate(H1.c, Y1)

# Second-stage ------------------------------------------------------------

#AUG: -1, INT: 1
A2 <- A2_f(N)
A2[R == 1] <- 0  # Fill NA at end

H2.c <- mutate(H1.c, A2)

# EOS Outcome -----------------------------------------------------------------

## Baseline
Y2 <- linearMult(Y2.baseline, H2.c)

## First-stage OC: BMOD is marginally better in the long run. For those on PriorMed, Med(-1) is better. For those on PriorMed and
# Responded, Med(-1) is MUCH better.
tx1 <- linearMult(Y2.tx1, H2.c)  # small main effect!
Y2 <- Y2 + tx1

## Intermediate Associations
n1 <- linearMult(Y2.n1, H2.c)
Y2 <- Y2 + n1

## Second-stage: For non-responders OC: AUG(-1) is marginally better. The effect of AUG(-1) is neg for those that
## started with Med(-1). For those Adherent, INT is MUCH better.
tx2 <- (1 - R) * linearMult(Y2.tx2, H2.c)
Y2 <- Y2 + tx2

## Error
Y2 <- Y2 + rnorm(N, 0, 0.5) # determine what error is needed for sig effects

H3.c <- mutate(H2.c, Y2)

A2[R == 1] <- NA
event_time[R == 1] <- NA


# Compute Marginal Means --------------------------------------------------
getMarginalMeans <- function(a1, a2) {
  A1 <-  a1
  A2 <-  a2
  odd <-  0
  severity <-  0
  priormed <- 0
  race <- 0
  adherence.resid <- 0 # centered residuals
  R.resid <- 0 # centered residuals
  R <- 0
  adherence <- ES_a1(a1) # difference in adherence between levels of a1 since we are grand mean centering all moderators
  event_time <- 0
  U <- 0

  EH2.c <- data.frame(odd, severity, priormed, race, A1, U, R, R.resid, adherence, adherence.resid, A2)
  DTR <- linearMult(Y2.baseline, EH2.c) + linearMult(Y2.tx1, EH2.c) + (1 - ER_a1(a1)) * linearMult(Y2.tx2, EH2.c)
  DTR
}

DTRmeans <- expand.grid(a1 = c(-1,1), a2 = c(-1,1)) %>%
  rowwise() %>%
  mutate(mean = getMarginalMeans(a1, a2))

data <- data.frame(ID = 1:N, odd = odd.nc, severity = severity.nc,
                   priormed = priormed.nc, race = race.nc, Y0,
                   A1, R, event_time, adherence, Y1, A2, Y2)
return(list(data = data, DTRmeans = DTRmeans))

}




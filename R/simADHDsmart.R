#' simADHDsmart
#'
#' @author Mason Ferlic, December 2022, University of Michigan, D3C
#'
#' @description
#' This function simulates synthetic data to mirror desired characteristics of the Pelham ADHD SMART study.
#' Created specifically for CATIE 2023.
#'
#'
#' @param N number of observations to generate
#' @param baseline.params A list specifying the probability `odd`, mean `severity`, prob. `priormed`, prob. `race`
#' @param Y0.coef A named vector specifying the linear model coefficients for baseline school performance. Can be a function of any past variables. See *Named Vectors*
#' @param U.params A vector specifying the normal distribution parameters `mu` and `sd`
#' @param R.coef A named vector specifying the linear **probit model** coefficients for the probability of being a responder to first-stage treatment. Can be a function of any past variables.
#' @param adherence.coef A named vector specifying the linear **probit model** coefficients for the probability of being adherent to first-stage treatment. Can be a function of any past variables.
#' @param NRtime.coef not implemented
#' @param Y1.coef A named vector of linear model coefficients specifying the first-stage treatment causal effect on \eqn{Y_1} school performance
#' @param Y2.baseline A named vector of linear model coefficients specifying the baseline covariate associations on end-of-study \eqn{Y_2} school performance
#' @param Y2.tx1 A named vector of linear model coefficients specifying the first-stage treatment causal effect on end-of-study \eqn{Y_2} school performance. Can be a function of any baseline moderators.
#' @param Y2.n1 A named vector of linear model coefficients specifying the nuisance associations on end-of-study \eqn{Y_2} school performance. Can be a function of any prior moderators. NOTE: must specify R.resid and adherence.resid which are the direct associations in the SNMM.
#' @param Y2.tx2 A named vector of linear model coefficients specifying the second-stage treatment causal effect, among non-responders, on end-of-study \eqn{Y_2} school performance. Can be a function of any baseline moderators. NOTE: all moderators are grand mean centered
#' @param sigma gaussian noise added to \eqn{Y_{0,1,2}}
#'
#' @return A list with components
#' \describe{
#'  \item{data}{data.frame of observed variables}
#'  \item{DTRmean}{marginal embedded DTR means}
#'  }
#'
#' @details # Named Vectors
#' The named vectors of coefficients use formula notation to specify interactions i.e. "`priormed:A1`", which must be quoted. The order does not matter but spelling does. Note: if the first element is unamed it is treated as the intercept term, else if all elements are named the intercept is ommited.
#'
#' # Operating Characteristics
#' All baseline and moderator variable are grand mean centered. The nuisance terms `R.resid` and `adherence.resid` are residualized using the true probabilities.
#'
#' ## Default Structural Nested Mean Model:
#' ### School performance after first-stage
#' \deqn{Y_{1}(a_1, a_2) = 2.5 - 0.3a_1 + 0.4U + \epsilon}
#' ### School performance after second-stage
#' \deqn{Y_{2}(a_1, a_2) = 3 - 0.2odd - 0.3severity + 0*priormed + .5race + \\ (0.3 - 1.4priormed)a_1 + \\ 0.8R.resid(a_1, u) + 1 adherence.resid(a_1, u) \\ +(1-R)(-0.3 + 0.1a_1 + 1.2adherence)a_2 + 0.4U + \epsilon}
#'
#' ## Baseline Covariates:
#' -  `odd`: binary, centered
#' -  `severity`: standard normal
#' -  `priormed`: binary, centered
#' -  `race`: binary, centered
#' -  `U`: an unknown, common cause of `R`, `adherence`, `Y1`, and `Y2`. Induces collider bias if naively conditioning on the time-varying covariates.
#'
#' ## First-Stage:
#' -  Marginally, `BMOD(1)` is initially worse compared to `MED(-1)` after the first-stage \eqn{Y_1}, but better in the long run \eqn{Y_2}
#' -  IF on `priormed`, starting with `MED` is better
#'
#' ## Nuisance Associations
#' -  `R`: positively associated with outcome
#' -  `adherence`: positively associated with outcome
#' -  `NRtime`: or time to non-response has no effect
#'
#' ## Second-Stage (among non-responders only):
#' -  Marginally, `AUG(-1)` is better compared to `INT(1)`
#' -  Positive interaction for `AUG` if given `MED`, and `INT` if given `BMOD`
#' -  IF adherent to first-stage, much better to `INT`; if non-adherent to first-stage, much better to `AUG`
#'
#' @details # Internal functions
#' -  linearMult: takes a named vector of coefficients and a data.frame and creates the necessary design matrix to multiply and return the function values.
#' - sampleProbitMean: generates observations from the posterior in order to empirically evaluate probability the integral.
#' -  getMarginalMeans: takes the specified coefficients and returns the marginal means of the four embedded adaptive interventions, averaging over response and adherence.
#'
#' @import dplyr
#'
#' @export
#'


simADHDsmart <- function(N = 150,
                         baseline.params = list(p.odd = 0.4, m.severity = 0, p.priormed = 0.3, p.race = 0.8),
                         Y0.coef = c(2, odd = -0.2, severity = -0.3),
                         U.params = c(mu = 0, sd = 1), # do not change
                         R.coef = c(-0.4, A1 = -0.1, "priormed:A1" = -0.2, U = 0.2), # probit model
                         adherence.coef = c(-0.1, "priormed:A1" = -0.2, U = 0.1), # probit model
                         NRtime.coef = NULL,
                         Y1.coef = c(2.5, A1 = -0.3, U = 0.1),
                         Y2.baseline = c(3, odd = -0.3, severity = -0.4, priormed = 0, race = 0.5),
                         Y2.tx1 = c(A1 = 0.3, "priormed:A1" = -1.4),
                         Y2.n1 = c(R.resid = 0.8, adherence.resid = 1 , U = 0.4),
                         Y2.tx2 = c(A2 = -0.3, "A1:A2" = 0.1, "adherence:A2" = 1.2),
                         sigma = 1) {


# Generator functions -----------------------------------------------------
odd_f <- function(N){rbinom(N, 1, baseline.params$p.odd) -  baseline.params$p.odd}
severity_f <- function(N){rnorm(N, baseline.params$m.severity, 1) - baseline.params$m.severity}
priormed_f <- function(N){rbinom(N, 1, baseline.params$p.priormed) - baseline.params$p.priormed}
race_f <- function(N){rbinom(N, 1, baseline.params$p.race) - baseline.params$p.race}

U_f <- function(N){rnorm(N, U.params[1], U.params[2])} # common cause of R, Adh, Y1, Y2
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

Y0 <- linearMult(Y0.coef, H0.c) + rnorm(N, 0, sigma)

# First-stage -------------------------------------------------------------
#Med: -1, BMOD: 1
A1 <- A1_f(N)

## Unknown
U <- U_f(N)

H1.c <- dplyr::mutate(H0.c, A1, U)

## Response OC: Med (-1) positively associated with Response, priorMed + association
## R Probit model
R.lat <- linearMult(R.coef, H1.c) + rnorm(N)
R <- as.numeric(R.lat > 0)
p.R <- pnorm(R.lat)
R.resid <- R - p.R

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

ES.H <- sampleProbitMean(adherence.coef) # sample from the posterior, tibble
ES.A1 <- ES.H %>% group_by(A1) %>% summarize(mean = mean(D))
ES_a1 <- Vectorize(function(a1) {
  ES.A1$mean[ES.A1$A1 == a1]
}) # function if we want to center by first stage a1
ES <- mean(ES.H$D)
adherence.c = adherence - ES # grand mean centered
adherence.a1 = adherence - ES_a1(A1)

## Time to Non-response (Event-time) OC: No association for now
NRtime <- sample.int(8, N, replace = TRUE)  # uniform; for now
NRtime[R == 1] <- NA

H1.c <- dplyr::mutate(H1.c, R, R.resid, adherence = adherence.c, adherence.resid, NRtime) # use grand mean centered Adherence!

## First-stage Outcome OC: Small change from baseline, Med (-1) is better initially
Y1 <- linearMult(Y1.coef, H1.c) + rnorm(N, 0, sigma)

H1.c <- dplyr::mutate(H1.c, Y1)

# Second-stage ------------------------------------------------------------

#AUG: -1, INT: 1
A2 <- A2_f(N)
A2[R == 1] <- 0  # Fill NA at end

H2.c <- dplyr::mutate(H1.c, A2)

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
Y2 <- Y2 + rnorm(N, 0, sigma) # determine what error is needed for sig effects

H3.c <- dplyr::mutate(H2.c, Y2)

A2[R == 1] <- NA
NRtime[R == 1] <- NA


# Compute Marginal Means --------------------------------------------------
getMarginalMeans <- function(a1, a2) {
  A1 <-  a1
  A2 <-  a2

  # Mean values
  odd <-  0
  severity <-  0
  priormed <- 0
  race <- 0
  adherence.resid <- 0 # centered residuals
  R.resid <- 0 # centered residuals
  R <- 0
  adherence <- ES_a1(a1) - ES # since we are grand mean centering all moderators, we need the expectation conditioned on A1... basically 0
  NRtime <- 0
  U <- 0

  EH2.c <- data.frame(odd, severity, priormed, race, A1, U, R, R.resid, adherence, adherence.resid, A2)
  DTR <- linearMult(Y2.baseline, EH2.c) + linearMult(Y2.tx1, EH2.c) + (1 - ER_a1(a1)) * linearMult(Y2.tx2, EH2.c)
  DTR
}

DTRmeans <- expand.grid(a1 = c(-1,1), a2 = c(-1,1)) %>%
  rowwise() %>%
  dplyr::mutate(mean = getMarginalMeans(a1, a2))

data <- data.frame(ID = 1:N, odd = odd.nc, severity = severity.nc,
                   priormed = priormed.nc, race = race.nc, Y0,
                   A1, R, NRtime, adherence, Y1, A2, Y2)
return(list(data = data, DTRmeans = DTRmeans))

}




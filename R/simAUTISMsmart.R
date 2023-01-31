#' simAUTISMsmart
#'
#' @author Mason Ferlic, December 2022, University of Michigan, D3C
#'
#' @description
#' This function simulates synthetic data to roughly mirror desired characteristics of the AUTISM SMART study by Kasari.
#' Created specifically for CATIE 2023.
#'
#'
#' @param N number of observations to generate.
#' @param baseline.params A named list specifying the `unif` parameters for O11 and `beta` parameters for O12.
#' @param Y0.coef A named vector specifying the linear model coefficients for score. Can be a function of any past variables. See *Named Vectors*.
#' @param U.params A vector specifying the normal distribution parameters `mu` and `sd`.
#' @param R.coef A named vector specifying the **linear probit model** coefficients for the probability of being a responder to first-stage treatment. Can be a function of any past variables.
#' @param O21.coef A named vector specifying the **linear model** coefficients for O21. Can be a function of any past variables.
#' @param O22.coef A named vector specifying the **linear model** coefficients for O22. Can be a function of any past variables.
#' @param Y1.coef A named vector of **linear model** coefficients specifying the first-stage treatment causal effect on \eqn{Y_1}.
#' @param Y2.baseline A named vector of **linear model** coefficients specifying the baseline covariate associations on end-of-study \eqn{Y_2}.
#' @param Y2.tx1 A named vector of **linear model** coefficients specifying the first-stage treatment causal effect on end-of-study \eqn{Y_2} school performance. Can be a function of any baseline moderators.
#' @param Y2.n1 A named vector of **linear model** coefficients specifying the nuisance associations on end-of-study \eqn{Y_2}. Can be a function of any prior moderators. NOTE: must specify O21.resid and O22.resid which are the residualized direct associations in the SNMM.
#' @param Y2.tx2 A named vector of **linear model** coefficients specifying the second-stage treatment causal effect, among non-responders to a1 == 1, on end-of-study \eqn{Y_2}. Can be a function of any baseline moderators. NOTE: all moderators are residualized centered
#' @param sigma Gaussian noise added to \eqn{Y_{0,1,2}}
#'
#'
#' @return A list with components
#' \describe{
#'  \item{data}{data.frame of observed variables}
#'  \item{DTRmean}{marginal embedded DTR means}
#'  }
#'
#' @details # Named Vectors
#' The named vectors of coefficients use formula notation to specify interactions i.e. "`A1:O11`", which must be quoted. The order does not matter but spelling does. Note: if the first element is unamed it is treated as the intercept term, else if all elements are named the intercept is omitted.
#'
#' # Operating Characteristics
#' All baseline variables are grand mean centered. All moderators are centered conditional on the expectation given the past (residuals).
#'
#' ## Default Structural Nested Mean Model:
#' ### Y1 after first-stage
#' \deqn{Y_{1,i}(a_1, a_2) = 50 + 1a_1 + 4U_i + \epsilon_i}
#' ### Y2 after second-stage
#' \deqn{Y_{2,i}(a_1, a_2) = 60 + 0.3\bar{O}_{11} + (1 + 0.3\bar{O}_{11})a_1 + \\ 18\tilde{R} + 5\tilde{O}_{21} \\ +I(a_1 =1)(1-R)(-4 + 6\tilde{O}_{21})a_2 + 6U_i + \epsilon_i}
#'
#'
#' # Internal functions
#' -  linearMult: takes a named vector of coefficients and a data.frame and creates the necessary design matrix to multiply and return the function values.
#' - sampleProbitMean: generates observations from the posterior in order to empirically evaluate probability the integral.
#' -  getMarginalMeans: takes the specified coefficients and returns the marginal means of the four embedded adaptive interventions, averaging over response and adherence.
#'
#' @import dplyr
#'
#' @export
#'


simAUTISMsmart <- function(N = 200,
                         baseline.params = list(O11 = c(lo = 0, up = 80), O12 = c(2,4)), #unif, beta
                         Y0.coef = c(40, O11 = 0.4, O12 = 0), # linear
                         U.params = c(mu = 0, sd = 1),
                         R.coef = c(0.25, A1 = -0.3, U = 0.1), # probit model
                         O21.coef = c(5.6, A1 = -0.4, R = 2.5, U = 5), # linear model
                         O22.coef = c(48, A1 = 3, R = 8, "A1:R" = -4, "O21"= -2, U = 4), # linear model f(R, O21, O22)
                         Y1.coef = c(50, A1 = 1, U = 4),
                         Y2.baseline = c(60, O11 = 0.3),
                         Y2.tx1 = c(A1 = 1, "A1:O11" = 0.3), # grand mean centered covariates
                         Y2.n1 = c(R.resid = 18, O21.resid = 5, U = 6),
                         Y2.tx2 = c(A2 = -4, "A2:O21" = 6), # residual centered covariates
                         sigma = 8) {


# Generator functions -----------------------------------------------------
O11_f <- function(N) {runif(N, baseline.params$O11[1], baseline.params$O11[2])}
O12_f <- function(N) {50 * rbeta(N, baseline.params$O12[1], baseline.params$O12[2])}

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
  O11 <- O11_f(N)
  O11.c <- O11 - mean(O11)
  O12 <- O12_f(N)
  O12.c <- O12 - mean(O12)
  U <- U_f(N)
  H <- data.frame(A1, O11 = O11.c, O12 = O12.c, U)
  lat <- linearMult(named.coefs, H) + rnorm(N)
  D <- as.numeric(lat > 0)
  data.frame(H, D)
}

# Baseline Covariates -----------------------------------------------------
# centered
O11 <- O11_f(N)
EO11 <- (baseline.params$O11[1] + baseline.params$O11[2])/2 # mean unif
O11.c <- O11 - EO11
O12 <- O12_f(N)
EO12 <- 50 * baseline.params$O12[1] / (baseline.params$O12[1] + baseline.params$O12[2]) # mean beta
O12.c <- O12 - EO12

H0.c <- data.frame(O11 = O11.c, O12 = O12.c) # use centered baseline

Y0 <- linearMult(Y0.coef, H0.c) + rnorm(N, 0, sigma)

# First-stage -------------------------------------------------------------
A1 <- A1_f(N)

## Unknown
U <- U_f(N)

H1.c <- dplyr::mutate(H0.c, A1, U)

## Response OC
## R Probit model: using centered covariates
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

H1.c <- dplyr::mutate(H0.c, A1, U, R = R.c, R.resid) # use centered R

## O21
O21e <- rnorm(N, 0, 2)
O21 <- linearMult(O21.coef, H1.c) + O21e
O21.resid <- O21e
EO21_a1 <- function(a1) {linearMult(O21.coef, data.frame(A1 = a1, R = ER_a1(a1) - ER, U = 0))}
EO21_a1r <- function(a1, r) {linearMult(O21.coef, data.frame(A1 = a1, R = r, U = 0))}
EO21 <- O21.coef[1] # grand mean (if using centered coeffecients!)
O21.c <- O21 - EO21

H1.c <- dplyr::mutate(H0.c, A1, U, R = R.c, R.resid, O21 = O21.c, O21.resid) # use centered O21.c for O22

## O22
O22e <- rnorm(N, 0, 17)
O22 <- linearMult(O22.coef, H1.c) + O22e
O22.resid <- O22e
EO22_a1 <- function(a1) {linearMult(O22.coef, data.frame(A1 = a1, R = ER_a1(a1) - ER, O21 = EO21_a1(a1) - EO21, U = 0))}
EO22 <- O22.coef[1]
O22.c <- O22 - EO22

H1.c <- dplyr::mutate(H0.c, A1, U, R = R.c, R.resid, O21 = O21.c, O21.resid, O22 = O22.c, O22.resid) # use grand mean centered covariates for Y

## First-stage Outcome OC:
Y1 <- linearMult(Y1.coef, H1.c) + rnorm(N, 0, sigma)

H1.c <- dplyr::mutate(H1.c, Y1, O21 = O21.resid, O22 = O22.resid) # Use residuals for 2nd stage moderators

# Second-stage ------------------------------------------------------------

#AUG: -1, INT: 1
A2 <- A2_f(N)
Rerand <- (1 - R) * (A1 == 1) # Only A1==1, non-responders Rerand to second stage option
A2[!Rerand] <- 0  # Fill NA at end

H2.c <- dplyr::mutate(H1.c, A2, Rerand)

# EOS Outcome -----------------------------------------------------------------

## Baseline
Y2 <- linearMult(Y2.baseline, H2.c)

## First-stage OC:
tx1 <- linearMult(Y2.tx1, H2.c)
Y2 <- Y2 + tx1

## Intermediate Associations
n1 <- linearMult(Y2.n1, H2.c)
Y2 <- Y2 + n1

## Second-stage:
tx2 <- Rerand * linearMult(Y2.tx2, H2.c)
Y2 <- Y2 + tx2

## Noise
Y2 <- Y2 + rnorm(N, 0, sigma) # determine what noise is needed for sig effects

H3.c <- dplyr::mutate(H2.c, Y2)

A2[!Rerand] <- NA

# Compute Marginal Means --------------------------------------------------
getMarginalMeans <- function(a1, a2) {
  EH2.c <- data.frame(O11=0, O12=0, A1=a1, U=0, O21=0, O22=0, A2=a2)
  DTR <- linearMult(Y2.baseline, EH2.c) + linearMult(Y2.tx1, EH2.c) + (a1 == 1) * ER_a1(a1) * linearMult(Y2.tx2, EH2.c)
  DTR
}

DTRmeans <- data.frame(a1 = c(1,1,-1), a2 = c(1,-1,0)) %>%
  rowwise() %>%
  dplyr::mutate(mean = getMarginalMeans(a1, a2))

data <- data.frame(ID = 1:N, O11, O12, Y0, A1, R, O21, O22, Y1, A2, Y2)
return(list(data = data, DTRmeans = DTRmeans))

}




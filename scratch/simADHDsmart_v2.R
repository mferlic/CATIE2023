
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



simADHDsmart_v2 <- function(N = 100,
                         baseline.params = list(p.odd = 0.4, m.severity = 0, p.priormed = 0.3, p.race = 0.8),
                         Y0.coef = c(2, odd = -0.2, severity = -0.3),
                         U.params = list(mu = 0, sd = 1),
                         R.coef = c(-0.7, A1 = -0.2, "priormed:A1" = -0.4, U = 0.5),
                         adherence.coef = c(-0.2, "priormed:A1" = -0.3, U = 0.2),
                         ET.coef = NULL,
                         Y1.coef = c(2.5, A1 = -0.3, U = 0.4),
                         Y2.baseline = c(3, odd = -0.3, severity = -0.4, priormed = 0, race = 0.5),
                         Y2.tx1 = c(A1 = 0.2, "priormed:A1" = -0.6, "priormed:R:A1" = -1.2),
                         Y2.n1 = c(R = 0.4, U = 0.4),
                         Y2.tx2 = c(A2 = -0.3, "A1:A2" = -0.2, "adherence:A2" = 1.2)) {

  # Baseline Covariates -----------------------------------------------------
  # generator functions; centered
  odd_f <- function(N){rbinom(N, 1, baseline.params$p.odd) -  baseline.params$p.odd}
  severity_f <- function(N){rnorm(N, baseline.params$m.severity, 1) - baseline.params$m.severity}
  priormed_f <- function(N){rbinom(N, 1, baseline.params$p.priormed) - baseline.params$p.priormed}
  race_f <- function(N){rbinom(N, 1, baseline.params$p.race) - baseline.params$p.race}

  # Internals
  linMult <- function(named.coefs) {
    ff <- reformulate(names(named.coefs)[-1], env = rlang::caller_env()) # generate formula; look in caller environment. exclude first name which is intercept
    mt <- terms.formula(ff, keep.order = TRUE) # model terms
    mf <- model.frame.default(mt) # model frame; where variables are created
    mm <- model.matrix(mt, data = mf) # model matrix
    xb <- mm %*% named.coefs
    xb
  }

  sampleLogisticMean <- function(x, named.coefs, N = 1e6) {
    A1 <- rep(x, N)
    odd <- odd_f(N)
    severity <- severity_f(N)
    priormed <- priormed_f(N)
    race <- race_f(N)
    U <- U_f(N)

    lo <- linMult(named.coefs)
    p <- plogis(lo) # link function
    mean(p)
  }

  # centered
  odd <- odd_f(N)
  odd.nc <- odd +  baseline.params$p.odd
  severity <- severity_f(N)
  severity.nc <- severity + baseline.params$m.severity
  priormed <- priormed_f(N)
  priormed.nc <- priormed + baseline.params$p.priormed
  race <- race_f(N)
  race.nc <- race + baseline.params$p.race

  X.c <- data.frame(odd, severity, priormed, race)

  Y0 <- linMult(Y0.coef) + rnorm(N)

  #Y02 <- cbind(1, odd, severity) %*% c(2, -0.2, -0.3)

  # First-stage -------------------------------------------------------------
  #Med: -1, BMOD: 1
  A1_f <- function(N){2 * rbinom(N, 1, 0.5) - 1}
  A1 <- A1_f(N)

  ## Unknown
  U_f <- function(N){rnorm(N, U.params$mu, U.params$sd)} # common cause of R, Adh, Y1, Y2
  U <- U_f(N)

  ## Response OC: Med (-1) positively associated with Response, priorMed + association
  lo.R <- linMult(R.coef)

  lo.R2 <- cbind(1, A1, priormed*A1, U) %*% c(Int = -0.7, A1 = -0.2, "pm:A1" = -0.4, U = 0.5)  # log-odds scale
  p.R <- plogis(lo.R)
  R.nc <- rbinom(N, 1, p.R)
  R <- R.nc - p.R

  # E[R | A1]
  a1 <- c(-1, 1)
  ER.A1 <- sapply(a1, sampleLogisticMean, named.coefs = R.coef)

  ## Adherence OC: For those on PriorMed, Med(-1) positively associated with Adherence
  lo.adherence <- linMult(adherence.coef)

  #lo.adherence2 <- cbind(1, priormed * A1, U) %*% c(Int = -0.2, "pm:A1" = -0.3, U = 0.2)  # log-odds scale
  p.adherence <- plogis(lo.adherence)
  adherence.nc <- rbinom(N, 1, p.adherence)
  adherence <- adherence.nc - p.adherence

  # E[Adh | A1]
  EAdh.A1 <- sapply(a1, sampleLogisticMean, named.coefs = adherence.coef)

  ## Time to Non-response (Event-time) OC: No association for now
  ET <- sample.int(8, N, replace = TRUE)  # uniform; for now
  ET[R.nc == 1] <- NA

  ## First-stage Outcome OC: Small change from baseline, Med (-1) is better initially
  Y1 <- linMult(Y1.coef) + rnorm(N)


  # Second-stage ------------------------------------------------------------

  #AUG: -1, INT: 1

  A2_f <- function(N){2 * rbinom(N, 1, 0.5) - 1}
  A2 <- A2_f(N)
  A2[R.nc == 1] <- 0  # Fill NA at end


  # EOS Outcome -----------------------------------------------------------------

  ## Baseline
  Y2 <- linMult(Y2.baseline)
  #Y22 <- cbind(1, odd, severity, priormed, race) %*% c(3, -0.3, -0.4, 0, 0.5)

  ## First-stage OC: BMOD is marginally better in the long run. For those on PriorMed, Med(-1) is better. For those on PriorMed and
  # Responded, Med(-1) is MUCH better.
  tx1 <- linMult(Y2.tx1)  # small main effect!
  Y2 <- Y2 + tx1

  ## Intermediate Associations
  n1 <- linMult(Y2.n1)
  Y2 <- Y2 + n1

  ## Second-stage: For non-responders OC: AUG(-1) is marginally better. The effect of ADD(-1) is neg for those that
  ## started with Med(-1). For those Adherent, INT is MUCH better.
  #tx2 <- (1 - R) * cbind(A2, A1 * A2, adherence * A2) %*% c(-0.3, -0.2, 1.2)
  tx2 <- (1 - R.nc) * linMult(Y2.tx2)
  Y2 <- Y2 + tx2

  ## Error
  Y2 <- Y2 + rnorm(N, 0, 0.5) # determine what error is needed for sig effects

  A2[R.nc == 1] <- NA
  ET[R.nc == 1] <- NA
  browser()
  ## Marginal Model
  marginalMeanModel <- function(a1, a2) {
    A1 <-  a1
    A2 <-  a2
    adherence <-  EAdh.A1[-a1]
    R <- ER.A1[-a1]
    odd <-  0
    severity <-  0
    priormed <- 0
    race <- 0
    linMult(Y2.baseline) + linMult(Y2.tx1) + (1 - ER.A1[-a1]) * linMult(Y2.tx2)
  }

  data <- data.frame(ID = 1:N, odd = odd.nc, severity = severity.nc,
                     priormed = priormed.nc, race = race.nc, Y0,
                     A1, R, ET, adherence = adherence.nc, Y1, A2, Y2)
  return(list(data))

}



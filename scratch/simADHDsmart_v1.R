
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

simADHDsmart <- function(N = 100) {

  # Baseline Parameters (from old synthetic data) ----------------------------
  p.odd <- 0.353  # binary, O11
  mean.severity <- 0  # continuous; O12 mean centered and scaled
  sd.severity <- 1
  p.priormed <- 0.313  # binary, O13
  p.race <- 0.807  # binary, O14


  # Baseline Covariates -----------------------------------------------------

  odd <- rbinom(N, 1, p.odd)
  odd.c <- odd - p.odd  # centered
  severity <- rnorm(N)
  severity.c <- severity
  priormed <- rbinom(N, 1, p.priormed)
  priormed.c <- priormed - p.priormed
  race <- rbinom(N, 1, p.race)
  race.c <- race - p.race

  X <- data.frame(odd.c, severity.c, priormed.c, race.c)

  Y0 <- cbind(1, odd.c, severity.c) %*% c(2, -0.2, -0.3) + rnorm(N)


  # First-stage -------------------------------------------------------------

  #Med: -1, BMOD: 1
  A1 <- 2 * rbinom(N, 1, 0.5) - 1

  ## Unknown
  U <- rnorm(N) # common cause of R, Adh, Y1, Y2

  ## Response OC: Med (-1) positively associated with Response, priorMed + association
  coef.R <- c(Int = -0.7, A1 = -0.2, "pm:A1" = -0.4, U = 0.5)
  lo.R <- cbind(1, A1, priormed.c*A1, U) %*% coef.R  # log-odds scale
  p.R <- plogis(lo.R)
  R <- rbinom(N, 1, p.R)
  R.c <- R - p.R

  # Marginal Response E[R | A1]
  # use sampling

  marginalR <- function(N = 1e-6, C = coef.R) {
    priormed <- rbinom(N, 1, p.priormed)
    priormed.c <- priormed - p.priormed
    U <- rnorm(N)
    out <- rep(0,2)
    a1 <- c(-1,1)
    for (i in 1:2) {
      out[i] <- mean(plogis(cbind(1, a1[i], priormed.c*a1[i], U) %*% t(C)))
      }
  }

  ## Adherence OC: For those on PriorMed, Med(-1) positively associated with Adherence
  lo.adherence <- cbind(1, priormed.c * A1, U) %*% c(Int = -0.2, "pm:A1" = -0.3, U = 0.2)  # log-odds scale
  p.adherence <- plogis(lo.adherence)
  adherence <- rbinom(N, 1, p.adherence)
  adherence.c <- adherence - p.adherence

  ## Time to Non-response (Event-time) OC: No association for now
  ET <- sample.int(8, N, replace = TRUE)  # uniform; for now
  ET[R == 1] <- NA

  ## First-stage Outcome OC: Small change from baseline, Med (-1) is better initially
  Y1 <- cbind(1, A1, U) %*% c(Int = 2.5, A1 = -0.3, U = 0.4) + rnorm(N)


  # Second-stage ------------------------------------------------------------

  #AUG: -1, INT: 1

  A2 <- 2 * rbinom(N, 1, 0.5) - 1
  A2[R == 1] <- 0  # Fill NA at end


  # EOS Outcome -----------------------------------------------------------------

  ## Baseline
  Y2 <- cbind(1, odd.c, severity.c, priormed.c, race.c) %*% c(3, -0.3, -0.4, 0, 0.5)

  ## First-stage OC: BMOD is marginally better in the long run. For those on PriorMed, Med(-1) is better. For those on PriorMed and
  # Responded, Med(-1) is MUCH better.
  tx1 <- cbind(A1, priormed.c * A1, priormed.c * R.c * A1) %*% c(0.2, -0.6, -1.2)  # small main effect!
  Y2 <- Y2 + tx1

  ## Intermediate Associations
  Y2 <- Y2 + 0.4 * R.c + 0.4 * U

  ## Second-stage: For non-responders OC: AUG(-1) is marginally better. The effect of ADD(-1) is neg for those that
  ## started with Med(-1). For those Adherent, INT is MUCH better.
  tx2 <- (1 - R) * cbind(A2, A1 * A2, adherence.c * A2) %*% c(-0.3, -0.2, 1.2)
  Y2 <- Y2 + tx2

  ## Error
  Y2 <- Y2 + rnorm(N, 0, 0.5) # determine what error is needed for sig effects

  A2[R == 1] <- NA
  ET[R == 1] <- NA

  data <- data.frame(ID = 1:N, odd, severity, priormed, race, Y0,
                     A1, R, ET, adherence, Y1, A2, Y2)
  return(data)

}


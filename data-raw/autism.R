# Generate realization

set.seed(2023-01-21)
autism <- simAUTISMsmart(N = 200,
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
                         sigma = 8)

readr::write_csv(autism, here::here("data-raw","autism-simulated-2023.csv"))

usethis::use_data(autism, overwrite = TRUE)

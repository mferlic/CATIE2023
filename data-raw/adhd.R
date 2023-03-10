# Generate realization

set.seed(1234)
adhd <- simADHDsmart(N = 150,
                     baseline.params = list(p.odd = 0.4, m.severity = 5, p.priormed = 0.3, p.race = 0.8),
                     Y0.coef = c(2, odd = -0.4, severity = -0.1, U = 0.1),
                     U.params = c(mu = 0, sd = 1),
                     R.coef = c(-0.4, A1 = -0.3, "priormed:A1" = -0.2, U = 0.2, Y0 = 0.1), # probit model
                     adherence.coef = c(-0.1, "priormed:A1" = -0.2, U = 0.2), # probit model
                     NRtime.coef = NULL,
                     Y1.coef = c(2.5, A1 = -0.3, U = 0.5, Y0 = 0),
                     Y2.baseline = c(3, odd = -0.5, severity = -0.1, priormed = 0, race = 0.4),
                     Y2.tx1 = c(A1 = 0.4, "priormed:A1" = -2),
                     Y2.n1 = c(R.resid = 1, adherence.resid = 2, U = 0.2, Y1.resid = 1.2),
                     Y2.tx2 = c(A2 = -0.3, "A1:A2" = 0.1, "adherence:A2" = 2),
                     sigma = 0.4,
                     rho = 0.8)

readr::write_csv(adhd, here::here("data-raw","adhd-simulated-2023.csv"))

usethis::use_data(adhd, overwrite = TRUE)

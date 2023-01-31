# Generate realization

set.seed(2023-01-21)
autObj <- simAUTISMsmart()
autism <- autObj$data
readr::write_csv(autism, here::here("data-raw","autism-simulated-2023.csv"))
usethis::use_data(autism, overwrite = TRUE)

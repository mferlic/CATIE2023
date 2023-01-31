# Generate realization

set.seed(2023-01-01)
adhdObj <- simADHDsmart()
adhd <- adhdObj$data
readr::write_csv(adhd, here::here("data-raw","adhd-simulated-2023.csv"))
usethis::use_data(adhd, overwrite = TRUE)

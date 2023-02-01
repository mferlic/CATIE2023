# Generate realization

set.seed(2023-01-01)
adhdObj <- simADHDsmart()
adhd <- adhdObj$data
adhdDTRmeans <- adhdObj$DTRmeans

readr::write_csv(adhd, here::here("data-raw","adhd-simulated-2023.csv"))

usethis::use_data(adhd, overwrite = TRUE)
usethis::use_data(adhdDTRmeans, overwrite = TRUE)

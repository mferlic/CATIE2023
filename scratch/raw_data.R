
set.seed(2023-01-01)
adhdObj <- simADHDsmart()
adhd <- adhdObj$data
use_data(adhd, overwrite = TRUE)

set.seed(2023-01-21)
autObj <- simAUTISMsmart()
autism <- autObj$data
use_data(autism, overwrite = TRUE)



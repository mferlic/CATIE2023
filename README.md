
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CATIE2023

<!-- badges: start -->
<!-- badges: end -->

The goal of CATIE2023 is to bundle the data sets and code used for the
CATIE 2023 demo and practicum days.

## Installation

You can install the development version of CATIE2023 from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mferlic/CATIE2023")
```

## Example

This is how you access the internal data sets

``` r
library(CATIE2023)

adhd <- CATIE2023::adhd # data.frame

# read documentation
?adhd
```

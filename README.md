
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CATIE2023

<!-- badges: start -->

[![R-CMD-check](https://github.com/mferlic/CATIE2023/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mferlic/CATIE2023/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of CATIE2023 is to bundle the datasets and code used for the
CATIE 2023 demo and practicum days.

## Installation

You can install the development version of CATIE2023 from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mferlic/CATIE2023")
```

## Example

This is a basic example which shows you how to access the ADHD dataset

``` r
library(CATIE2023)

adhd <- CATIE2023::adhd
```

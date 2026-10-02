
# aprlab04

<!-- badges: start -->
[![R-CMD-check](https://github.com/totoposeb/aprlab04/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/totoposeb/aprlab04/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This repository contains an implementation of the multiple linear regression algorithm in R, using QR decomposition and a reference class (RC): `linreg()`

## Installation

You can install the development version of aprlab04 from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("totoposeb/aprlab04", build_vignettes = TRUE)
```

## Example

This is a basic example which shows you how to fit a linear regression model with `linreg`:

``` r
library(aprlab04)
m <- linreg$new(Petal.Length ~ Sepal.Width + Sepal.Length, data = iris)
m$print()
m$coef()
m$plot()
```


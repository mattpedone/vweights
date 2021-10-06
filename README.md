
# vweights

<!-- badges: start -->
<!-- badges: end -->

The goal of `vweights` is to compute non-negative weights for cohesion function used to define a law for random partition, that is of Gibbs type.

The implementation has been done in `C++` through the use of `Rcpp` and `RcppArmadillo`.

Authors: Matteo Pedone, Raffaele Argiento.

Maintainer: Matteo Pedone.

## Installation

You can install the development version from [GitHub](https://CRAN.R-project.org), using the **devtools** package with:

``` r
devtools::install_github("mattpedone/vweights")
```

**NOTE** that this package depends on [gsl](https://www.gnu.org/software/gsl/). It has only been tested on a PC running Ubuntu 20.04.2 LTS.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(vweights)
vweights::computev(5, .1, 10.5)
```


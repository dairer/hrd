
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hrd

<!-- badges: start -->
<!-- badges: end -->

The goal of hrd is to …

## Installation

You can install the development version of hrd like so:

``` r
install.packages('devtools')
devtools::install_github("dairer/hrd")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
## sample data from a bivariate Hüsler-Reiss copula
library(hrd)
rhr(n=10, lambda = 1)
#>                u         v
#>  [1,] 0.47175837 0.2948507
#>  [2,] 0.51771378 0.5014396
#>  [3,] 0.71541700 0.4621651
#>  [4,] 0.72278432 0.5397291
#>  [5,] 0.11767617 0.3935401
#>  [6,] 0.77555017 0.8174187
#>  [7,] 0.57016371 0.7556169
#>  [8,] 0.06499856 0.3435320
#>  [9,] 0.24603109 0.7828829
#> [10,] 0.08167499 0.2832719
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
# summary(cars)
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

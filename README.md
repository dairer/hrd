
# hrd

<!-- badges: start -->
<!-- badges: end -->

`hrd` is a collection of functions to help implement, develop and fit
Hüsler-Reiss distribution based models.

**Introduction**

> -   [Getting set up](#installation)

**Wokring with bivariate Hüsler-Reiss copula**

> -   [Sampling data](#example)
> -   [Fitting
>     (MLE)](#esimating-the-dependence-parameter-of-a-hüsler-reiss-copula)

## Installation

You can install the development version of hrd like so:

``` r
install.packages('devtools')
devtools::install_github("dairer/hrd")
```

``` r
library(hrd)
library(ggplot2) 
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.1
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1
    ## ✔ purrr   1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

## Example

### Sampling bivariate Hüsler-Reiss copula

``` r
# sample data from a bivariate Hüsler-Reiss copula
data = rhr(n=1000, lambda = 1.5)

# visualise samples
data %>%
  as.data.frame() %>%
  ggplot()+
  geom_point(aes(u,v), shape = 4)+
  theme_minimal(12)
```

![](README_files/figure-gfm/example1-1.png)<!-- -->

### Esimating the dependence parameter of a Hüsler-Reiss copula

``` r
fit_hrc(u = data[,1], v = data[,2], initial_est = 1)
```

    ## $estimate
    ## [1] 1.465131
    ## 
    ## $ci
    ## [1] 1.373155 1.557108

### Esimating generalised pareto marings and dependence parameter of a Hüsler-Reiss copula jointly

``` r
# transform margins to be GP distributed
x = qgp(data[,1], scale = 2, shape = -0.1)
y = qgp(data[,2], scale = 1.3, shape = -0.1)

# visualise the multivariate distribution
data.frame(x,y) %>%
  ggplot()+
  geom_point(aes(x,y), shape = 4)+
  theme_minimal()
```

![](README_files/figure-gfm/example3-1.png)<!-- -->

``` r
model_fit = fit_hrc_gp(x,y,initial_est = c(1,0,1,0,1))

# visualise estimates
data.frame(actual = c(2, -0.1, 1.3, -0.1, 1.5),
       estimates = model_fit$estimate,
       lower = model_fit$ci[,1],
       upper = model_fit$ci[,2],
       params = c('σ1', 'ξ1', 'σ2', 'ξ2', 'λ')) %>%
  ggplot()+
  geom_point(aes(params, estimates))+
  geom_point(aes(params, actual),position = position_nudge(x = -0.075), shape = 4, col = 'red')+
  geom_segment(aes(x = params, xend = params, y = lower, yend = upper))+
  labs(x = "Parameters",
       y = "Parameter estimates")+
  theme_minimal(12)
```

![](README_files/figure-gfm/example3-2.png)<!-- -->

In each case the red ‘x’ is the true value. The black point along with
the segment are the point estimates and 95% confidence intervals. λ is
the dependence parameter of the Hüsler-Reiss copula. σ1 and ξ1 are the
scale and shape parameters of GPD margin 1. Equivalent for margin 2.

### We can compare estimates from margin-dependence indpendent analysis (in magenta) to joint modelling (in black).

``` r
dependence_mod = fit_hrc(data[,1], data[,2], initial_est = 1)

marg1_mod = fit_gpd(x, initial_est = c(1,0))
marg2_mod = fit_gpd(y, initial_est = c(1,0))

data.frame(actual = c(2, -0.1, 1.3, -0.1, 1.5),
       estimates_joint = model_fit$estimate,
       lower_joint = model_fit$ci[,1],
       upper_joint = model_fit$ci[,2],
       estimates_indep = c(marg1_mod$estimate, marg2_mod$estimate, dependence_mod$estimate),
       lower_indep = c(marg1_mod$ci[,1],marg2_mod$ci[,1],dependence_mod$ci[1]),
       upper_indep = c(marg1_mod$ci[,2],marg2_mod$ci[,2],dependence_mod$ci[2]),
       params = c('σ1', 'ξ1', 'σ2', 'ξ2', 'λ')) %>%
  ggplot()+
  geom_point(aes(params, actual),position = position_nudge(x = -0.075), shape = 4, col = 'red')+
  geom_point(aes(params, estimates_joint))+
  geom_segment(aes(x = params, xend = params, y = lower_joint, yend = upper_joint))+
  geom_point(aes(params, estimates_indep),position = position_nudge(x = 0.075), col = "magenta")+
  geom_segment(aes(x = params, xend = params, y = lower_indep, yend = upper_indep),position = position_nudge(x = 0.075), col = "magenta")+
  labs(x = "Parameters",
       y = "Parameter estimates")+
  theme_minimal(12)
```

![](README_files/figure-gfm/example4-1.png)<!-- -->

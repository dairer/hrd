% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bivariate_hr.R
\name{fit_gpd_bay}
\alias{fit_gpd_bay}
\title{Fitting a GPD in a bayesian framework}
\usage{
fit_gpd_bay(
  x,
  scale_prior_param = 1,
  shape_prior_lower = -0.5,
  shape_prior_upper = 0.5,
  chains = 2,
  iter = 500,
  warmup = floor(iter/2),
  cores = 2,
  thin = 1
)
}
\arguments{
\item{x}{(numeric) data}

\item{scale_prior_param}{(numeric) rate of exp prior for scale parameter}

\item{shape_prior_lower}{(numeric) lower bound for uniform prior for shape parameter}

\item{shape_prior_upper}{(numeric) upper bound for uniform prior for shape parameter}

\item{chains}{(numeric) see ?rstan::sampling}

\item{iter}{(numeric) see ?rstan::sampling}

\item{warmup}{(numeric) see ?rstan::sampling}

\item{cores}{(numeric) see ?rstan::sampling}

\item{thin}{(numeric) see ?rstan::sampling}

\item{v}{(numeric) Uniform margin of copula}
}
\value{
(list(estimate (numeric), post (numeric), fitted_model (stan model))) List of 3 elements,
first element is the estimate of the GPD parameters, calculated as the mean of the posterior distribution.
Second element is the posterior distribution samples. Third element is the full fitted stan model.
}
\description{
Fitting a GPD in a bayesian framework
}
\examples{
x = rgp(1000, scale = 2, shape = -0.1)
fit_gpd_bay(x)
}

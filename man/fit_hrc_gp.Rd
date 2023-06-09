% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bivariate_hr.R
\name{fit_hrc_gp}
\alias{fit_hrc_gp}
\title{Jointly fitting a bi-variate Hüsler-Reiss copula and generlaised Pareto margins}
\usage{
fit_hrc_gp(x, y, initial_est)
}
\arguments{
\item{x}{(numeric) margin 1, assumed to be generalised Pareto distributed}

\item{y}{(numeric) margin 2, assumed to be generalised Pareto distributed}

\item{initial_est}{(numeric) Vector of initial parameter estimates}
}
\value{
(list(estimate (numeric), ci (matrix))) List of 2 elements,
first element is the estimates of the generalised Pareto distribution parameters and Hüsler Reiss copula parameter,
second element is the estimated 95% confidence interval from the Hessian, the i-th row corresponds to the i-th estimate.
}
\description{
Jointly fitting a bi-variate Hüsler-Reiss copula and generlaised Pareto margins
}
\examples{
# sample from a Hüsler-Reiss copula
copula_sample = rhr(1000, 2)

# Transform margins to be GPD
margin_1 = qgp(copula_sample[,1], scale = 2, shape = -0.1)
margin_2 = qgp(copula_sample[,2], scale = 1.3, shape = -0.1)

# fit margins and copula jointly
fit_hrc_gp(margin_1, margin_2, initial_est = c(1,0.1,1,-0.12,1))
}

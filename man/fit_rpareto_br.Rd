% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brown_resnick.R
\name{fit_rpareto_br}
\alias{fit_rpareto_br}
\title{Estimate variogram parameters for Brown-Resnick process with changing or constant dimensions.}
\usage{
fit_rpareto_br(dt, lcs, vr, conditioned.site = 1, initial_est = c(1, 1))
}
\arguments{
\item{dt}{(list) Observations of extreme process, with standarised margins.}

\item{lcs}{(list) List of matrices that describe locations of observations in corresponding position in list dt.}

\item{vr}{(numeric) Variogram function, evaluates with "pars"}

\item{conditioned.site}{(numeric) Site that is "conditioned" on.}

\item{initial_est}{(numeric) Starting point of variogram parameter optimisation}
}
\value{
(numeric) Log-liklihood
}
\description{
Estimate variogram parameters for Brown-Resnick process with changing or constant dimensions.
}

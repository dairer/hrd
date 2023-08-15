#' Density function of the generalised Pareto distribution
#'
#' @param x (numeric) data
#' @param scale (numeric)
#' @param shape (numeric)
#' @param log (logical) If true, returns log
#'
#' @return (numeric)
#' @export
#'
#' @examples
#' dgp(x = 0.5, scale = 1, shape = 0)
dgp = function(x , scale, shape, log = F){
  # --- error control
  if(scale < 0)
    stop("Unsupported generalised Pareto distribution")

  if((shape < 0) & any(x > (-scale/shape)))
    stop("Unsupported generalised Pareto distribution")

  if(shape == 0) res = exp(-x/scale) * (1/scale)
  else res = 1/scale*(1+(shape*x)/scale)^(-1/shape - 1)

  if(log) log(res)
  else res
}

#' Cumulative distribution function of the generalised Pareto distribution
#'
#' @param x (numeric) data
#' @param scale (numeric)
#' @param shape (numeric)
#'
#' @return (numeric)
#' @export
#'
#' @examples
#' pgp(x = 0.5, scale = 1, shape = 0)
pgp = function(x, scale, shape){
  # --- error control
  if(scale < 0)
    stop("Unsupported generalised Pareto distribution")

  if((shape < 0) & any(x > (-scale/shape)))
    stop("Unsupported generalised Pareto distribution")

  if(shape == 0) 1-exp(-x/scale)
  else 1-(1+((shape*x)/scale))^(-1/shape)
}

#' Quantile function of the generalised Pareto distribution (inverse cumulative distribution function)
#'
#' @param x (numeric) data
#' @param scale (numeric)
#' @param shape (numeric)
#'
#' @return (numeric)
#' @export
#'
#' @examples
#' pgp(x = 0.5, scale = 1, shape = 0)
qgp = function(x, scale, shape){
  # --- error control
  if(scale < 0)
    stop("Unsupported generalised Pareto distribution")

  if((shape < 0) & any(x > (-scale/shape)))
    stop("Unsupported generalised Pareto distribution")

  if (shape == 0)  -scale*log(1-x)
  else scale * ((1-x)^(-shape) - 1)/shape
}

#' Generate random samples from the generalised Pareto distribution
#'
#' @param n (numeric) Number of samples
#' @param scale (numeric)
#' @param shape (numeric)
#'
#' @return (numeric)
#' @export
#'
#' @examples
#' rgp(n = 100, scale = 1, shape = 0)
rgp = function(n, scale, shape) stats::runif(n) %>% qgp(scale, shape)

# negative log likelihood function of the generalised Pareto distribution
gp_ngll = function(par, x){
  if(par[1]<0) return(100^100)
  if((par[2] < 0) & any(x > (-par[1]/par[2]))) return(100^100)
  -sum(dgp(x, scale = par[1], shape = par[2], log = T))
}

#' Fitting a generalised Pareto distribution
#'
#' @param x (numeric) Vector of observations
#' @param initial_est (numeric) Initial estimates of distribution parameter (scale, shape)
#'
#' @return (list(estimate (numeric), ci (matrix))) List of 2 elements,
#' first element is the estimates of the generalised Pareto distribution parameters,
#' second element is the estimated 95% confidence interval from the Hessian, the i-th row corresponds to the i-th estimate.
#' @export
#'
#' @examples
#' x = rgp(500, scale = 1, shape = -0.3)
#' fit_gpd(x)
fit_gpd = function(x, initial_est = c(1,0)){
  this_fit = stats::optim(fn=gp_ngll,
                   par = initial_est,
                   x = x,
                   hessian = T)

  if(det(this_fit$hessian)==0){
    cat("Could not compute confidence interval as det(hessian) = 0. Try different initial_est.")
    list(Estimate = this_fit$par)
  }else{
    # estimate confidence interval from Hessian matric
    this_se = calc_se(this_fit$hessian)
    # return estimate and CI
    list(estimate = this_fit$par,
         ci = matrix(c(this_fit$par-this_se, this_fit$par+this_se), nrow = 2))
  }
}





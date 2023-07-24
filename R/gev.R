#' Density function of the generalised extreme value distribution
#'
#' @param x (numeric) data
#' @param loc (numeric)
#' @param scale (numeric)
#' @param shape (numeric)
#' @param log (logical) If true, returns log
#'
#' @return (numeric)
#' @export
#'
#' @examples
#' dgev(x = 0.5, loc = 0, scale = 1, shape = 0)
dgev = function(x , loc, scale, shape, log = F){
  # --- error control
  if(scale <= 0)
    stop("Unsupported geneneralised extreme value distribution")

  if((shape < 0) & any(x > (loc - (scale/shape))))
    stop("Unsupported geneneralised extreme value distribution")

  if((shape > 0) & any(x < (loc - (scale/shape))))
    stop("Unsupported geneneralised extreme value distribution")

  tm = rep(NA, length(x))
  if(shape != 0){
    tm = (1 + shape*((x - loc)/scale))^(-1/shape)
  }else{
    tm = exp(-(x - loc)/scale)
  }

  res =  (1/scale)*(tm^(shape+1))*exp(-tm)

  if(log) log(res)
  else res
}




#' Cumulative distribution function of the generalised extreme value distribution
#'
#' @param x (numeric) data
#' @param loc (numeric)
#' @param scale (numeric)
#' @param shape (numeric)
#'
#' @return (numeric)
#' @export
#'
#' @examples
#' pgev(x = 0.5, loc = 0, scale = 1, shape = 0)
pgev = function(x, loc, scale, shape){
  # --- error control
  if(scale <= 0)
    stop("Unsupported geneneralised extreme value distribution")

  if((shape < 0) & any(x > (loc - (scale/shape))))
    stop("Unsupported geneneralised extreme value distribution")

  if((shape > 0) & any(x < (loc - (scale/shape))))
    stop("Unsupported geneneralised extreme value distribution")

  tm = rep(NA, length(x))
  if(shape != 0){
    tm = pmax((1 + shape*((x - loc)/scale)), 0)^(-1/shape)
  }else{
    tm = exp(-(x - loc)/scale)
  }

  exp(-tm)
}


#' Quantile function of the generalised extreme value distribution (inverse cumulative distribution function)
#'
#' @param x (numeric) data
#' @param loc (numeric)
#' @param scale (numeric)
#' @param shape (numeric)
#'
#' @return (numeric)
#' @export
#'
#' @examples
#' qgev(x = 0.5, loc = 1, scale = 1, shape = 0)
qgev = function(x, loc, scale, shape){
  # --- error control
  if(scale <= 0)
    stop("Unsupported geneneralised extreme value distribution")

  if((shape < 0) & any(x > (loc - (scale/shape))))
    stop("Unsupported geneneralised extreme value distribution")

  if((shape > 0) & any(x < (loc - (scale/shape))))
    stop("Unsupported geneneralised extreme value distribution")

  if (shape == 0) return(loc - scale * log(-log(x)))
  else return(loc+scale*((-log(x))^(-shape)-1)/shape)
}




#' Generate random samples from the generalised extreme value distribution
#'
#' @param n (numeric) Number of samples
#' @param loc (numeric)
#' @param scale (numeric)
#' @param shape (numeric)
#'
#' @return (numeric)
#' @export
#'
#' @examples
#' rgev(n = 100, loc = 1, scale = 1, shape = 0)
rgev = function(n, loc, scale, shape) runif(n) %>% qgev(loc, scale, shape)

# negative log likelihood function of the generalised extreme value distribution
gev_ngll = function(par, x){

  loc = par[1]
  scale = par[2]
  shape = par[3]

    if(scale <= 0) return(100^100)
  if((shape < 0) & any(x > (loc - (scale/shape)))) return(100^100)
  if((shape > 0) & any(x < (loc - (scale/shape)))) return(100^100)

  -sum(dgev(x, loc = loc, scale =scale, shape = shape, log = T))
}

#' Fitting a generalised extreme value distribution
#'
#' @param x (numeric) Vector of observations
#' @param initial_est (numeric) Initial estimates of distribution parameter (scale, shape)
#'
#' @return (list(estimate (numeric), ci (matrix))) List of 2 elements,
#' first element is the estimates of the generalised extreme value distribution parameters,
#' second element is the estimated 95% confidence interval from the Hessian, the i-th row corresponds to the i-th estimate.
#' @export
#'
#' @examples
#' x = rgev(500, loc = 5, scale = 1, shape = -0.2)
#' fit_gevd(x)
fit_gevd = function(x, initial_est = c(0,1,0)){
  this_fit = optim(fn=gev_ngll,
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

#' Cumulative distribution function of the bi-variate Hüsler-Reiss copula
#'
#' @param u (numeric) Uniform margin of copula
#' @param v (numeric) Uniform margin of copula
#' @param lambda (numeric) Dependence parameter of the Hüsler-Reiss copula
#' @param log (logical) If true returns log-Likelihood, default is False
#'
#' @return (numeric)
#'
#' @export
#'
#' @examples
#' dhr(u = c(0.1, 0.2, 0.4), v = c(0.6, 0.7, 0.2), lambda = 1)
phr = function(u,v, theta){
  exp(log(u) * pnorm((1/theta) + (theta/2)*log(log(u)/log(v))) +
        log(v)*pnorm((1/theta) + (theta/2)*log(log(v)/log(u))))
}

#' Density function of the bi-variate Hüsler-Reiss copula
#'
#' @param u (numeric) Uniform margin of copula
#' @param v (numeric) Uniform margin of copula
#' @param lambda (numeric) Dependence parameter of the Hüsler-Reiss copula
#' @param log (logical) If true returns log
#'
#' @return (numeric)
#' @export
#'
#' @examples
#' dhr(u = c(0.1, 0.2, 0.4), v = c(0.6, 0.7, 0.2), lambda = 1)
dhr = function(u, v, lambda, log = F){
  z = log(log(u)/log(v))
  a = 1/lambda + 0.5*lambda*z
  # note that log(x) = -log(1/x)
  likelihood = 1/(u*v) * phr(u,v,lambda) *
    (pnorm(1/lambda - 0.5*lambda*z) *pnorm(a) + 0.5*lambda * -1/log(v)*dnorm(a))

  if(log) log(likelihood)
  else likelihood
}

#' Generate random samples from the bi-variate Hüsler-Reiss copula
#'
#' @param n (numeric) Number of random samples
#' @param lambda (numeric) Dependence parameter of the Hüsler-Reiss copula
#'
#' @return (Matrix) Random samples from bi-variate Hüsler-Reiss copula
#' @export
#'
#' @examples
#' rhr(n = 10, lambda = 1)
rhr = function(n, lambda){

  u = runif(n)
  solved_v <- numeric(n) # for storing simulated r.v.

  # dC/du, from Joe, 1997
  conditional_density <- function(u, v, y) {
    phr(u,v,lambda) * 1/u * pnorm(1/lambda + 0.5*lambda*log(log(u)/log(v))) - y
  }

  # find root of this fn (finding inverse)
  # see https://stackoverflow.com/questions/10081479/solving-for-the-inverse-of-a-function-in-r
  solved_v = seq(n) %>%
    purrr::map(~{
      uniroot(conditional_density,
              interval = c(1e-20, 1-(1e-7)),
              y=runif(1),
              u=u[.x])$root
    }) %>%
    unlist

  cbind(u, v=solved_v)
}

calc_se = function(hess){
  se = hess %>% solve() %>% diag %>% sqrt()
  1.96*se
}

hrc_ll = function(pars, u, v){
  LL = -sum(dhr(u,v,lambda = pars, log = T))
  if(is.infinite(LL)) return(100^100)
  LL
}

#' Fitting a bi-variate Hüsler-Reiss copula
#'
#' @param u (numeric) Uniform margin of copula
#' @param v (numeric) Uniform margin of copula
#' @param lambda (numeric) Dependence parameter of the Hüsler-Reiss copula
#' @param log (logical) If true returns log-Likelihood, default is False
#'
#' @return (numeric)
#'
#' @export
#'
#' @examples
#' dat = rhr(1000, 1.2)
#' fit_hrc(dat[,1], dat[,2],1)
fit_hrc = function(u, v, initial_est){

  # --- Error handling
  # check variables are same length
  if(length(u) != length(v)){
    stop("Uniform variables must be the same length")
  }

  # check if values are in the range [0,1]
  unif_range_check = range(u,v)
  if(any(unif_range_check < 0) | any(unif_range_check>1)){
    stop("Uniform variables must be in the range [0,1]")
  }

  this_fit = optim(fn=hrc_ll,
                   par = initial_est,
                   u = u,v = v,
                   method = 'Brent',
                   lower = 0.0000001,
                   upper = 15,
                   hessian = T)

  if(det(this_fit$hessian)==0){
    cat("Could not compute confidence interval as det(hessian) = 0. Try different initial_est.")
    list(Estimate = this_fit$par)
  }else{

    # estimate confidence interval from Hessian matric
    this_se = calc_se(this_fit$hessian)
    # return estimate and CI
    list(Estimate = this_fit$par,
         ci = c(this_fit$par-this_se, this_fit$par+this_se))
  }
}






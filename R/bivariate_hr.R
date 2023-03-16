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
#' @param log (logical) If true, returns log
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
  se = hess %>% solve(tol = 1e-1000) %>% diag %>% sqrt()
  1.96*se
}

hrc_ngll = function(pars, u, v){
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
#' @return (list(estimate (numeric), ci (numeric))) List of 2 elements,
#' first element is the estimate of the Hüsler-Reiss dependence parameter,
#' second element is the estimated 95% confidence interval from the Hessian.
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

  this_fit = optim(fn=hrc_ngll,
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
    list(estimate = this_fit$par,
         ci = c(this_fit$par-this_se, this_fit$par+this_se))
  }
}

# negative log likelihood of copula and margins
hrc_gp_nll = function(pars, x, y){

  if(any(pars[1] <= 0)) return(100^100)
  if(pars[2]<0){
    if(any(pars[1] > -1/pars[2])) return(100^100)
    if(any((1+pars[2]*x/pars[1])< 0)) return(100^100)
  }

  if(any(pars[3] <= 0)) return(100^100)
  if(pars[4]<0){
    if(any(pars[3] > -1/pars[4])) return(100^100)
    if(any((1+pars[4]*y/pars[3])< 0)) return(100^100)
  }

  if(pars[5]<=0 | pars[5]>10) return(100^100)

  # likelihood of margin 1
  llm1 = dgp(x, scale = pars[1], shape = pars[2])

  # likelihood of margin 2
  llm2 = dgp(y, scale = pars[3], shape = pars[4])

  # likelihood of copula
  unif_1 = pgp(x, scale = pars[1], shape = pars[2])
  unif_2 = pgp(y, scale = pars[3], shape = pars[4])
  llc = dhr(unif_1, unif_2, lambda = pars[5])

  # full = likelihood
  LL =  -sum(log(llm1*llm2*llc))

  if(is.infinite(LL)) return(100^100)
  if(is.na(LL)) return(100^100)
  LL
}

#' Jointly fitting a bi-variate Hüsler-Reiss copula and generlaised Pareto margins
#'
#' @param x (numeric) margin 1, assumed to be generalised Pareto distributed
#' @param y (numeric) margin 2, assumed to be generalised Pareto distributed
#' @param initial_est (numeric) Vector of initial parameter estimates
#'
#' @return (list(estimate (numeric), ci (matrix))) List of 2 elements,
#' first element is the estimates of the generalised Pareto distribution parameters and Hüsler Reiss copula parameter,
#' second element is the estimated 95% confidence interval from the Hessian, the i-th row corresponds to the i-th estimate.
#' @export
#'
#' @examples
#' # sample from a Hüsler-Reiss copula
#' copula_sample = rhr(1000, 2)
#'
#' # Transform margins to be GPD
#' margin_1 = qgp(copula_sampe[,1], scale = 2, shape = -0.1)
#' margin_2 = qgp(copula_sampe[,2], scale = 1.3, shape = -0.1)
#'
#' # fit margins and copula jointly
#' fit_hrc_gp(margin_1, margin_2, initial_est = c(1,0.1,1,-0.12,1))
fit_hrc_gp = function(x, y, initial_est){
  # --- Error handling
  # check variables are same length
  if(length(x) != length(y)){
    stop("Variables must be the same length")
  }


  this_fit = optim(fn=hrc_gp_nll,
                   par = initial_est,
                   x = x, y = y,
                   hessian = T)

  if(det(this_fit$hessian)==0){
    cat("Could not compute confidence interval as det(hessian) = 0. Try different initial_est.")
    list(Estimate = this_fit$par)
  }else{

    # estimate confidence interval from Hessian matric
    this_se = calc_se(this_fit$hessian)
    # return estimate and CI
    list(estimate = this_fit$par,
         ci = matrix(c(this_fit$par-this_se, this_fit$par+this_se), nrow = 5))
  }
}

#' Fitting a bi-variate Hüsler-Reiss copula in a bayesian framework
#'
#' @param u (numeric) Uniform margin of copula
#' @param v (numeric) Uniform margin of copula
#' @param prior_mean (numeric) Dependence parameter of the Hüsler-Reiss copula
#' @param prior_sd (numeric) If true returns log-Likelihood, default is False
#' @param chains (numeric) see ?rstan::sampling
#' @param iter (numeric) see ?rstan::sampling
#' @param cores (numeric) see ?rstan::sampling
#' @param warmup (numeric) see ?rstan::sampling
#' @param thin (numeric) see ?rstan::sampling
#'
#' @return (list(estimate (numeric), post (numeric), fitted_model (stan model))) List of 3 elements,
#' first element is the estimate of the Hüsler-Reiss dependence parameter, calculated as the mean of the posterior distribution.
#' Second element is the posterior distribution samples. Third element is the full fitted stan model.
#' @export
#'
#' @examples
#' dat = rhr(1000, 1.2)
#' fit_hrc_bay(dat[,1], dat[,2], chains = 2)
fit_hrc_bay = function(u, v,
                       prior_mean = 0,
                       prior_sd = 1,
                       chains = 4,
                       iter = 1000,
                       warmup = floor(iter/2),
                       cores=2,
                       thin = 1){

  # load compiled stan model
  hr_bivar_stan = readRDS("stan/compiled_bivar_hrc")

  # prepare data for stan model
  stan_data = list(N=length(u),
                   u = u,
                   v = v,
                   prior_mean = prior_mean,
                   prior_sd = prior_sd)

  # sample from stan model
  fitted_model <- rstan::sampling(hr_bivar_stan,
                                  data=stan_data,
                                  chains = chains,
                                  iter = iter,
                                  cores = cores,
                                  warmup = warmup,
                                  thin = thin)

  # extract posterior distribution
  posterior_samples = rstan::extract(fitted_model)

  list(estimate = mean(posterior_samples$lambda),
       post = posterior_samples$lambda,
       fitted_model = fitted_model)
}

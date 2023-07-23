#' Cumulative distribution function of the bi-variate Hüsler-Reiss copula
#'
#' @param u (numeric) Uniform margin of copula
#' @param v (numeric) Uniform margin of copula
#' @param theta (numeric) dependence parameter of Hüsler-Reiss distribution
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
#'@export
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
#' margin_1 = qgp(copula_sample[,1], scale = 2, shape = -0.1)
#' margin_2 = qgp(copula_sample[,2], scale = 1.3, shape = -0.1)
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
         ci = matrix(c(this_fit$par-this_se, this_fit$par+this_se), nrow = 7))
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
  # hr_bivar_stan = readRDS("R/stan/compiled_bivar_hrc")
  # hr_bivar_stan = stanmodels$bivar_hrc
  # prepare data for stan model
  stan_data = list(N=length(u),
                   u = u,
                   v = v,
                   prior_mean = prior_mean,
                   prior_sd = prior_sd)

  # sample from stan model
  fitted_model <- rstan::sampling(stanmodels$bivarhrc,
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


#' Fitting a GPD in a bayesian framework
#'
#' @param x (numeric) data
#' @param v (numeric) Uniform margin of copula
#' @param scale_prior_param (numeric) rate of exp prior for scale parameter
#' @param shape_prior_lower (numeric) lower bound for uniform prior for shape parameter
#' @param shape_prior_upper (numeric) upper bound for uniform prior for shape parameter
#' @param chains (numeric) see ?rstan::sampling
#' @param iter (numeric) see ?rstan::sampling
#' @param cores (numeric) see ?rstan::sampling
#' @param warmup (numeric) see ?rstan::sampling
#' @param thin (numeric) see ?rstan::sampling
#'
#' @return (list(estimate (numeric), post (numeric), fitted_model (stan model))) List of 3 elements,
#' first element is the estimate of the GPD parameters, calculated as the mean of the posterior distribution.
#' Second element is the posterior distribution samples. Third element is the full fitted stan model.
#' @export
#'
#' @examples
#' x = rgp(1000, scale = 2, shape = -0.1)
#' fit_gpd_bay(x)
fit_gpd_bay = function(x,
                       scale_prior_param = 1,
                       shape_prior_lower = -0.5,
                       shape_prior_upper = 0.5,
                       chains = 2,
                       iter = 2000,
                       warmup = floor(iter/2),
                       cores=2,
                       thin = 1){

  stan_data = list(y = x,
                   N = length(x),
                   scale_prior_param = scale_prior_param,
                   shape_prior_upper = shape_prior_upper,
                   shape_prior_lower = shape_prior_lower)

  # sample from stan model
  fitted_model <- rstan::sampling(stanmodels$gpd,
                                  data=stan_data,
                                  chains = chains,
                                  iter = iter,
                                  cores = cores,
                                  warmup = warmup,
                                  thin = thin)

  # extract posterior distribution
  posterior_samples = rstan::extract(fitted_model)

  list(estimate = c(mean(posterior_samples$scale), mean(posterior_samples$shape)),
       post = list(scale = posterior_samples$scale,
                   shape = posterior_samples$shape),
       fitted_model = fitted_model)
}






#' Fitting a GEV in a bayesian framework
#'
#' @param x (numeric) data
#' @param prior_mean_loc (numeric) Mean of normal prior of location parameter of GEV distribution
#' @param prior_sd_loc (numeric) Standard deviation of normal prior of location parameter of GEV distribution
#' @param prior_mean_sig (numeric) Mean of normal prior of scale parameter of GEV distribution
#' @param prior_sd_sig (numeric) Standard deviation of normal prior of scale parameter of GEV distribution
#' @param prior_mean_xi (numeric) Mean of normal prior of shape parameter of GEV distribution
#' @param prior_sd_xi (numeric) Standard deviation of normal prior of shape parameter of GEV distribution
#' @param chains (numeric) see ?rstan::sampling
#' @param iter (numeric) see ?rstan::sampling
#' @param cores (numeric) see ?rstan::sampling
#' @param warmup (numeric) see ?rstan::sampling
#' @param thin (numeric) see ?rstan::sampling
#'
#' @return (list(estimate (numeric), post (numeric), fitted_model (stan model))) List of 3 elements,
#' first element is the estimate of the GEV parameters, calculated as the mean of the posterior distribution.
#' Second element is the posterior distribution samples. Third element is the full fitted stan model.
#' @export
#'
#' @examples
#' x = rgev(100, scale = 1, scale = 1, shape = -0.1)
#' fit_gev_bay(x)
fit_gev_bay = function(x,
                       prior_mean_loc = 1,
                       prior_sd_loc = 1,
                       prior_mean_xi = 0,
                       prior_sd_xi = 1,
                       prior_mean_sig = 1,
                       prior_sd_sig = 0,
                       chains = 2,
                       iter = 2000,
                       warmup = floor(iter/2),
                       cores=2,
                       thin = 1){

  stan_data = list(y = x,
                   N = length(x),
                   prior_mean_loc = prior_mean_loc,
                   prior_sd_loc = prior_sd_loc,
                   prior_mean_xi = prior_mean_xi,
                   prior_sd_xi = prior_sd_xi,
                   prior_mean_sig = prior_mean_sig,
                   prior_sd_sig = prior_sd_sig)

  # sample from stan model
  fitted_model <- rstan::sampling(stanmodels$gev,
                                  data=stan_data,
                                  chains = chains,
                                  iter = iter,
                                  cores = cores,
                                  warmup = warmup,
                                  thin = thin)

  # extract posterior distribution
  posterior_samples = rstan::extract(fitted_model)

  list(estimate = c(mean(posterior_samples$loc), mean(posterior_samples$scale), mean(posterior_samples$shape)),
       post = list(scale = posterior_samples$loc,
                   scale = posterior_samples$scale,
                   shape = posterior_samples$shape),
       fitted_model = fitted_model)
}






# negative log likelihood of copula and margins
hrc_gev_nll = function(pars, x, y){

  if(pars[2] <= 0) return(100^100)
  if((pars[3] < 0) & any(x > (pars[1] - (pars[2]/pars[3])))) return(100^100)
  if((pars[3] > 0) & any(x < (pars[1] - (pars[2]/pars[3])))) return(100^100)


  if(pars[5] <= 0) return(100^100)
  if((pars[6] < 0) & any(y > (pars[4] - (pars[5]/pars[6])))) return(100^100)
  if((pars[6] > 0) & any(y < (pars[4] - (pars[5]/pars[6])))) return(100^100)


  if(pars[7]<=0 | pars[7]>15) return(100^100)

  # likelihood of margin 1
  llm1 = dgev(x, loc = pars[1], scale = pars[2], shape = pars[3])

  # likelihood of margin 2
  llm2 = dgev(y, loc = pars[4], scale = pars[5], shape = pars[6])

  # likelihood of copula
  unif_1 = pgev(x, loc = pars[1], scale = pars[2], shape = pars[3])
  unif_2 = pgev(y, loc = pars[4], scale = pars[5], shape = pars[6])
  llc = dhr(unif_1, unif_2, lambda = pars[7])

  # full = likelihood
  LL =  -sum(log(llm1*llm2*llc))

  if(is.infinite(LL)) return(100^100)
  if(is.na(LL)) return(100^100)
  LL
}

#' Jointly fitting a bi-variate Hüsler-Reiss copula and GEV
#'
#' @param x (numeric) margin 1, assumed to be GEV distributed
#' @param y (numeric) margin 2, assumed to be GEV distributed
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
#' # Transform margins to be GEV
#' margin_1 = qgev(copula_sample[,1], loc = 1, scale = 2, shape = -0.1)
#' margin_2 = qgev(copula_sample[,2], loc = 1.5, scale = 1.3, shape = -0.1)
#'
#' # fit margins and copula jointly
#' fit_hrc_gp(margin_1, margin_2, initial_est = c(1,1,0,1,1,0,1))
fit_hrc_gev = function(x, y, initial_est){
  # --- Error handling
  # check variables are same length
  if(length(x) != length(y)){
    stop("Variables must be the same length")
  }

  this_fit = optim(fn=hrc_gev_nll,
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
         ci = matrix(c(this_fit$par-this_se, this_fit$par+this_se), nrow = 7))
  }
}


# copula_sample = hrd::rhr(5000, 2)
# #
# # Transform margins to be GPD
# margin_1 = hrd::qgp(copula_sample[,1], scale = 1, shape = -0.2)
# margin_2 = hrd::qgp(copula_sample[,2], scale = 1, shape = -0.2)
# #
# # # fit margins and copula jointly
# stan_data = list(x = margin_1,
#                  y = margin_2,
#                  prior_mean = 1,
#                  prior_sd = 1,
#                  N = length(margin_1),
#                  prior_mean_xi_1 = 0,
#                  prior_sd_xi_1 = 1,
#                  prior_mean_sig_1 = 1,
#                  prior_sd_sig_1 = 1,
#                  prior_mean_xi_2 = 0,
#                  prior_sd_xi_2 = 1,
#                  prior_mean_sig_2 = 1,
#                  prior_sd_sig_2 = 1)
# #
# fitted_model = rstan::stan(
#   file = "inst/stan/bivarhrcgpd.stan",
#   data = stan_data,
#   iter = 1000,
#   cores = 2,
#   chains = 2)
# #
# # # extract posterior distribution
# # posterior_samples = rstan::extract(fitted_model)
# #
# # list(estimate = mean(posterior_samples$lambda),
# #      post = posterior_samples$lambda,
# #      fitted_model = fitted_model)



#' Fitting a bi-variate Hüsler-Reiss copula in a bayesian framework jointly with GPD margins estimates
#'
#' @param x (numeric) margin 1, assumed to be generalised Pareto distributed
#' @param y (numeric) margin 2, assumed to be generalised Pareto distributed
#' @param prior_mean (numeric) Mean of normal prior of dependence parameter of the Hüsler-Reiss copula
#' @param prior_sd (numeric) Standard deviation of normal prior of dependence parameter of the Hüsler-Reiss copula
#' @param prior_mean_sig_1 (numeric) Mean of normal prior on scale parameter of GPD of margin 1
#' @param prior_sd_sig_1 (numeric)  Standard deviation of normal prior on scale parameter of GPD of margin 1
#' @param prior_mean_sig_2 (numeric)  Mean of normal prior on scale parameter of GPD of margin 2
#' @param prior_sd_sig_2 (numeric)  Standard deviation of normal prior on scale parameter of GPD of margin 2
#' @param prior_mean_xi_1 (numeric)  Mean of normal prior on shape parameter of GPD of margin 1
#' @param prior_sd_xi_1 (numeric)  Standard deviation of normal prior on shape parameter of GPD of margin 1
#' @param prior_mean_xi_2 (numeric)  Mean of normal prior on shape parameter of GPD of margin 2
#' @param prior_sd_xi_2 (numeric)  Standard deviation of normal prior on shape parameter of GPD of margin 2
#'
#' @param chains (numeric) see ?rstan::sampling
#' @param iter (numeric) see ?rstan::sampling
#' @param cores (numeric) see ?rstan::sampling
#' @param warmup (numeric) see ?rstan::sampling
#' @param thin (numeric) see ?rstan::sampling
#'
#' @return (list(estimate (numeric), post (list), fitted_model (stan model))) List of 3 elements,
#' first element is the estimate of the GPD margin parameters as well as the Hüsler-Reiss dependence parameter, calculated as the mean of the posterior distribution.
#' Second element are the posterior distribution samples. Third element is the full fitted stan model.
#' @export
#'
#' @examples
#' # sample from a Hüsler-Reiss copula
#' copula_sample = rhr(1000, 2)
#'
#' # Transform margins to be GPD
#' margin_1 = qgp(copula_sample[,1], scale = 1.2, shape = -0.11)
#' margin_2 = qgp(copula_sample[,2], scale = 1.3, shape = -0.09)
#'
#' # fit margins and copula jointly
#' fit_hrc_gp_bay(margin_1, margin_2)
fit_hrc_gp_bay = function(x, y,
                          prior_mean = 1,
                          prior_sd = 1,
                          prior_mean_xi_1 = 0,
                          prior_sd_xi_1 = 1,
                          prior_mean_sig_1 = 1,
                          prior_sd_sig_1 = 1,
                          prior_mean_xi_2 = 0,
                          prior_sd_xi_2 = 1,
                          prior_mean_sig_2 = 1,
                          prior_sd_sig_2 = 1,
                          chains = 2,
                          iter = 1000,
                          warmup = floor(iter/2),
                          cores=2,
                          thin = 1){

  # load compiled stan model
  # prepare data for stan model
  stan_data = list(x = x,
                   y = y,
                   prior_mean = prior_mean,
                   prior_sd = prior_sd,
                   N = length(x),
                   prior_mean_xi_1 = prior_mean_xi_1,
                   prior_sd_xi_1 = prior_sd_xi_1,
                   prior_mean_sig_1 = prior_mean_sig_1,
                   prior_sd_sig_1 = prior_sd_sig_1,
                   prior_mean_xi_2 = prior_mean_xi_2,
                   prior_sd_xi_2 = prior_sd_xi_2,
                   prior_mean_sig_2 = prior_mean_sig_2,
                   prior_sd_sig_2 = prior_sd_sig_2)

  # sample from stan model
  fitted_model <- rstan::sampling(stanmodels$bivarhrcgpd,
                                  data=stan_data,
                                  chains = chains,
                                  iter = iter,
                                  cores = cores,
                                  warmup = warmup,
                                  thin = thin)

  # extract posterior distribution
  posterior_samples = rstan::extract(fitted_model)

  list(estimate = c(mean(posterior_samples$scale_1),
                    mean(posterior_samples$shape_1),
                    mean(posterior_samples$scale_2),
                    mean(posterior_samples$shape_2),
                    mean(posterior_samples$lambda)),
       post = list(scale_1 = posterior_samples$scale_1,
                   shape_1 = posterior_samples$shape_1,
                   scale_2 = posterior_samples$scale_2,
                   shape_2 = posterior_samples$shape_2,
                   lambda = posterior_samples$lambda),
       fitted_model = fitted_model)
}




#' Fitting a bi-variate Hüsler-Reiss copula in a bayesian framework jointly with GEV margins estimates
#'
#' @param x (numeric) margin 1, assumed to be generalised extreme value distributed
#' @param y (numeric) margin 2, assumed to be generalised extreme value distributed
#' @param prior_mean (numeric) Mean of normal prior of dependence parameter of the Hüsler-Reiss copula
#' @param prior_sd (numeric) Standard deviation of normal prior of dependence parameter of the Hüsler-Reiss copula
#' @param prior_mean_loc_1 (numeric) Mean of normal prior on location parameter of GPD of margin 1
#' @param prior_sd_loc_1 (numeric)  Standard deviation of normal prior on location parameter of GPD of margin 1
#' @param prior_mean_loc_2 (numeric)  Mean of normal prior on location parameter of GPD of margin 2
#' @param prior_sd_loc_2 (numeric)  Standard deviation of normal prior on location parameter of GPD of margin 2
#' @param prior_mean_sig_1 (numeric) Mean of normal prior on scale parameter of GPD of margin 1
#' @param prior_sd_sig_1 (numeric)  Standard deviation of normal prior on scale parameter of GPD of margin 1
#' @param prior_mean_sig_2 (numeric)  Mean of normal prior on scale parameter of GPD of margin 2
#' @param prior_sd_sig_2 (numeric)  Standard deviation of normal prior on scale parameter of GPD of margin 2
#' @param prior_mean_xi_1 (numeric)  Mean of normal prior on shape parameter of GPD of margin 1
#' @param prior_sd_xi_1 (numeric)  Standard deviation of normal prior on shape parameter of GPD of margin 1
#' @param prior_mean_xi_2 (numeric)  Mean of normal prior on shape parameter of GPD of margin 2
#' @param prior_sd_xi_2 (numeric)  Standard deviation of normal prior on shape parameter of GPD of margin 2
#'
#' @param chains (numeric) see ?rstan::sampling
#' @param iter (numeric) see ?rstan::sampling
#' @param cores (numeric) see ?rstan::sampling
#' @param warmup (numeric) see ?rstan::sampling
#' @param thin (numeric) see ?rstan::sampling
#'
#' @return (list(estimate (numeric), post (list), fitted_model (stan model))) List of 3 elements,
#' first element is the estimate of the GPD margin parameters as well as the Hüsler-Reiss dependence parameter, calculated as the mean of the posterior distribution.
#' Second element are the posterior distribution samples. Third element is the full fitted stan model.
#' @export
#'
#' @examples
#' # sample from a Hüsler-Reiss copula
#' copula_sample = rhr(1000, 2)
#'
#' # Transform margins to be GPD
#' margin_1 = qgev(copula_sample[,1], loc = 1, scale = 1.2, shape = -0.11)
#' margin_2 = qgev(copula_sample[,2], loc = 1.1, scale = 1.3, shape = -0.09)
#'
#' # fit margins and copula jointly
#' fit_hrc_gev_bay(margin_1, margin_2)
fit_hrc_gev_bay = function(x, y,
                          prior_mean = 1,
                          prior_sd = 1,
                          prior_mean_xi_1 = 0,
                          prior_sd_xi_1 = 1,
                          prior_mean_sig_1 = 1,
                          prior_sd_sig_1 = 1,
                          prior_mean_xi_2 = 0,
                          prior_sd_xi_2 = 1,
                          prior_mean_sig_2 = 1,
                          prior_sd_sig_2 = 1,
                          prior_mean_loc_1 = 1,
                          prior_sd_loc_1 = 1,
                          prior_mean_loc_2 = 1,
                          prior_sd_loc_2 = 1,
                          chains = 2,
                          iter = 1000,
                          warmup = floor(iter/2),
                          cores=2,
                          thin = 1){

  # load compiled stan model
  # prepare data for stan model
  stan_data = list(x = x,
                   y = y,
                   prior_mean = prior_mean,
                   prior_sd = prior_sd,
                   N = length(x),
                   prior_mean_xi_1 = prior_mean_xi_1,
                   prior_sd_xi_1 = prior_sd_xi_1,
                   prior_mean_sig_1 = prior_mean_sig_1,
                   prior_sd_sig_1 = prior_sd_sig_1,
                   prior_mean_xi_2 = prior_mean_xi_2,
                   prior_sd_xi_2 = prior_sd_xi_2,
                   prior_mean_sig_2 = prior_mean_sig_2,
                   prior_sd_sig_2 = prior_sd_sig_2,
                   prior_mean_loc_1 = prior_mean_loc_1,
                   prior_sd_loc_1 = prior_sd_loc_1,
                   prior_mean_loc_2 = prior_mean_loc_2,
                   prior_sd_loc_2 = prior_sd_loc_2)

  # sample from stan model
  fitted_model <- rstan::sampling(stanmodels$bivarhrcgev,
                                  data=stan_data,
                                  chains = chains,
                                  iter = iter,
                                  cores = cores,
                                  warmup = warmup,
                                  thin = thin)

  # extract posterior distribution
  posterior_samples = rstan::extract(fitted_model)

  list(estimate = c(mean(posterior_samples$loc_1),
                    mean(posterior_samples$scale_1),
                    mean(posterior_samples$shape_1),
                    mean(posterior_samples$loc_2),
                    mean(posterior_samples$scale_2),
                    mean(posterior_samples$shape_2),
                    mean(posterior_samples$lambda)),
       post = list(loc_1 = posterior_samples$loc_1,
                   scale_1 = posterior_samples$scale_1,
                   shape_1 = posterior_samples$shape_1,
                   loc_2 = posterior_samples$loc_2,
                   scale_2 = posterior_samples$scale_2,
                   shape_2 = posterior_samples$shape_2,
                   lambda = posterior_samples$lambda),
       fitted_model = fitted_model)
}

#
#
# copula_sample = hrd::rhr(1000, 2)
# #
# # Transform margins to be GPD
# margin_1 = hrd::qgev(copula_sample[,1], loc = 1.2,  scale = 1, shape = -0.2)
# margin_2 = hrd::qgev(copula_sample[,2], loc = 0.9, scale = 1, shape = -0.2)
# # #
#
# hrd::fit_gevd(margin_1)
# # # # fit margins and copula jointly
# stan_data = list(x = margin_1,
#                  y = margin_2,
#                  prior_mean = 1,
#                  prior_sd = 1,
#                  N = length(margin_1),
#                  prior_mean_xi_1 = 0,
#                  prior_sd_xi_1 = 0.5,
#                  prior_mean_sig_1 = 1,
#                  prior_sd_sig_1 = 1,
#                  prior_mean_xi_2 = 0,
#                  prior_sd_xi_2 = 0.5,
#                  prior_mean_sig_2 = 1,
#                  prior_sd_sig_2 = 1,
#                  prior_mean_loc_1 = 1,
#                  prior_sd_loc_1 = 0.2,
#                  prior_mean_loc_2 = 1,
#                  prior_sd_loc_2 = 0.2)
#
# fitted_model = rstan::stan(
#   file = "inst/stan/bivarhrcgev.stan",
#   data = stan_data,
#   # init = 0,
#   iter = 1000,
#   cores = 1,
#   chains = 1)
#
#
#
#
#
# evd::pgev(margin_1, 1, 1.87918, 0.166335)
#
# evd::pgev(margin_2, -89.6998, 5.71705, -0.0620918)
#
#
#
#
# loc_1 =4.99994
# loc_2 =-89.6998
# scale_1 =0.319944
# scale_2 =5.71705
# shape_1 =-0.239771
# shape_2 =-0.0620918
# lambda =3.43341
#
#
#
# # # # extract posterior distribution
# # # posterior_samples = rstan::extract(fitted_model)
# # #
# # # list(estimate = mean(posterior_samples$lambda),
# # #      post = posterior_samples$lambda,
# # #      fitted_model = fitted_model)
#
#
#
#
#
#
#
# # Transform margins to be GPD
# x = hrd::rgev(5000, loc = 1.2,  scale = 1, shape = -0.2)
# # #
# # # # fit margins and copula jointly
# stan_data = list(x = x,
#                  N = length(x),
#                  prior_mean_xi = 0,
#                  prior_sd_xi = 1,
#                  prior_mean_sig = 1,
#                  prior_sd_sig = 1,
#                  prior_mean_loc = 1,
#                  prior_sd_loc = 1)
#
#
#
# # #
# fitted_model = rstan::stan(
#   file = "inst/stan/gev.stan",
#   data = stan_data,
#   iter = 1000,
#   init=0,
#   cores = 2,
#   chains = 2)
#
#
#
#
#
#
#
#
# margin_1
#
# v =0
# y =1.50848
# loc_1 =-1.2955
# loc_2 =2.93369
# scale_1 =1.89384
# scale_2 =0.200606
# shape_1 =-0.302694
# shape_2 =0.0659788
# lambda =7.35677
# evd::pgev(1.50848, 2.93369, 0.200606, 0.0659788)
#
#
# margin_2
# -43.0397 - (1.92204/-0.0402527)
#
# # #
# # # # extract posterior distribution
# # # posterior_samples = rstan::extract(fitted_model)
# # #
# # # list(estimate = mean(posterior_samples$lambda),
# # #      post = posterior_samples$lambda,
# # #      fitted_model = fitted_model)
#

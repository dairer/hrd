functions{
    // likelihood function of Hüsler-Reiss bivariate copula
    real HuslerReiss(real u, real v, real lambda){
      real z = log(log(u) / log(v));
      real a = (1/lambda) + (lambda/2)*z;
      real log_p = log(1 / (u * v));

    log_p += log(exp(log(u) * normal_cdf(a, 0, 1) +
                 log(v) * normal_cdf((1/lambda) + (lambda/2) * log(log(v) / log(u)), 0, 1)) *
                 (normal_cdf(1/lambda - (lambda/2)*z, 0, 1) * normal_cdf(a, 0, 1) +
                 (lambda/2) * -1 / log(v) * exp(-0.5 * a^2) / sqrt(2 * pi())));

    return log_p;
    }
}

data {
  int <lower = 0> N;
  real u[N];
  real v[N];
  real <lower=0> prior_mean;
  real <lower=0> prior_sd;
}

parameters {
  real <lower=0> lambda;
}

model {
  lambda ~ normal(prior_mean,prior_sd); // prior on HR dependence parameter
  for(n in 1:N){
    target += HuslerReiss(u[n], v[n], lambda);
  }
}




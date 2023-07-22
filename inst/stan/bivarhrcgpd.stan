functions{
   real HuslerReiss_likelihood(real u, real v, real lambda){
      real z = log(log(u) / log(v));
      real a = (1/lambda) + (lambda/2)*z;
      real log_p = log(1 / (u * v));

    log_p += log(exp(log(u) * normal_cdf(a, 0, 1) +
                 log(v) * normal_cdf((1/lambda) + (lambda/2) * log(log(v) / log(u)), 0, 1)) *
                 (normal_cdf(1/lambda - (lambda/2)*z, 0, 1) * normal_cdf(a, 0, 1) +
                 (lambda/2) * -1 / log(v) * exp(-0.5 * a^2) / sqrt(2 * pi())));

    return log_p;
    }


  // density of gpd
  real gpd_likelihood(real z, real shape, real scale){
    if (shape<0 && shape < (-scale/z))
      reject(shape, scale)
    if (scale<=0)
      reject(scale)
    if(shape == 0)
      return (exp(-z/scale) * (1/scale));
    return (1 / (scale * ( 1 + shape * z/scale)^(1/shape + 1)));
  }

  // gpd cdf
  real gpd_distribution(real z, real shape, real scale){
    if (shape<0 && (z > (-scale/shape)))
      reject(shape, scale)
    if (scale<=0)
      reject(scale)
    if(shape == 0)
      return (1-exp(-z/scale));
    return (1-(1+((shape*z)/scale))^(-1/shape));
  }

    // likelihood function of Hüsler-Reiss bivariate copula
    real HuslerReiss(real x, real y,
                     real scale_1, real scale_2,
                     real shape_1, real shape_2,
                     real lambda){

      real u;
      real v;
      real density_marg_1;
      real density_marg_2;
      real z;
      real a;
      real density_HR;
      real LL;

      if(lambda<=0){
        reject(lambda);
      }

      if(shape_1<0 && shape_1 < (-scale_1/x)){
        reject(shape_1, scale_1);
      }

      if (scale_1<=0){
         reject(scale_1);
      }


      if (shape_2<0 && shape_2 < (-scale_2/y)){
        reject(shape_2, scale_2);
      }

      if (scale_2<=0){
        reject(scale_2);
      }


      // likelihood margin 1
      density_marg_1 = gpd_likelihood(x, shape_1, scale_1);
      if(density_marg_1 <=0){
        reject(shape_1, scale_1);
      }


      // likelihood margin 2
      density_marg_2 = gpd_likelihood(y, shape_2, scale_2);
      if(density_marg_2 <=0){
        reject(shape_2, scale_2);
      }


      // transform data to uniform
      u = gpd_distribution(x, shape_1, scale_1);
      v = gpd_distribution(y, shape_2, scale_2);

      if(u<=0){
        reject(shape_1, scale_1);
      }
      if(u>=1){
        reject(shape_1, scale_1);
      }
      if(v<=0){
        reject(shape_2, scale_2);
      }
      if(v>=1){
        reject(shape_2, scale_2);
      }



      z = log(log(u) / log(v));
      a = (1/lambda) + (lambda/2)*z;


      density_HR = (1 / (u * v))*(exp(log(u) * normal_cdf(a, 0, 1) +
      log(v) * normal_cdf((1/lambda) + (lambda/2) * log(log(v) / log(u)), 0, 1)) *
      (normal_cdf(1/lambda - (lambda/2)*z, 0, 1) * normal_cdf(a, 0, 1) +
      (lambda/2) * -1 / log(v) * exp(-0.5 * a^2) / sqrt(2 * pi())));


    if(density_HR<=0){
      reject(shape_1, scale_1, shape_2, scale_2, lambda);
    }

    LL = log(density_HR*density_marg_1*density_marg_2);

    if(LL < (-10^50)){
      reject(shape_1, scale_1, shape_2, scale_2);
    }

    return(LL);
    }
}

data {
  int <lower = 0> N;
  real<lower = 0> x[N];
  real<lower = 0> y[N];
  real <lower=0> prior_mean;
  real <lower=0> prior_sd;
  real <lower=-0.4, upper = 0.4> prior_mean_xi_1;
  real <lower=0>prior_sd_xi_1;
  real <lower=0>prior_mean_sig_1;
  real <lower=0>prior_sd_sig_1;
  real <lower=-0.4, upper = 0.4>prior_mean_xi_2;
  real <lower=0>prior_sd_xi_2;
  real <lower=0>prior_mean_sig_2;
  real <lower=0>prior_sd_sig_2;
  }

transformed data {
  real maxObs1 = max(x);
  real maxObs2 = max(y);
}


// parameters {
//   real <lower=0, upper = 5> lambda;
//   real <lower=0.01, upper = 5> scale_1;
//   real <lower=-0.49, upper = 0.49> shape_1;
//   real <lower=0.01, upper = 5> scale_2;
//   real <lower=-0.49, upper = 0.49> shape_2;
// }
//


parameters {
  real <lower=0, upper = 10> lambda;
  real <lower=0> scale_1;
  real <lower=0> scale_2;
  real <lower=-scale_1/(maxObs1)> shape_1;
  real <lower=-scale_2/(maxObs2)> shape_2;

}



model {
  lambda ~ normal(prior_mean,prior_sd); // prior on HR dependence parameter
  shape_1 ~ normal(prior_mean_xi_1,prior_sd_xi_1);
  scale_1 ~ normal(prior_mean_sig_1,prior_sd_sig_1);
  shape_2 ~ normal(prior_mean_xi_2,prior_sd_xi_2);
  scale_2 ~ normal(prior_mean_sig_2,prior_sd_sig_2);

  for(n in 1:N){
    target += HuslerReiss(x[n], y[n], scale_1,scale_2,shape_1,shape_2,lambda);
  }
}

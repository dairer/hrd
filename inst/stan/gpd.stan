functions{
  // log of the density function of a GEV
  real gpd_likelihood(real y, real shape, real scale){
    if (shape<0 && shape < (-scale/y))
      reject(shape, scale)
    if (scale<=0)
      reject(scale)
    if(shape == 0)
      return log(exp(-y/scale) * (1/scale));
    return log(1 / (scale * ( 1 + shape * y/scale)^(1/shape + 1)));
  }
}


data {
  int<lower = 0> N;
  real<lower = 0> y[N];
  real <lower = 0> scale_prior_param;
  real shape_prior_upper;
  real shape_prior_lower;

}

transformed data {
  real maxObs = max(y);
}

parameters {
  real <lower=0> scale;
  real <lower=-scale/(maxObs)> shape;
}


model {
  // priors
  shape ~ uniform(shape_prior_lower, shape_prior_upper);
  scale ~ exponential(scale_prior_param);

  // likelihood (can i vectorise this?)
  for(n in 1:N){
    target += gpd_likelihood(y[n], shape, scale);
  }
}

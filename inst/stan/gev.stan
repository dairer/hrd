functions{

  // density of gev
  real gev_likelihood(real z, real loc, real shape, real scale){

    real tm;

    if(scale <= 0)
      reject(scale)

    if((shape < 0) && (z > (loc - (scale/shape))))
      reject(loc, scale, shape)

    if((shape > 0) && (z < (loc - (scale/shape))))
          reject(loc, scale, shape)

    if(shape != 0){
      tm = (1 + shape*((z - loc)/scale))^(-1/shape);
    }else{
      tm = exp(-(z - loc)/scale);
    }

    return log((1/scale)*(tm^(shape+1))*exp(-tm));
  }

}

data {
  int <lower = 0> N;
  real x[N];
  real <lower=-0.5, upper = 0.5> prior_mean_xi;
  real <lower=0>prior_sd_xi;
  real <lower=0>prior_mean_sig;
  real <lower=0>prior_sd_sig;
  real <lower=0>prior_mean_loc;
  real <lower=0>prior_sd_loc;
  }

transformed data {
  real maxObs = max(x);
  real minObs = min(x);
}


parameters {
  real <lower=0> scale;
  real <lower=-0.5, upper = 0.5> shape;

// //equiv to if_else(shape < 0,  ((scale/shape) + maxObs1), negative_infinity())
 real<lower= ((shape < 0) ?  ((scale/shape) + maxObs) :  negative_infinity()),
      upper=((shape > 0) ?   ((scale/shape) + minObs): positive_infinity()) > loc;

 //
 // real<lower=if_else( xi < 0, min_y + sigma / xi, negative_infinity() ),
 //       upper=if_else( xi > 0, positive_infinity(), max_y + sigma / xi )> mu;


}


model {
  loc ~ normal(prior_mean_loc,prior_sd_loc);
  shape ~ normal(prior_mean_xi,prior_sd_xi);
  scale ~ normal(prior_mean_sig,prior_sd_sig);

  for(n in 1:N){
    target += gev_likelihood(x[n], loc, shape, scale);
  }
}


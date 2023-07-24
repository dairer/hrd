functions{
  // density of gev
  real gev_lpdf(vector z, real loc, real scale, real shape){

    int n = rows(z);
    vector[n] a;
    vector[n] b;
    vector[n] c;

    real LL = 0;

    if(shape > 0.00001 || shape < -0.00001){
      for(i in 1:n){
        a[i] = 1+shape*((z[i] - loc)/scale);
        b[i] = pow(a[i], (-1/shape));
        c[i] = log(a[i]);
      }
      LL = ((-1*n*log(scale)) - ((1+(1/shape))*sum(c)) - sum(b));
    }else{
      for(i in 1:n){
        a[i] = ((z[i] - loc)/scale);
        b[i] = exp(-1*a[i]);
      }
     LL = (-n*log(scale) - sum(a) - sum(b));
    }
    return(LL);
  }
}



data {
  int <lower = 0> N;
  vector[N] x;
  real <lower=0>prior_mean_sig;
  real <lower=0>prior_sd_sig;
  real prior_mean_loc;
  real <lower=0>prior_sd_loc;
  real shape_lower;
  real shape_upper;

  }

transformed data {
  real maxObs = max(x);
  real minObs = min(x);
}


parameters {
  real <lower=0> scale;
  real <lower=-0.5, upper = 0.5> shape;
  // real  loc;
  // real<lower=if_else( shape > 0, minObs, negative_infinity()),
  //      upper=if_else( shape > 0, positive_infinity(), maxObs )> loc;
  real<lower=(shape > 0 ? minObs : negative_infinity()),
       upper=(shape > 0 ? positive_infinity() : maxObs )> loc;




// // //equiv to if_else(shape < 0,  ((scale/shape) + maxObs1), negative_infinity())
//  real<lower= ((shape < 0) ?  ((scale/shape) + maxObs) :  negative_infinity()),
//       upper=((shape > 0) ?   ((scale/shape) + minObs): positive_infinity()) > loc;

 // real<lower= ((shape < -0.00001) ? (maxObs + (scale/shape)) :  negative_infinity()),
 //      upper=((shape > 0.00001) ?  (minObs + (scale/shape)): positive_infinity()) > loc;
}


model {
  loc ~ normal(prior_mean_loc,prior_sd_loc);
  shape ~ uniform(shape_lower,shape_upper)T[shape_lower,shape_upper];
  scale ~ normal(prior_mean_sig,prior_sd_sig);
  x ~ gev(loc, scale, shape);
}

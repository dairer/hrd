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

    return ((1/scale)*(tm^(shape+1))*exp(-tm));
  }

  // gev cdf
  real gev_distribution(real z, real loc, real shape, real scale){
    real tm;

    if(scale <= 0)
      reject(scale)

    if((shape < 0) && (z > (loc - (scale/shape))))
      reject(loc, scale, shape)

    if((shape > 0) && (z < (loc - (scale/shape))))
        reject(loc, scale, shape)

    if(shape != 0){
      if((1 + shape*((z - loc)/scale)) > 0){
        tm = (1 + shape*((z - loc)/scale))^(-1/shape);
      }
      else{
        tm = 0^(-1/shape);
      }
    }else{
      tm = exp(-(z - loc)/scale);
    }

    return(exp(-tm));
  }

    // likelihood function of Hüsler-Reiss bivariate copula
  real HuslerReiss(real x, real y, real loc_1, real loc_2, real scale_1,
                   real scale_2, real shape_1, real shape_2, real lambda){

      real u;
      real v;
      real density_marg_1;
      real density_marg_2;
      real z;
      real a;
      real density_HR;
      real LL;

//   print("loc_1 =", loc_1);
//   print("loc_2 =", loc_2);
//   print("scale_1 =", scale_1);
//   print("scale_2 =", scale_2);
//   print("shape_1 =", shape_1);
//   print("shape_2 =", shape_2);
//   print("lambda =", lambda);

      if(lambda<=0){
        reject(lambda);
      }

      if(scale_1 <= 0)
        reject(loc_1, scale_1, shape_1)

      if((shape_1 < 0) && (x > (loc_1 - (scale_1/shape_1))))
        reject(loc_1, scale_1, shape_1)

      if((shape_1 > 0) && (x < (loc_1 - (scale_1/shape_1))))
            reject(loc_1, scale_1, shape_1)

      if(scale_2 <= 0)
        reject(loc_2, scale_2, shape_2)

      if((shape_2 < 0) && (y > (loc_2 - (scale_2/shape_2))))
        reject(loc_2, scale_2, shape_2)

      if((shape_2 > 0) && (y < (loc_2 - (scale_2/shape_2))))
            reject(loc_2, scale_2, shape_2)

      // likelihood margin 1
      density_marg_1 = gev_likelihood(x, loc_1, shape_1, scale_1);
      if(density_marg_1 < 0 || density_marg_2 >=1){
        density_marg_1 = 0;
//
//           print("loc_1 =", loc_1);
//           print("loc_2 =", loc_2);
//           print("scale_1 =", scale_1);
//           print("scale_2 =", scale_2);
//           print("shape_1 =", shape_1);
//           print("shape_2 =", shape_2);
//           print("lambda =", lambda);


        // reject(loc_1, shape_1, scale_1);
      }

      // likelihood margin 2
      density_marg_2 = gev_likelihood(y, loc_2, shape_2, scale_2);
      if(density_marg_2 <= 0 || density_marg_2 >=1){

          density_marg_2 = 0;
          // print("loc_1 =", loc_1);
          // print("loc_2 =", loc_2);
          // print("scale_1 =", scale_1);
          // print("scale_2 =", scale_2);
          // print("shape_1 =", shape_1);
          // print("shape_2 =", shape_2);
          // print("lambda =", lambda);
//
//         reject(loc_2, shape_2, scale_2);
      }


      if(density_marg_2 > 0 && density_marg_1 > 0){

              // transform data to uniform
      u = gev_distribution(x, loc_1, shape_1, scale_1);
      v = gev_distribution(y, loc_2, shape_2, scale_2);

      // if(u<=0){
      //
      //   u = 0 + 1e-9;
      //     // print("u =", u);
      //     // print("x =", x);
      //     //
      //     // print("loc_1 =", loc_1);
      //     // print("loc_2 =", loc_2);
      //     // print("scale_1 =", scale_1);
      //     // print("scale_2 =", scale_2);
      //     // print("shape_1 =", shape_1);
      //     // print("shape_2 =", shape_2);
      //     // print("lambda =", lambda);
      //
      //   // reject(loc_1, shape_1, scale_1);
      // }
      // if(u>=1){
      //   u = 1 - 1e-9;
      //
      //   // print("u =", u);
      //   //           print("x =", x);
      //   //
      //   // print("loc_1 =", loc_1);
      //   // print("loc_2 =", loc_2);
      //   // print("scale_1 =", scale_1);
      //   // print("scale_2 =", scale_2);
      //   // print("shape_1 =", shape_1);
      //   // print("shape_2 =", shape_2);
      //   // print("lambda =", lambda);
      //
      //   // reject(loc_1, shape_1, scale_1);
      // }
      // if(v<=0){
      //   v = 0 + 1e-9;
      //
      //   // print("v =", v);
      //   //           print("y =", y);
      //   //
      //   // print("loc_1 =", loc_1);
      //   // print("loc_2 =", loc_2);
      //   // print("scale_1 =", scale_1);
      //   // print("scale_2 =", scale_2);
      //   // print("shape_1 =", shape_1);
      //   // print("shape_2 =", shape_2);
      //   // print("lambda =", lambda);
      //
      //   // reject(loc_2, shape_2, scale_2);
      // }
      // if(v>=1){
      //   v = 1- 1e-9;
      //
      //   // print("v =", v);
      //   //           print("y =", y);
      //   //
      //   // print("loc_1 =", loc_1);
      //   // print("loc_2 =", loc_2);
      //   // print("scale_1 =", scale_1);
      //   // print("scale_2 =", scale_2);
      //   // print("shape_1 =", shape_1);
      //   // print("shape_2 =", shape_2);
      //   // print("lambda =", lambda);
      //
      //   // reject(loc_2, shape_2, scale_2);
      // }


      z = log(log(u) / log(v));
      a = (1/lambda) + (lambda/2)*z;


      density_HR = (1 / (u * v))*(exp(log(u) * normal_cdf(a, 0, 1) +
      log(v) * normal_cdf((1/lambda) + (lambda/2) * log(log(v) / log(u)), 0, 1)) *
      (normal_cdf(1/lambda - (lambda/2)*z, 0, 1) * normal_cdf(a, 0, 1) +
      (lambda/2) * -1 / log(v) * exp(-0.5 * a^2) / sqrt(2 * pi())));


      LL = log(density_HR*density_marg_1*density_marg_2);
        return(LL);
      }else{
        return(negative_infinity());
      }

      // }
//
//     if(LL < (-1^10)){
//       LL = -1^10;
//       // print("LL =", LL);
//       // print("u =", u);
//       // print("v =", v);
//       //
//       // print("loc_1 =", loc_1);
//       // print("loc_2 =", loc_2);
//       // print("scale_1 =", scale_1);
//       // print("scale_2 =", scale_2);
//       // print("shape_1 =", shape_1);
//       // print("shape_2 =", shape_2);
//       // print("lambda =", lambda);
//       //
//       // reject(loc_1, loc_2, shape_1, scale_1, shape_2, scale_2);
//     }
//
//     return(LL);
    }
}

data {
  int <lower = 0> N;
  real x[N];
  real y[N];
  real <lower=0> prior_mean;
  real <lower=0> prior_sd;
  real <lower=-0.4, upper = 0.4> prior_mean_xi_1;
  real <lower=0>prior_sd_xi_1;
  real <lower=0>prior_mean_sig_1;
  real <lower=0>prior_sd_sig_1;
  real <lower=0>prior_mean_loc_1;
  real <lower=0>prior_sd_loc_1;
  real <lower=-0.4, upper = 0.4>prior_mean_xi_2;
  real <lower=0>prior_sd_xi_2;
  real <lower=0>prior_mean_sig_2;
  real <lower=0>prior_sd_sig_2;
  real <lower=0>prior_mean_loc_2;
  real <lower=0>prior_sd_loc_2;
  }




transformed data {
  real maxObs1 = max(x);
  real maxObs2 = max(y);
  real minObs1 = min(x);
  real minObs2 = min(y);

}

// parameters {
//   real <lower=0, upper = 5> lambda;
//   real <lower=0.01, upper = 5> scale_1;
//   real <lower=-0.49, upper = 0.49> shape_1;
//   real <lower=0.01, upper = 5> scale_2;
//   real <lower=-0.49, upper = 0.49> shape_2;
// }


parameters {
  real <lower=0, upper = 10> lambda;
  real <lower=0> scale_1;
  real <lower=0> scale_2;
  real <lower=-0.5, upper = 0.5> shape_1;
  real <lower=-0.5, upper = 0.5> shape_2;

// //equiv to if_else(shape < 0,  ((scale/shape) + maxObs1), negative_infinity())
 real<lower= ((shape_1 < 0) ? ((scale_1/shape_1) + maxObs1) :  negative_infinity()),
      upper= ((shape_1 > 0) ? ((scale_1/shape_1) + minObs1) : positive_infinity()) > loc_1;

 real<lower= ((shape_2 < 0) ? ((scale_2/shape_2) + maxObs2) :  negative_infinity()),
      upper= ((shape_2 > 0) ? ((scale_2/shape_2) + minObs2) : positive_infinity()) > loc_2;
}

model {
  lambda ~ normal(prior_mean,prior_sd); // prior on HR dependence parameter

  loc_1 ~ normal(prior_mean_loc_1,prior_sd_loc_1);
  shape_1 ~ normal(prior_mean_xi_1,prior_sd_xi_1);
  scale_1 ~ normal(prior_mean_sig_1,prior_sd_sig_1);

  loc_2 ~ normal(prior_mean_loc_2,prior_sd_loc_2);
  shape_2 ~ normal(prior_mean_xi_2,prior_sd_xi_2);
  scale_2 ~ normal(prior_mean_sig_2,prior_sd_sig_2);

  for(n in 1:N){
    target += HuslerReiss(x[n], y[n], loc_1, loc_2, scale_1,scale_2,shape_1,shape_2,lambda);
  }
}

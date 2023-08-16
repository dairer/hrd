
#' Log likelihood of spectral density of Brown-Resnick process as in Engelke (2015).
#'
#' @param pars (numeric) Parameter values of variogram which parameterises the Brown-Resnick process.
#' @param dt (list) Observations of extreme process, with standarised margins.
#' @param lcs (matrix) Locations of observations in dt.
#' @param vr (numeric) Variogram function, evaluates with "pars"
#' @param conditioned.site (numeric) Site that is "conditioned" on.
#'
#' @return (numeric) Log-liklihood
#' @export
br_ll = function(pars, dt, lcs, vr, conditioned.site = 1){
  nlocs = nrow(lcs) # number of sites
  loc.id.pairs = expand.grid(seq(nlocs) ,seq(nlocs))
  loc.pairs = cbind(lcs[loc.id.pairs[,1],], lcs[loc.id.pairs[,2],]) # all pairs of locations

  # get distance between all pairs
  dstncs = loc.pairs %>%
    apply(MARGIN = 1, FUN = function(x) sqrt((x[3] - x[1])^2 + (x[4] - x[2])^2))

  Lambda = (vr(dstncs, pars)) %>% matrix(nrow = nlocs, byrow = T)
  conditioned.site.lambda = Lambda[conditioned.site,]
  Lambda = Lambda[-conditioned.site, -conditioned.site]

  psi = outer(conditioned.site.lambda[-1], conditioned.site.lambda[-1], "+") - Lambda # "covarianve" matrix
  inversePsi = psi %>% solve()  # inverse of "covarianve" matrix

  detPsi = determinant(psi)$modulus[1] # calculates log det by default

  omegas = log(t(sapply(dt, function(x) x[-conditioned.site]))/sapply(dt, "[[", conditioned.site)) %>% sweep(2, conditioned.site.lambda[-1], "+")

  if(nlocs == 2) omegas = t(omegas)

  summ = apply(omegas, MARGIN = 1, function(x){t(x) %*% inversePsi %*% (x)}) %>% sum
  0.5*length(dt)*detPsi+ 0.5*summ # Log liklihood
}

#' Log likelihood of spectral density of Brown-Resnick process with changing dimentions. Relies of function br_ll()
#'
#' @param pars (numeric) Parameter values of variogram which parameterises the Brown-Resnick process.
#' @param dt (list) Observations of extreme process, with standarised margins.
#' @param lcs (list) List of matrices that describe locations of observations in corresponding position in list dt.
#' @param vr (numeric) Variogram function, evaluates with "pars"
#' @param conditioned.site (numeric) Site that is "conditioned" on.
#'
#' @return (numeric) Log-liklihood
#' @export
change_dim_br_ll = function(pars, dt, lcs, vr, conditioned.site = 1){
  LL = 0
  R = length(dt)
  for(r in seq(R)){
    this_dat = dt[r]
    this_locs = lcs[r][[1]]
    LL =  LL+br_ll(pars, this_dat, this_locs, vr, conditioned.site = conditioned.site)
  }
  return(LL)
}


ngll = function(pars, lcs, dt, vr, conditioned.site){
  if(class(lcs) != "list"){ # no missing observations
    br_ll(pars, lcs = lcs, dt = dt, vr = vr, conditioned.site = conditioned.site)
  }else{
    change_dim_br_ll(pars, dt, lcs, vr = vr, conditioned.site = conditioned.site)
  }
}

#' Estimate variogram parameters for Brown-Resnick process with changing or constant dimensions.
#'
#' @param dt (list) Observations of extreme process, with standarised margins.
#' @param lcs (list) List of matrices that describe locations of observations in corresponding position in list dt.
#' @param vr (numeric) Variogram function, evaluates with "pars"
#' @param conditioned.site (numeric) Site that is "conditioned" on.
#'
#' @return (numeric) Log-liklihood
#' @export
fit_rpareto_br = function(dt, lcs, vr, conditioned.site = 1){
  optim(par = c(1 ,1), fn = ngll, hessian = T, lcs = lcs, dt = dt, vr = vr, conditioned.site = conditioned.site)
}

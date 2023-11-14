# internal
# this function fits the husler reiss to pair wise combinations transformed to have uniform scale.
#' @noRd
hrd_exploritory = function(dat, data_scale = T){
  lambda = NULL
  if(data_scale){
    get_emp_rank = function(x) rank(x)/(length(x)+1)

    dat = dat %>%
      base::apply(MARGIN = 2, FUN = get_emp_rank) %>%
      dplyr::as_tibble()
  }

  combs = base::names(dat) %>% utils::combn(2)
  combs = ncol(dat) %>% seq %>% utils::combn(2)

  all_pairs = data.frame()

  for(i in seq(ncol(combs))){
    dat1 = dat[,combs[,i][1]]
    dat2 = dat[,combs[,i][2]]

    all_pairs = rbind(all_pairs, c(base::names(dat1), base::names(dat2), hrd::fit_hrc(unlist(dat1), unlist(dat2), initial_est = 0)$estimate))
  }

  base::names(all_pairs) = c("s1", "s2", "lambda")

  all_pairs %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(lambda = as.numeric(lambda)) %>%
    dplyr::mutate(chi = 2*stats::dnorm(1/lambda))
}

#' Estimate and plot matrix of dependence parameter of Hüsler-Reiss
#'
#' This function plots the estimated dependence matrix of the Hüsler-Reiss copula where data is on any scale. The function transforms the data to uniform margins using an empirical ranking.
#'
#' @param dat (data.frame or matrix) Each column corresponds to a different site/variable, with component-wise maxima in each row.
#' @param lower_diag If TRUE, the default, plot lower diagonal elements of dependence matrix
#'
#' @return (ggplot2::ggplot) Plotted dependence matrix
#' @export
explore_lambda = function(dat, lower_diag = T){

  s1=s2=lambda=NULL

  all_pairs = hrd_exploritory(dat)

  if(!lower_diag){
    all_pairs = rbind(
      all_pairs,
      all_pairs %>%
        dplyr::rename(s1 = s2,
               s2 = s1) %>%
        dplyr::select(s1, s2, lambda))
  }

  all_pairs$s1 = base::factor(all_pairs$s1, levels = base::unique(all_pairs$s1) %>% sort)
  all_pairs$s2 = base::factor(all_pairs$s2, levels = base::unique(all_pairs$s2) %>% sort %>% rev)


  all_pairs%>%
    ggplot2::ggplot()+
    ggplot2::geom_tile(ggplot2::aes(s1, s2, fill = lambda))+
    ggplot2::labs(fill = expression(lambda))+
    ggplot2::theme_minimal(12)+
    ggplot2::scale_fill_gradientn(colours = RColorBrewer::brewer.pal(10,"Spectral"))+
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust=0.5, hjust=1))

}

#' Estimate and plot matrix pairwise dependence measure
#'
#' This function plots the estimated extremal coefficient (chi) matrix derived from an estimated dependence matrix of the Hüsler-Reiss copula where data is on any scale. The function transforms the data to uniform margins using an empirical ranking.
#'
#' @param dat (data.frame or matrix) Each column corresponds to a different site/variable, with component-wise maxima in each row.
#' @param lower_diag If TRUE, the default, plot lower diagonal elements of extremal coefficient (chi) matrix
#' @param lims (numeric) Vector of length 2 controlling the limits of the scale for plotting extremal coefficient. Default plot has limit c(0,1).
#'
#' @return (ggplot2::ggplot) Plotted extremal coefficient (chi) matrix
#' @export
explore_chi = function(dat, lower_diag = T, lims = NA){
  s1=s2=lambda=chi=NULL

  all_pairs = hrd_exploritory(dat)

  if(!lower_diag){
    all_pairs = rbind(
      all_pairs,
      all_pairs %>%
        dplyr::rename(s1 = s2,
               s2 = s1) %>%
        dplyr::select(s1, s2, lambda, chi))
  }

  all_pairs$s1 = base::factor(all_pairs$s1, levels = base::unique(all_pairs$s1) %>% sort)
  all_pairs$s2 = base::factor(all_pairs$s2, levels = base::unique(all_pairs$s2) %>% sort %>% rev)


  plt = all_pairs%>%
    ggplot2::ggplot()+
    ggplot2::geom_tile(ggplot2::aes(s1, s2, fill = chi))+
    ggplot2::labs(fill = expression(chi))+
    ggplot2::theme_minimal(12)+
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust=0.5, hjust=1))

  if(length(lims) == 1){
    plt + ggplot2::scale_fill_gradientn(colours = RColorBrewer::brewer.pal(10,"Spectral"))
  }else{
    plt + ggplot2::scale_fill_gradientn(colours = RColorBrewer::brewer.pal(10,"Spectral"), limits = lims)
  }
}

#' Calculate pairwise dependence measure
#'
#' Internal function used to estimate emprical pairwise extremal dependance coefficient
#'
#' @param dat (data.frame or matrix) Each column corresponds to a different site/variable, with observations in each row.
#' @param u (numeric) Quantile at which to estimate extremal dependence coefficient
#' @param lower_diag If TRUE, the default, plot lower diagonal elements of extremal coefficient (chi) matrix
#' @param data_scale If TRUE, the default, function transforms the data to uniform margins using an empirical ranking.
#' @NoRd
emp_chi = function(dat, u, data_scale = T, lower_diag = T){
  s1=s2=lambda=NULL
  if(data_scale){
    get_emp_rank = function(x) rank(x)/(length(x)+1)

    dat = dat %>%
      base::apply(MARGIN = 2, FUN = get_emp_rank) %>%
      dplyr::as_tibble()
  }

  combs = base::names(dat) %>% utils::combn(2)
  combs = ncol(dat) %>% seq %>% utils::combn(2)

  all_pairs = data.frame()

  for(i in seq(ncol(combs))){

    dat1 = dat[,combs[,i][1]]
    dat2 = dat[,combs[,i][2]]
    all_pairs = rbind(all_pairs, c(base::names(dat1), base::names(dat2), sum(unlist(dat1) > u & unlist(dat2) > u)/sum(unlist(dat1) > u)))
  }

  base::names(all_pairs) = c("s1", "s2", "emp_chi")

  if(!lower_diag){
    all_pairs = rbind(
      all_pairs,
      all_pairs %>%
        dplyr::rename(s1 = s2,
               s2 = s1) %>%
        dplyr::select(s1, s2, lambda))
  }

  all_pairs$s1 = base::factor(all_pairs$s1, levels = base::unique(all_pairs$s1) %>% sort)
  all_pairs$s2 = base::factor(all_pairs$s2, levels = base::unique(all_pairs$s2) %>% sort %>% rev)



  all_pairs %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(emp_chi = as.numeric(emp_chi))
}

#' Estimate (empirically) and plot matrix pairwise dependence measure
#'
#' This function plots the extremal coefficient (chi) matrix estimated empirically.
#'
#' @param dat (data.frame or matrix) Each column corresponds to a different site/variable, each row is an observation.
#' @param u (numeric) Quantile at which to estimate extremal dependence coefficient
#' @param data_scale If TRUE, the default, function transforms the data to uniform margins using an empirical ranking.
#' @param lower_diag If TRUE, the default, plot lower diagonal elements of extremal coefficient (chi) matrix
#' @param lims (numeric) Vector of length 2 controlling the limits of the scale for plotting extremal coefficient. Default plot has limit c(0,1).
#'
#' @return (ggplot2::ggplot) Plotted extremal coefficient (chi) matrix
#' @export
explore_emp_chi = function(dat, u = 0.8, data_scale = T, lower_diag = T, lims = NA){
  s1=s2=NULL
  plt = emp_chi(dat, u, data_scale = T) %>%
    ggplot2::ggplot()+
    ggplot2::geom_tile(ggplot2::aes(s1, s2, fill = emp_chi))+
    ggplot2::labs(fill = expression(emp_chi))+
    ggplot2::theme_minimal(12)+
    ggplot2::labs(fill = expression(chi[u]))+
    ggplot2::theme(axis.title = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 90,
                                     vjust=0.5, hjust=1))
  if(length(lims) == 1){
    plt + ggplot2::scale_fill_gradientn(colours = RColorBrewer::brewer.pal(10,"Spectral"))
  }else{
    plt + ggplot2::scale_fill_gradientn(colours = RColorBrewer::brewer.pal(10,"Spectral"), limits = lims)
  }
}





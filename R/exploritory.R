#
# get_emp_rank = function(x) rank(x)/(length(x)+1)
#
# emp_chi_univar = function(x,y,u) sum(x>u & y>u)/sum(x>u)


#' Calculate exmpirical extremal dependence coefficient
#'
#' @param dat (data.frame or matrix) Each column corresponds to a different site/ variable, if ncol(dat) == 2, then the function returns an estimate of the extremal coefficient over a range of u. If ncol(dat)>2, the function returns an estimate of the extremal coefficient for one value of u, for all pairs of sites/variables.
#' @param u (vector or numeric) Probability to estimate the extremal coefficient for. If ncol(dat) == 2 then u should be a vector. If ncol(dat) > 2 then u should be a single numeric value.
#'
#' @return (data.frame)
#' @export
# emp_chi = function(dat, u = 0.9){
#   # can give matrix of 2 or more
#
#   # data frame to store chi
#   all_pairs = data.frame()
#
#   # transform data to uniform margins empirically
#   dat = dat %>% apply(MARGIN = 2, FUN = get_emp_rank) %>% as_tibble()
#
#   if(ncol(dat) > 2){
#     # here we have more than one sites,
#     # look at all pairs chi at single u
#
#     if(length(u) > 1) stop("u must be unviraite")
#
#     # all pairs of sites
#     combs = ncol(dat) %>% seq %>% combn(2)
#
#     for(i in seq(ncol(combs))){
#       dat1 = dat[,combs[,i][1]]
#       dat2 = dat[,combs[,i][2]]
#       all_pairs = rbind(all_pairs, c(names(dat1),
#                                      names(dat2),
#                                      emp_chi_univar(unlist(dat1), unlist(dat2), u)))
#     }
#
#     names(all_pairs) = c('s1', 's2', 'chi')
#     all_pairs %>% mutate(chi = as.numeric(chi))
#
#   }else{
#     # just two sites, look at chi over range of u
#     if(length(u) <= 1) u = seq(0.1, 0.95, by = 0.05)
#
#     dat1 = dat[,1] %>% unlist
#     dat2 = dat[,2] %>% unlist
#
#     for(i in u){
#       all_pairs = rbind(all_pairs, c(i, emp_chi_univar(dat1, dat2, i)))
#     }
#
#     names(all_pairs) = c('u', 'chi')
#     all_pairs %>% mutate(chi = as.numeric(chi))
#   }
# }
#



#
# plt_lambda_matrix = function(dat){
#
#   dat$s1 = factor(dat$s1, levels = unique(dat$s1) %>% sort)
#   dat$s2 = factor(dat$s2, levels = unique(dat$s2) %>% sort %>% rev)
#
#   dat%>%
#     ggplot()+
#     geom_tile(aes(s1, s2, fill = lambda))+
#     labs(fill = expression(lambda))+
#     theme_minimal(12)+
#     ggplot2::scale_fill_gradientn(colours = RColorBrewer::brewer.pal(10,"Spectral"))+
#     theme(axis.title = element_blank(),
#           axis.text.x = element_text(angle = 90,
#                                      vjust=0.5, hjust=1))
# }





# chi_matrix = function(dat, lower_diag = T){
#
#   # check if lambda is already calculated
#   if(!(ncol(dat) == 3 & names(dat) == c('s1', 's2', 'chi'))) dat = emp_chi(dat)
#
#   if(!lower_diag){
#     dat = rbind(
#       dat,
#       dat %>%
#         rename(s1 = s2,
#                s2 = s1) %>%
#         dplyr::select(s1, s2, lambda))
#   }
#
#   dat %>% plt_lambda_matrix()
# }
# chi_matrix(data)

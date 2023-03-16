rstan_options(auto_write = TRUE) # avoid recompilation of unchanged Stan programs

# ---- compile and save model
hr_bivar_stan = rstan::stan_model("stan/bivar_hrc.stan")
hr_bivar_stan %>% saveRDS("stan/compiled_bivar_hrc")

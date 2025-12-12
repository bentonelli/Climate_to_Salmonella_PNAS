
library(rstan)
library(MCMCvis)
library(shinystan)

data_mdl <- readRDS("data/3_model_in/model_full_data_2024.rds")

fitmast <- stan(
  file = "code/2_model_run/2a_cmid.stan",  # Stan program
  data = data_mdl,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2000,          # number of warmup iterations per chain
  iter = 8000,            # total number of iterations per chain
  cores = 4, # number of cores (could use one per chain)
  thin = 2, 
  refresh = 10,
  pars = c("dt_miss_else_mat","ii_miss_else","irr_ind_miss_e","irr_ind_miss_w"),
  include = FALSE,
  control = list(adapt_delta = .9)
)

saveRDS(fitmast,"data/4_model_out/mdl_fit.rds")


MCMCsummary(fitmast,params=c("psi1","psi2",
                             "beta1","beta2",
                             "lambda1","lambda2",
                             "zeta1","zeta2",
                             "kappa1","kappa2",
                             "mu_alpha",  
                             "sigma_alpha",
                             "mu_omega1","mu_omega2",
                             "sigma_omega1","sigma_omega2",
                             "gamma1","gamma2",
                             "mu_eta1","mu_eta2",
                             "sigma_eta1","sigma_eta2",
                             "theta","nu", "mast_total_mean",
                             "upsilon",
                             "sigma_phi"
),pg0 = TRUE,probs = c(0.055,0.5,0.945))

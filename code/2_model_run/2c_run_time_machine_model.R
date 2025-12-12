library(rstan)
library(MCMCvis)
library(shinystan)

#Blinded (1) or not blinded (0)?
blinded <- 0

for(yrs_back in 21:21){
  print(yrs_back)
  
  if (blinded == 0){
    data_mdl <- readRDS(paste("data/3_model_in/time_machine/not_blinded/time_machine_minus_",yrs_back,"_yrs.rds",sep=""))
  } else {
    data_mdl <- readRDS(paste("data/3_model_in/time_machine/blinded/time_machine_minus_blinded_",yrs_back,"_yrs.rds",sep=""))
  }
  
  fitmast <- stan(
    file = "code/2_model_run/2a_cmid.stan",  # Stan program
    data = data_mdl,    # named list of data
    chains = 4,             # number of Markov chains
    warmup = 1000,          # number of warmup iterations per chain
    iter = 2000,            # total number of iterations per chain
    cores = 4, # number of cores (could use one per chain)
    thin = 1, 
    refresh = 500,
    pars = c("dt_miss_else_mat","ii_miss_else","irr_ind_miss_e","irr_ind_miss_w","dt_else"),
    include = FALSE,
    #init = .1,
    control = list(adapt_delta = 0.9,
                   max_treedepth = 10#,
                   #stepsize=.005
    )
  )
  
  get_bfmi(fitmast)
  check_energy(fitmast)
  
  if (blinded == 0){
    saveRDS(fitmast,paste("data/4_model_out/time_machine/not_blinded/tm_",yrs_back,".rds",sep=""))
  } else {
    saveRDS(fitmast,paste("data/4_model_out/time_machine/blinded/tm_",yrs_back,".rds",sep=""))
  }
  fitmast <- NULL
}

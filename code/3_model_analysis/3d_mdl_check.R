#Check model fits - time machine
library(MCMCvis)
library(shinystan)
library(rstan)

for (nn in 21:1){
  print(nn)
  file_name <- paste("data/4_model_out/time_machine/not_blinded/tm_",nn,".rds",sep="")
  model_in <- readRDS(file_name)
  mdl_sum <- MCMCsummary(model_in,excl = c("pre_data_mast","phi",
                                           "dis_miss_e","dis_miss_w","yy_else",#"dt_else",
                                           "year_sum_x_fc_ef_else","year_sum_x_fc_wf_else",
                                           "mast_yr_means","phi","irr_e","irr_w",
                                           "disease_e","disease_w","Rho1","Rho2","lp__","mast_year_sd"))
  print(paste("Time Machine Model",nn))
  print(check_divergences(model_in))
  print("Energy warnings")
  print(check_energy(model_in))
  print("Rhat Warnings:")
  print(rownames(mdl_sum)[which(mdl_sum$Rhat > 1.03)])
  print("neff Warnings:")
  print(rownames(mdl_sum)[which(mdl_sum$n.eff < 400)])
  print("Problematic Params:")
  print(mdl_sum[unique(c(which(mdl_sum$n.eff < 400),which(mdl_sum$Rhat > 1.03))),])
}



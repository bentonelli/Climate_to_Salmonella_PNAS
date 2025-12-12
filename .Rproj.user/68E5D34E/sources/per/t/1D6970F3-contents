# Code to look at some effect sizes and check model params
library(dplyr)
library(MCMCvis)
library(ggplot2)
library(rstan)

mdl_data <- readRDS("data/3_model_in/model_full_data_2024.rds")
mdl_in <- readRDS("data/4_model_out/mdl_fit.rds")

mdl_sum <- MCMCsummary(mdl_in,excl = c("pre_data_mast","phi","yy_miss_else",
                                         "dis_miss_e","dis_miss_w","yy_else","dt_else",
                                         "year_sum_x_fc_ef_else","year_sum_x_fc_wf_else",
                                         "phi","irr_e","irr_w",
                                         "disease_e","disease_w","Rho1","Rho2","lp__"))

print(check_divergences(mdl_in))
print("Rhat Warnings:")
print(rownames(mdl_sum)[which(mdl_sum$Rhat > 1.03)])
print("neff Warnings:")
print(rownames(mdl_sum)[which(mdl_sum$n.eff < 400)])
print("Problematic Params:")
print(mdl_sum[unique(c(which(mdl_sum$n.eff < 400),which(mdl_sum$Rhat > 1.03))),])

MCMCsummary(mdl_in,params=c("psi1","psi2",
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


#Parameters for supplemental table
MCMCsummary(mdl_in,params=c("mu_omega1","omega1"),pg0 = TRUE,probs = c(0.055,0.5,0.945))
MCMCsummary(mdl_in,params=c("mu_eta1","eta1"),pg0 = TRUE,probs = c(0.055,0.5,0.945))

MCMCsummary(mdl_in,params=c("mu_omega2","omega2"),pg0 = TRUE,probs = c(0.055,0.5,0.945))
MCMCsummary(mdl_in,params=c("mu_eta2","eta2"),pg0 = TRUE,probs = c(0.055,0.5,0.945))
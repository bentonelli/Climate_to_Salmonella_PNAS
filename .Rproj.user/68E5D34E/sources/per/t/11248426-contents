# Posterior predictive checks for masting values
library(shinystan)
library(dplyr)
library(MCMCvis)
library(readr)
library(ggplot2)
library(MASS)

main_model_data <- readRDS("data/3_model_in/model_full_data_2024.rds")
main_model_fit <-  readRDS("data/4_model_out/mdl_fit.rds")

print("NS")

known_masting_indices <- which(main_model_data$yy_obs_mast!=.00001)
known_masting <- main_model_data$yy_obs_mast[which(main_model_data$yy_obs_mast!=.00001)]

yy_else_chains <- MCMCchains(main_model_fit,params = c("yy_else"))

par(mfrow=c(1,1))

plot(density(known_masting,bw=.5),main="NS/Mu_lambda",ylim=c(0,1),xlim=c(-.1,10))
for (n in 1:100){
  points(density(yy_else_chains[n,known_masting_indices],bw=.5),type="l",col=alpha("firebrick",.1))
}


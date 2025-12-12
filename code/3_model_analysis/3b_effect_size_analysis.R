# Code to look at some effect sizes
library(dplyr)
library(MCMCvis)
library(ggplot2)
library(boot)

mdl_data <- readRDS("data/3_model_in/model_full_data_2024.rds")
mdl_in <- readRDS("data/4_model_out/mdl_fit.rds")

MCMCsummary(mdl_in,params=c("omega1[1]","omega2[1]","eta1[1]","eta2[1]"),pg0 = TRUE,ISB=FALSE,probs = c(0.055,0.5,0.945))

MCMCplot(mdl_in,params=c("WRM"),horiz=FALSE)
MCMCplot(mdl_in,params=c("ERM"),horiz=FALSE)


### Cone Production and Delta-T ####
# First, look at the effect  on cone production of a X degree drop in temperatures 
# compared to a X degree increase.
## Model code reference: 
## exp(alpha[i] + varphi * yy_else[i-1] + theta * dt_else[i]+phi[i]*sigma_phi);
## yy_else[i] ~ gamma(upsilon,upsilon/mu_mast[i]);


summary(as.numeric(mdl_data$dt_obs))
sd(as.numeric(mdl_data$dt_obs))

upsilon_chains <- MCMCchains(mdl_in,params="upsilon")
mu_alpha_chains <- MCMCchains(mdl_in,params="mu_alpha")
nu_chains <- MCMCchains(mdl_in,params="nu")
theta_chains <- MCMCchains(mdl_in,params="theta")
yy_else_chains <- MCMCchains(mdl_in,params="yy_else")

avg_masting <- mean(yy_else_chains)
quant_masting <- quantile(yy_else_chains,c(0.05,.1,.5,.9,.95))

par(pty="s")
xx <- seq(-4,4,by=.1)
low_bound <- which(xx == -3)
upper_bound <- which(xx == 3)
mean_mast <- which(xx == 0)

#Randomly sample 1000 iterations
chain_sample <- 1:length(upsilon_chains)
plot(NULL,xlim=c(-5,5),ylim=c(0,2.5),ylab="Predicted Cone Production",xlab="Delta-T",cex.lab=1.6,cex.axis=1.6)
prop_mean_masting_pred_rec <- c()
for(nn in chain_sample){
  lp_iter <- mu_alpha_chains[nn] + nu_chains[nn] * avg_masting + theta_chains[nn] * xx
  
  # Mean of gamma is alpha/beta
  prop_mean_masting_pred <- upsilon_chains[nn]/(upsilon_chains[nn]/exp(lp_iter))
  prop_mean_masting_pred_rec <- rbind(prop_mean_masting_pred_rec,prop_mean_masting_pred)
  
}

mean_pred_response <- colMeans(prop_mean_masting_pred_rec)
upper_pred_response <- apply(prop_mean_masting_pred_rec,2,quantile,.945)
lower_pred_response <- apply(prop_mean_masting_pred_rec,2,quantile,.055)

points(xx,mean_pred_response,type="l",lwd=6,col="grey10")
polygon(c(xx,rev(xx)),c(lower_pred_response,rev(upper_pred_response)),
        col=alpha("grey10",.5),border=NA)


pred_drop_all <- (prop_mean_masting_pred_rec[,low_bound]/prop_mean_masting_pred_rec[,mean_mast])

quantile(pred_drop_all,c(0.055,.5,0.945))

### Look at autoregressive effects here ####
mean_dt <- mean(mdl_data$dt_obs)

par(pty="s")
xx <- seq(0.01,2.6,by=.01)
low_bound <- which(xx == .03)
upper_bound <- which(xx == 2.53)
mean_mast <- which(xx == 1)

chain_sample <- 1:length(upsilon_chains)
plot(NULL,xlim=c(0,3),ylim=c(0,3),ylab="Cone Production (t)",xlab="Cone Production (t-1)",cex.lab=1.6,cex.axis=1.6)
prop_mean_masting_pred_rec <- c()
for(nn in chain_sample){
  lp_iter <- mu_alpha_chains[nn] + nu_chains[nn] * xx + theta_chains[nn] * mean_dt
  
  # Mean of gamma is alpha/beta
  prop_mean_masting_pred <- upsilon_chains[nn]/(upsilon_chains[nn]/exp(lp_iter))
  prop_mean_masting_pred_rec <- rbind(prop_mean_masting_pred_rec,prop_mean_masting_pred)
  
}

mean_pred_response <- colMeans(prop_mean_masting_pred_rec)
upper_pred_response <- apply(prop_mean_masting_pred_rec,2,quantile,.945)
lower_pred_response <- apply(prop_mean_masting_pred_rec,2,quantile,.055)

points(xx,mean_pred_response,type="l",lwd=6,col="grey10")
polygon(c(xx,rev(xx)),c(lower_pred_response,rev(upper_pred_response)),
        col=alpha("grey10",.5),border=NA)


pred_drop_all <- 1 - (prop_mean_masting_pred_rec[,upper_bound]/prop_mean_masting_pred_rec[,mean_mast])

quantile(pred_drop_all,c(0.055,.5,0.945))

### Regional Masting Trends ####

# Here, calculate the variance in interannual regional cone production trends
par(pty="m")
par(mfrow=c(2,1))
par(mar=c(4,6,4,4))

mast_yrs <- mdl_data$mast_years
west_regional_masting_chains <- MCMCchains(mdl_in,params="WRM")

med_mast_est_wrm <- apply(west_regional_masting_chains,2,quantile,.5)
upp_mast_est_wrm <- apply(west_regional_masting_chains,2,quantile,.945)
low_mast_est_wrm <- apply(west_regional_masting_chains,2,quantile,.055)

plot(NULL,xlim=c(1982,2024),ylim=c(-1,2),ylab="Cone Production",xlab="Year",cex.lab=1.6,cex.axis=1.6)
points(mast_yrs,med_mast_est_wrm,type="l",col="#648DE5",lwd=6)
polygon(c(mast_yrs,rev(mast_yrs)),c(low_mast_est_wrm,rev(upp_mast_est_wrm)),
        col=alpha("#648DE5",.5),border=NA)

east_regional_masting_chains <- MCMCchains(mdl_in,params="ERM")

med_mast_est_erm <- apply(east_regional_masting_chains,2,quantile,.5)
upp_mast_est_erm <- apply(east_regional_masting_chains,2,quantile,.945)
low_mast_est_erm <- apply(east_regional_masting_chains,2,quantile,.055)

plot(NULL,xlim=c(1982,2024),ylim=c(-1,2),ylab="Cone Production",xlab="Year",cex.lab=1.6,cex.axis=1.6)
points(mast_yrs,med_mast_est_erm,type="l",col="#E98A15",lwd=6)
polygon(c(mast_yrs,rev(mast_yrs)),c(low_mast_est_erm,rev(upp_mast_est_erm)),
        col=alpha("#E98A15",.5),border=NA)

# Get year to year variance from one year to next
interannual_diff <- function(chain_in){
  int_dif <- abs(chain_in[2:length(chain_in)] - chain_in[1:(length(chain_in)-1)])
  return(mean(int_dif))
}


west_est_diff <- apply(west_regional_masting_chains,1,interannual_diff)
mean(west_est_diff)

east_est_diff <- apply(east_regional_masting_chains,1,interannual_diff)
mean(east_est_diff)


### Effect of irruptions on outbreaks ####

#Reference:
#outbreak_yn_west[i] ~ bernoulli_logit(psi1 + beta1 * irr_w[i,1]);
#outbreak_yn_east[i] ~ bernoulli_logit(psi2 + beta2 * irr_e[i,1]);

#if (outbreak_yn_west[i] == 1){
#  disease_w[i] ~ gamma(lambda1,(lambda1/exp(zeta1 + kappa1 * irr_w[i,1])));
#}

#if (outbreak_yn_east[i] == 1){
#  //Model outbreak liklihood
#  disease_e[i] ~ gamma(lambda2,(lambda2/exp(zeta2 + kappa2 * irr_e[i,1])));
#}

sd_PISI_w <- sd(mdl_data$irr_ind_w[,1])*1.5 + mean(mdl_data$irr_ind_w[,1])
sd_PISI_e <- sd(mdl_data$irr_ind_e[,1])*1.5 + mean(mdl_data$irr_ind_e[,1])

par(pty="s")
irr_xx <- seq(-2,2,by=.1)

psi1_chains <- MCMCchains(mdl_in,params = "psi1")
beta1_chains <- MCMCchains(mdl_in,params = "beta1")

lambda1_chains <- MCMCchains(mdl_in,params = "lambda1")
zeta1_chains <- MCMCchains(mdl_in,params = "zeta1")
kappa1_chains <- MCMCchains(mdl_in,params = "kappa1")

#Set seed for replicability
set.seed(111)
rec_ob_low <- c()
rec_ob_hi <- c()
chain_sample <- 1:length(lambda1_chains)
yn_trial_rec <- c()
for (nn in chain_sample){
  outbreak_yn_west_low <- inv.logit(psi1_chains[nn] + beta1_chains[nn] * -sd_PISI_w)
  outbreak_yn_west_hi <- inv.logit(psi1_chains[nn] + beta1_chains[nn] * sd_PISI_w)

  yn_trial_low <- rbinom(1,size=1,prob = outbreak_yn_west_low)
  yn_trial_hi <- rbinom(1,size=1,prob = outbreak_yn_west_hi)
  
  yn_trial_rec <- rbind(yn_trial_rec,c(yn_trial_low,yn_trial_hi))
  
  if (yn_trial_low == 1){
    ob_low <- exp(yn_trial_low * rgamma(1,shape = lambda1_chains[nn],rate = lambda1_chains[nn]/exp(zeta1_chains[nn]+kappa1_chains[nn]*(-1*sd_PISI_w))))
  } else {
    ob_low = 0
  }
  
  if(yn_trial_hi == 1){
    ob_hi <- exp(yn_trial_hi * rgamma(1,shape = lambda1_chains[nn],rate = lambda1_chains[nn]/exp(zeta1_chains[nn]+kappa1_chains[nn]*(1*sd_PISI_w))))
  } else {
    ob_hi <- 0
  }
  
  rec_ob_low <- c(rec_ob_low,ob_low)
  rec_ob_hi <- c(rec_ob_hi,ob_hi)
}
print("Western")
colSums(yn_trial_rec)/nrow(yn_trial_rec)

sum(rec_ob_low>=500)/length(chain_sample)
sum(rec_ob_hi>=500)/length(chain_sample)

#East
psi2_chains <- MCMCchains(mdl_in,params = "psi2")
beta2_chains <- MCMCchains(mdl_in,params = "beta2")

lambda2_chains <- MCMCchains(mdl_in,params = "lambda2")
zeta2_chains <- MCMCchains(mdl_in,params = "zeta2")
kappa2_chains <- MCMCchains(mdl_in,params = "kappa2")

rec_ob_low <- c()
rec_ob_hi <- c()
chain_sample <- 1:length(lambda2_chains)
yn_trial_rec <- c()
for (nn in chain_sample){
  outbreak_yn_west_low <- inv.logit(psi2_chains[nn] + beta2_chains[nn] * -sd_PISI_w)
  outbreak_yn_west_hi <- inv.logit(psi2_chains[nn] + beta2_chains[nn] * sd_PISI_w)
  
  yn_trial_low <- rbinom(1,size=1,prob = outbreak_yn_west_low)
  yn_trial_hi <- rbinom(1,size=1,prob = outbreak_yn_west_hi)
  
  yn_trial_rec <- rbind(yn_trial_rec,c(yn_trial_low,yn_trial_hi))
  
  if (yn_trial_low == 1){
    ob_low <- exp(yn_trial_low * rgamma(1,shape = lambda2_chains[nn],rate = lambda2_chains[nn]/exp(zeta2_chains[nn]+kappa2_chains[nn]*(-1*sd_PISI_e))))
  } else {
    ob_low = 0
  }
  
  if(yn_trial_hi == 1){
    ob_hi <- exp(yn_trial_hi * rgamma(1,shape = lambda2_chains[nn],rate = lambda2_chains[nn]/exp(zeta2_chains[nn]+kappa2_chains[nn]*(1*sd_PISI_e))))
  } else {
    ob_hi <- 0
  }
  
  rec_ob_low <- c(rec_ob_low,ob_low)
  rec_ob_hi <- c(rec_ob_hi,ob_hi)
}
print("Eastern")
colSums(yn_trial_rec)/nrow(yn_trial_rec)
sum(rec_ob_low>=500)/length(chain_sample)
sum(rec_ob_hi>=500)/length(chain_sample)


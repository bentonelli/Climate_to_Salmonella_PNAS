#Script to plot regional masting estimates versus irruptions

library(MCMCvis)
library(ggplot2)
library(dplyr)
library(boot)

mdl_data <- readRDS("data/3_model_in/model_full_data_2024.rds")
mdl_in <- readRDS("data/4_model_out/mdl_fit.rds")


MCMCsummary(mdl_in,params=c("psi1","psi2",
                             "beta1","beta2",
                             "lambda1","lambda2",
                             "zeta1","zeta2",
                             "kappa1","kappa2",
                             "mu_alpha1","mu_alpha2",    
                             "sigma_alpha1","sigma_alpha2",
                             "mu_omega1","mu_omega2",
                             "sigma_omega1","sigma_omega2",
                             "mu_gamma1","mu_gamma2",
                             "mu_eta1","mu_eta2",
                             "sigma_eta1","sigma_eta2",
                             "theta","varphi", "varphi_total_mean",
                             "upsilon","tau","nu","sigma_dt",
                             "sigma_phi1","sigma_phi2"
),pg0 = TRUE,probs = c(0.055,0.5,0.945))

### Cone Production x Irruption Strength ####

# Reference

masting_w <- MCMCchains(mdl_in,params="WRM")
masting_e <- MCMCchains(mdl_in,params="ERM")

#West, subplot A
quants_w <- quantile(as.numeric(masting_w),c(.055,.945))
xx_w <- seq(quants_w[1],quants_w[2],by=.1)
#mu_gamma1 + omega1 * WRM[i] + eta1 .* pre_data_irr_w

mu_gamma1_chains <- MCMCchains(mdl_in,"gamma1")
omega1_chains <- MCMCchains(mdl_in,"omega1[1]",ISB=FALSE)

par(las=1,pty="s")
plot(NULL, xlim=c(-1,1.5),ylim=c(-3,3),cex.lab=1.8,cex.axis=1.5,ylab="",xlab="")
points(colMeans(masting_w)[mdl_data$irr_obs_yrs],mdl_data$irr_ind_w[mdl_data$irr_obs_yrs,1],
       col=alpha("#648DE5",.85),pch=19,cex=1.75)

rand_samp <- sample(length(mu_gamma1_chains),1000)
samp_eff_rec <- c()
for (nn in rand_samp){
  samp_eff <- mu_gamma1_chains[nn] + omega1_chains[nn] * xx_w
  samp_eff_rec <- rbind(samp_eff_rec,samp_eff)
}

mean_pred_response <- colMeans(samp_eff_rec)
upper_pred_response <- apply(samp_eff_rec,2,quantile,.945)
lower_pred_response <- apply(samp_eff_rec,2,quantile,.055)

points(xx_w,mean_pred_response,type="l",lwd=6,col="#648DE5")
polygon(c(xx_w,rev(xx_w)),c(lower_pred_response,rev(upper_pred_response)),
        col=alpha("#648DE5",.5),border=NA)

irr_all_w <- mdl_data$irr_ind_w[mdl_data$irr_obs_yrs,]
irr_all_e <- mdl_data$irr_ind_e[mdl_data$irr_obs_yrs,]


#East, subplot B
quants_e <- quantile(as.numeric(masting_e),c(.055,.945))
xx_e <- seq(quants_e[1],quants_e[2],by=.1)
#mu_gamma1 + omega1 * WRM[i] + eta1 .* pre_data_irr_w

mu_gamma2_chains <- MCMCchains(mdl_in,"gamma2")
omega2_chains <- MCMCchains(mdl_in,"omega2[1]",ISB=FALSE)

par(pty="s")
plot(NULL, xlim=c(-1,1.5),ylim=c(-3,3),cex.lab=1.8,cex.axis=1.5,ylab="",xlab="")
points(colMeans(masting_e)[mdl_data$irr_obs_yrs],mdl_data$irr_ind_e[mdl_data$irr_obs_yrs,1],
       col=alpha("#E98A15",.85),pch=19,cex=1.75)

rand_samp <- sample(length(mu_gamma2_chains),1000)
samp_eff_rec <- c()
for (nn in rand_samp){
  samp_eff <- mu_gamma2_chains[nn] + omega2_chains[nn] * xx_e
  samp_eff_rec <- rbind(samp_eff_rec,samp_eff)
}

mean_pred_response <- colMeans(samp_eff_rec)
upper_pred_response <- apply(samp_eff_rec,2,quantile,.945)
lower_pred_response <- apply(samp_eff_rec,2,quantile,.055)

points(xx_e,mean_pred_response,type="l",lwd=6,col="#E98A15")
polygon(c(xx_e,rev(xx_e)),c(lower_pred_response,rev(upper_pred_response)),
        col=alpha("#E98A15",.5),border=NA)

### Irruption Strength X Outbreak Probability #### 

#Reference
# //Disease
# if (dis_obs_vect[i] == 1){
#   
#   outbreak_yn_west[i] ~ bernoulli_logit(psi1 + rho1 * irr_w[i,1]);
#   outbreak_yn_east[i] ~ bernoulli_logit(psi2 + rho2 * irr_e[i,1]);
#   
#   if (outbreak_yn_west[i] == 1){
#     disease_w[i] ~ gamma(lambda1,(lambda1/exp(zeta1 + kappa1 * irr_w[i,1])));
#   }
#   
#   if (outbreak_yn_east[i] == 1){
#     //Model outbreak liklihood
#     disease_e[i] ~ gamma(lambda2,(lambda2/exp(zeta2 + kappa2 * irr_e[i,1])));
#   }
# }

#First plot, irruptions in the west compared to disease counts

irr_siskin_e <- mdl_data$irr_ind_e[mdl_data$irr_obs_yrs,1]
irr_siskin_w <- mdl_data$irr_ind_w[mdl_data$irr_obs_yrs,1]

dis_e <- mdl_data$east_disease_ts[mdl_data$disease_ts_observed]
dis_w <- mdl_data$west_disease_ts[mdl_data$disease_ts_observed]


# West - Plot C
xx_w <- seq(min(irr_siskin_w[mdl_data$disease_ts_observed]),max(irr_siskin_w[mdl_data$disease_ts_observed]),by=.1)

w_param_chains <- MCMCchains(mdl_in,params = c("zeta1","kappa1","lambda1","psi1","beta1")) %>% as.data.frame()

rand_bern_draws <- c()
for (n in rand_samp){
  lp_bern <- inv.logit(w_param_chains$psi1[n] + w_param_chains$beta1[n] * xx_w)
  rand_bern_draws <- rbind(rand_bern_draws,lp_bern)
}
bern_means <- colMeans(rand_bern_draws)

par(las=1,pty="s")
#Set no outbreak years as black dots, others as colored
col_fact <- as.numeric(as.factor(dis_w == 0))
col_options <- c("#648DE5","#648DE5")
pch_options <- c(19,1)

dis_w[which(dis_w != 0)] <- 1
plot(irr_siskin_w[mdl_data$disease_ts_observed],dis_w,pch=pch_options[col_fact],col=alpha("#648DE5",.8),
     ylab="",xlab="",cex.lab=1.8,cex.axis=1.5,cex=1.75,xlim=c(-3,3))
points(xx_w,bern_means, ylim=c(0,1),col=alpha("#648DE5",1),
       type = "l",lty=1,lwd=4,xlab="Irruption Strength",ylab="Outbreak Detection Lilklihood",
       cex.lab=1.6,cex.axis=1.25)
quants_bern <- apply(rand_bern_draws,2,quantile,c(.055,.945))
polygon(c(xx_w,rev(xx_w)),c(quants_bern[1,],rev(quants_bern[2,])),
        col=alpha("#648DE5",.25),border=NA)

#East - Plot D
xx_e <- seq(min(irr_siskin_e[mdl_data$disease_ts_observed]),max(irr_siskin_e[mdl_data$disease_ts_observed]),by=.1)

e_param_chains <- MCMCchains(mdl_in,params = c("zeta2","kappa2","lambda2","psi2","beta2")) %>% as.data.frame()

rand_bern_draws <- c()
for (n in rand_samp){
  lp_bern <- inv.logit(e_param_chains$psi2[n] + e_param_chains$beta2[n] * xx_e)
  rand_bern_draws <- rbind(rand_bern_draws,lp_bern)
}
bern_means <- colMeans(rand_bern_draws)

par(las=1,pty="s")
#Set no outbreak years as black dots, others as colored
col_fact <- as.numeric(as.factor(dis_e == 0))
col_options <- c("#E98A15","#E98A15")
pch_options <- c(19,1)

dis_e[which(dis_e != 0)] <- 1
plot(irr_siskin_e[mdl_data$disease_ts_observed],dis_e,pch=pch_options[col_fact],col=alpha("#E98A15",.8),
     ylab="",xlab="",cex.lab=1.8,cex.axis=1.5,cex=1.75,xlim=c(-3,3))
points(xx_e,bern_means, ylim=c(0,1),col=alpha("#E98A15",1),
       type = "l",lty=1,lwd=4,xlab="Irruption Strength",ylab="Outbreak Detection Lilklihood",
       cex.lab=1.6,cex.axis=1.25)
quants_bern <- apply(rand_bern_draws,2,quantile,c(.055,.945))
polygon(c(xx_e,rev(xx_e)),c(quants_bern[1,],rev(quants_bern[2,])),
        col=alpha("#E98A15",.25),border=NA)


### Irruptions Strength x Outbreak Size ####

# Reference: disease_w[i] ~ gamma(lambda1,(lambda1/exp(zeta1 + kappa1 * irr_w[i,1])));
#West - Plot E
rand_outbreak_size<- c()
for (n in 1:nrow(w_param_chains)){
  outbreak_size <- exp(w_param_chains$zeta1[n] + w_param_chains$kappa1[n]* xx_w)
  rand_outbreak_size <- rbind(rand_outbreak_size,w_param_chains$lambda1[n]/(w_param_chains$lambda1[n]/outbreak_size))
}

outbreak_size_mean <- colMeans(rand_outbreak_size)
quants <- apply(rand_outbreak_size,2,quantile,c(.055,.945))

plot(xx_w,outbreak_size_mean,type="l",col="#648DE5",lwd=5,ylim=c(0,10),
     xlim=c(-3,3),
     xlab="",ylab="",yaxt="n",
     cex.lab=1.8,cex.axis=1.5)
polygon(c(xx_w,rev(xx_w)),c(quants[1,],rev(quants[2,])),
        col=alpha("#648DE5",.5),border=NA)

#Plot points when disease data is known and > 0 
w_irr_dis <- data.frame(w_pisi_irr = irr_siskin_w[mdl_data$disease_ts_observed],
           w_dis = mdl_data$west_disease_ts[mdl_data$disease_ts_observed])
w_irr_dis <- w_irr_dis[which(w_irr_dis$w_dis != 0),]

points(w_irr_dis$w_pisi_irr,w_irr_dis$w_dis,col=alpha("#648DE5",.8),cex=1.75,pch=19)

y_log_scale <- c(0,10,100,1000,10000)+1
axis(side=2, at = log(y_log_scale),labels = y_log_scale-1,cex.lab=1.8,cex.axis=1.5)


#East - Plot F
rand_outbreak_size<- c()
for (n in 1:nrow(e_param_chains)){
  outbreak_size <- exp(e_param_chains$zeta2[n] + e_param_chains$kappa2[n]* xx_e)
  rand_outbreak_size <- rbind(rand_outbreak_size,e_param_chains$lambda2[n]/(e_param_chains$lambda2[n]/outbreak_size))
}

outbreak_size_mean <- colMeans(rand_outbreak_size)
quants <- apply(rand_outbreak_size,2,quantile,c(.055,.945))

plot(xx_e,outbreak_size_mean,type="l",col="#E98A15",lwd=5,ylim=c(0,10),
     xlim=c(-3,3),
     xlab="",ylab="",yaxt="n",
     cex.lab=1.8,cex.axis=1.5)
polygon(c(xx_e,rev(xx_e)),c(quants[1,],rev(quants[2,])),
        col=alpha("#E98A15",.5),border=NA)

#Plot points when disease data is known and > 0 
e_irr_dis <- data.frame(e_pisi_irr = irr_siskin_e[mdl_data$disease_ts_observed],
                        e_dis = mdl_data$east_disease_ts[mdl_data$disease_ts_observed])
e_irr_dis <- e_irr_dis[which(e_irr_dis$e_dis != 0),]

points(e_irr_dis$e_pisi_irr,e_irr_dis$e_dis,col=alpha("#E98A15",.8),cex=1.75,pch=19)

y_log_scale <- c(0,10,100,1000,10000)+1
axis(side=2, at = log(y_log_scale),labels = y_log_scale-1,cex.lab=1.8,cex.axis=1.5)


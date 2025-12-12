# Code to analyze predictions

library(MCMCvis)
library(dplyr)
library(scoringRules)
library(ggplot2)

full_data <- readRDS("data/3_model_in/model_full_data_2024.rds")

blinded_yn <- 0
final_year <- 2024

irr_pred_df <- data.frame()
dis_pred_df <- data.frame()
#dev.off()
if (blinded_yn == 1){
  pdf("data/5_predictive_analysis/pred_analysis_blinded.pdf")  
} else {
  pdf("data/5_predictive_analysis/pred_analysis_not_blinded.pdf")
}

for (nn in 21:1){
  par(mfrow=c(3,2))
  # Read in time machine model fit
  yrs_back <- nn
  
  if (blinded_yn == 1){
    pred_model <- readRDS(paste("data/4_model_out/time_machine/blinded/tm_",yrs_back,".rds",sep=""))
    pred_data <- readRDS(paste("data/3_model_in/time_machine/blinded/time_machine_minus_blinded_",yrs_back,"_yrs.rds",sep=""))
  } else {
    pred_model <- readRDS(paste("data/4_model_out/time_machine/not_blinded/tm_",yrs_back,".rds",sep=""))
    pred_data <- readRDS(paste("data/3_model_in/time_machine/not_blinded/time_machine_minus_",yrs_back,"_yrs.rds",sep=""))
  }
  
  
  if (nn == 1){
    next_year_data <- full_data
  } else {
    if (blinded_yn == 1){
      next_year_data <- readRDS(paste("data/3_model_in/time_machine/blinded/time_machine_minus_blinded_",yrs_back-1,"_yrs.rds",sep=""))
    } else {
      next_year_data <- readRDS(paste("data/3_model_in/time_machine/not_blinded/time_machine_minus_",yrs_back-1,"_yrs.rds",sep=""))
    }
  }
  
  ### Look at regional masting estimates ####
  MCMCplot(pred_model,params="WRM",horiz=FALSE)
  MCMCplot(pred_model,params="ERM",horiz=FALSE)
  
  #Save number of observed masting datapoints, 
  yr0_masting_cells <- pred_data$N_obs_mast[length(pred_data$N_obs_mast)]
  yr1_masting_cells <- pred_data$N_obs_mast[length(pred_data$N_obs_mast)-1]
  
  #### Look at irruption predictions for East and West ####
  irr_yr_length <- nrow(pred_data$irr_ind_w)
  p_ind_irr <- rep(1:8,irr_yr_length)
  
  irr_w_chains <- MCMCchains(pred_model,params="irr_w")
  irr_e_chains <- MCMCchains(pred_model,params="irr_e")
  
  irr_full_w <- as.numeric(unlist(pred_data$irr_ind_w[1:(nrow(pred_data$irr_ind_w)-1),1]))
  irr_full_e <- as.numeric(unlist(pred_data$irr_ind_e[1:(nrow(pred_data$irr_ind_e)-1),1]))
  
  pisi_w <- irr_w_chains[,p_ind_irr==1]
  pisi_e <- irr_e_chains[,p_ind_irr==1]
  
  pisi_w_pred <- pisi_w[,ncol(pisi_w)]
  pisi_e_pred <- pisi_e[,ncol(pisi_e)]
  
  pisi_historical_baseline_w <- irr_full_w
  pisi_historical_baseline_e <- irr_full_e
  
  actual_pisi_w <- as.numeric(unlist(next_year_data$irr_ind_w[nrow(next_year_data$irr_ind_w)-1,1]))
  actual_pisi_e <- as.numeric(unlist(next_year_data$irr_ind_e[nrow(next_year_data$irr_ind_w)-1,1]))
  
  pisi_w_score_cascade <- crps_sample(actual_pisi_w,pisi_w_pred)
  pisi_e_score_cascade <- crps_sample(actual_pisi_e,pisi_e_pred)
  
  pisi_w_score_historical <- crps_sample(actual_pisi_w,pisi_historical_baseline_w)
  pisi_e_score_historical <- crps_sample(actual_pisi_e,pisi_historical_baseline_e)
  
  print("West")
  print("Historical Cascade:")
  print(pisi_w_score_cascade)
  print("Historical Baseline:")
  print(pisi_w_score_historical)
  
  plot(density(pisi_historical_baseline_w,bw=.25),
       col="#648DE5",lwd=4,lty=3,xlim=c(-5,5),
       main=paste("Historical Prediction, West - ",final_year-yrs_back),
       xlab="Irruption Strength")
  points(density(pisi_w_pred,bw=.25),type="l",col="#648DE5",lwd=4,lty=1)
  abline(v=actual_pisi_w)
  abline(v=pisi_w[1,ncol(pisi_w)-1],lty=3)
  
  
  print("East")
  print("Historical Cascade:")
  print(pisi_e_score_cascade)
  print("Historical Baseline:")
  print(pisi_e_score_historical)
  
  plot(density(pisi_historical_baseline_e,bw=.25),
       col="#E98A15",lwd=4,lty=3,xlim=c(-5,5),
       main=paste("Historical Prediction, East - ",final_year-yrs_back),
       xlab="Irruption Strength")
  points(density(pisi_e_pred,bw=.25),type="l",col="#E98A15",lwd=4,lty=1)
  abline(v=actual_pisi_e)
  abline(v=pisi_e[1,ncol(pisi_e)-1],lty=3)
  
  # Save scores as df
  irr_pred_df <- rbind(irr_pred_df,c(final_year-yrs_back,"West",yr0_masting_cells,yr1_masting_cells,actual_pisi_w,pisi_w_score_cascade,pisi_w_score_historical))
  irr_pred_df <- rbind(irr_pred_df,c(final_year-yrs_back,"East",yr0_masting_cells,yr1_masting_cells,actual_pisi_e,pisi_e_score_cascade,pisi_e_score_historical))
  
  ### Get disease predictions #### 
  
  MCMCsummary(pred_model,params = c("rho1","kappa1","rho2","kappa2"))
  
  dis_w_chains <- MCMCchains(pred_model,params="dis_w_pred")
  dis_e_chains <- MCMCchains(pred_model,params="dis_e_pred")
  
  dis_w_pred <- dis_w_chains[,ncol(dis_w_chains)]
  dis_e_pred <- dis_e_chains[,ncol(dis_e_chains)]
  
  #dis_historical_baseline_w <- sample(pred_data$west_disease_ts[pred_data$disease_ts_observed],length(dis_w_pred),replace=TRUE)
  dis_historical_baseline_w <- pred_data$west_disease_ts[pred_data$disease_ts_observed]
  #dis_historical_baseline_e <- sample(pred_data$east_disease_ts[pred_data$disease_ts_observed],length(dis_e_pred),replace=TRUE)
  dis_historical_baseline_e <- pred_data$east_disease_ts[pred_data$disease_ts_observed]
  
  dis_act_w <- full_data$west_disease_ts[length(pred_data$west_disease_ts)]
  dis_act_e <- full_data$east_disease_ts[length(pred_data$east_disease_ts)]
  
  hist(dis_historical_baseline_w,freq=FALSE,xlim=c(-1,12),breaks = 0:50,col=alpha("black",.5),ylim=c(0,.5))
  hist(dis_w_pred,freq=FALSE,add=T,breaks = 0:50,col=alpha("#648DE5",.5))
  abline(v=dis_act_w,lwd=4)
  
  hist(dis_historical_baseline_e,freq=FALSE,xlim=c(-1,12),breaks = 0:50,col=alpha("black",.5),ylim=c(0,.5))
  hist(dis_e_pred,freq=FALSE,add=T,breaks = 0:50,col=alpha("#E98A15",.5))
  abline(v=dis_act_e,lwd=4)
  
  dis_w_score_cascade <- crps_sample(dis_act_w,dis_w_pred)
  dis_e_score_cascade <- crps_sample(dis_act_e,dis_e_pred)
  
  dis_w_score_historical <- crps_sample(dis_act_w,dis_historical_baseline_w)
  dis_e_score_historical <- crps_sample(dis_act_e,dis_historical_baseline_e)
  
  print("West")
  print("Historical Cascade:")
  print(dis_w_score_cascade)
  print("Historical Baseline:")
  print(dis_w_score_historical)
  
  print("East")
  print("Historical Cascade:")
  print(dis_e_score_cascade)
  print("Historical Baseline:")
  print(dis_e_score_historical)
  
  #Get the chance of a large outbreak (>=100 individuals)
  historical_pred_large_outbreak_w <- sum((exp(dis_historical_baseline_w)-1)>=100)/length(dis_historical_baseline_w)
  cascade_pred_large_outbreak_w <- sum((exp(dis_w_pred)-1)>=100)/length(dis_w_pred)
  
  historical_pred_large_outbreak_e <- sum((exp(dis_historical_baseline_e)-1)>=100)/length(dis_historical_baseline_e)
  cascade_pred_large_outbreak_e <- sum((exp(dis_e_pred)-1>=100))/length(dis_e_pred)
  
  dis_pred_df <- rbind(dis_pred_df,c(final_year-yrs_back,"West",yr0_masting_cells,
                                     yr1_masting_cells,dis_act_w,dis_w_score_cascade,
                                     dis_w_score_historical,mean(dis_historical_baseline_w),
                                     mean(dis_w_pred),
                                     historical_pred_large_outbreak_w,cascade_pred_large_outbreak_w))
  dis_pred_df <- rbind(dis_pred_df,c(final_year-yrs_back,"East",yr0_masting_cells,
                                     yr1_masting_cells,dis_act_e,dis_e_score_cascade,
                                     dis_e_score_historical,mean(dis_historical_baseline_e),
                                     mean(dis_e_pred),
                                     historical_pred_large_outbreak_e,cascade_pred_large_outbreak_e))
}
dev.off()
colnames(irr_pred_df) <- c("Year","Region","Mast_t0","Mast_tm1","Actual_Irr","crps_Cascade","crps_Historical")
colnames(dis_pred_df) <- c("Year","Region","Mast_t0","Mast_tm1","Actual_Dis","crps_Cascade","crps_Historical","mean_historical","mean_pred","large_epi_historical","lg_epi_mdl")

irr_pred_df$Year <- as.numeric(irr_pred_df$Year)
irr_pred_df$Actual_Irr <- as.numeric(irr_pred_df$Actual_Irr)
irr_pred_df$crps_Cascade <- as.numeric(irr_pred_df$crps_Cascade)
irr_pred_df$crps_Historical <- as.numeric(irr_pred_df$crps_Historical)

dis_pred_df$Year <- as.numeric(dis_pred_df$Year)
dis_pred_df$Actual_Dis <- as.numeric(dis_pred_df$Actual_Dis)
dis_pred_df$crps_Cascade <- as.numeric(dis_pred_df$crps_Cascade)
dis_pred_df$crps_Historical <- as.numeric(dis_pred_df$crps_Historical)


plot(irr_pred_df$Actual_Irr,irr_pred_df$crps_Cascade-irr_pred_df$crps_Historical)
abline(h=0)

1 - (sum(irr_pred_df$crps_Cascade)/sum(irr_pred_df$crps_Historical))
1 - (sum(dis_pred_df$crps_Cascade)/sum(dis_pred_df$crps_Historical))

plot(dis_pred_df$Actual_Dis,dis_pred_df$crps_Cascade-dis_pred_df$crps_Historical)
abline(h=0)
dev.off()

if(blinded_yn == 1){
  saveRDS(irr_pred_df,"data/5_predictive_analysis/irr_prediction_scoring_blinded.rds")
  saveRDS(dis_pred_df,"data/5_predictive_analysis/dis_prediction_scoring_blinded.rds")  
} else {
  saveRDS(irr_pred_df,"data/5_predictive_analysis/irr_prediction_scoring_not_blinded.rds")
  saveRDS(dis_pred_df,"data/5_predictive_analysis/dis_prediction_scoring_not_blinded.rds")
}


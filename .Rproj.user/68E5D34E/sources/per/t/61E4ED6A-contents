# Code to analyze predictions

library(MCMCvis)
library(dplyr)
library(scoringRules)
library(ggplot2)

full_data <- readRDS("data/3_model_in/model_full_data_2024.rds")

blinded_yn <- 0
final_year <- 2024
brier_score_record <- c()
par(mfrow=c(3,2))
for (nn in 21:1){
  print(nn)
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
  
  
  dis_w_chains <- MCMCchains(pred_model,params="dis_w_pred")
  dis_e_chains <- MCMCchains(pred_model,params="dis_e_pred")
  
  dis_w_pred <- dis_w_chains[,ncol(dis_w_chains)]
  dis_e_pred <- dis_e_chains[,ncol(dis_e_chains)]
  
  dis_historical_baseline_w <- pred_data$west_disease_ts[pred_data$disease_ts_observed]
  dis_historical_baseline_e <- pred_data$east_disease_ts[pred_data$disease_ts_observed]
  
  dis_act_w <- full_data$west_disease_ts[length(pred_data$west_disease_ts)]
  dis_act_e <- full_data$east_disease_ts[length(pred_data$east_disease_ts)]
  
  # WEST
  #Any outbreak
  dis_w_chance_binary_hist <- sum(dis_historical_baseline_w!=0)/length(dis_historical_baseline_w)
  dis_w_chance_binary_cascade <- sum(dis_w_pred!=0)/length(dis_w_pred)
  
  dis_act_w_binary <- as.numeric(dis_act_w>0)
  
  dis_w_brier_hist <- (dis_w_chance_binary_hist - dis_act_w_binary)^2
  dis_w_brier_cascade <- (dis_w_chance_binary_cascade - dis_act_w_binary)^2
  
  row_add <- c(nn,0,"West",dis_act_w_binary,dis_w_brier_hist,dis_w_brier_cascade)
  brier_score_record <- rbind(brier_score_record,row_add)
  
  # EAST 
  #Any outbreak
  dis_act_e_binary <- as.numeric(dis_act_e>0)
  
  dis_e_chance_binary_hist <- sum(dis_historical_baseline_e!=0)/length(dis_historical_baseline_e)
  dis_e_chance_binary_cascade <- sum(dis_e_pred!=0)/length(dis_e_pred)
  
  dis_e_brier_hist <- (dis_e_chance_binary_hist - dis_act_e_binary)^2
  dis_e_brier_cascade <- (dis_e_chance_binary_cascade - dis_act_e_binary)^2
  
  row_add <- c(nn,0,"East",dis_act_e_binary,dis_e_brier_hist,dis_e_brier_cascade)
  brier_score_record <- rbind(brier_score_record,row_add)
  
  # Big outbreak - sizes ultimately arbitrary for scoring purposes. CRPSS is a better
  # metric and is reported in the paper. 
  for(large_outbreak_size in c(100,500)){
    # West
    dis_act_w_binary_big <- as.numeric(dis_act_w>(log(large_outbreak_size+1)))
    
    dis_w_chance_binary_hist_big <- sum(dis_historical_baseline_w>log(large_outbreak_size+1))/length(dis_historical_baseline_w)
    dis_w_chance_binary_cascade_big <- sum(dis_w_pred>log(large_outbreak_size+1))/length(dis_w_pred)
    
    dis_w_brier_hist_big <- (dis_w_chance_binary_hist_big - dis_act_w_binary_big)^2
    dis_w_brier_cascade_big <- (dis_w_chance_binary_cascade_big - dis_act_w_binary_big)^2
    
    row_add <- c(nn,large_outbreak_size,"West",dis_act_w_binary_big,dis_w_brier_hist_big,dis_w_brier_cascade_big)
    brier_score_record <- rbind(brier_score_record,row_add)

    #East
    dis_act_e_binary_big <- as.numeric(dis_act_e>(log(large_outbreak_size+1)))
    
    dis_e_chance_binary_hist_big <- sum(dis_historical_baseline_e>log(large_outbreak_size+1))/length(dis_historical_baseline_e)
    dis_e_chance_binary_cascade_big <- sum(dis_e_pred>log(large_outbreak_size+1))/length(dis_e_pred)
    
    dis_e_brier_hist_big <- (dis_e_chance_binary_hist_big - dis_act_e_binary_big)^2
    dis_e_brier_cascade_big <- (dis_e_chance_binary_cascade_big - dis_act_e_binary_big)^2
    
    row_add <- c(nn,large_outbreak_size,"East",dis_act_e_binary_big,dis_e_brier_hist_big,dis_e_brier_cascade_big)
    brier_score_record <- rbind(brier_score_record,row_add)
  }
  
  
  hist(dis_historical_baseline_w,freq=FALSE,xlim=c(-1,12),breaks = 0:50,col=alpha("black",.5),ylim=c(0,.5),main=nn)
  hist(dis_w_pred,freq=FALSE,add=T,breaks = 0:50,col=alpha("#648DE5",.5))
  abline(v=dis_act_w,lwd=4)
  
  hist(dis_historical_baseline_e,freq=FALSE,xlim=c(-1,12),breaks = 0:50,col=alpha("black",.5),ylim=c(0,.5),main=nn)
  hist(dis_e_pred,freq=FALSE,add=T,breaks = 0:50,col=alpha("#E98A15",.5))
  abline(v=dis_act_e,lwd=4)
    
}

brier_score_record_df <- as.data.frame(brier_score_record)

brier_score_cat <- brier_score_record_df %>% 
  group_by(V2) %>% 
  summarise(brier_hist = mean(as.numeric(V5)),
            brier_cascade = mean(as.numeric(V6)))
colnames(brier_score_cat)[1] <- c("Outbreak_Size")

brier_score_cat$bss_historical_over_baseline <- 100*(1 - brier_score_cat$brier_hist/0.25)
brier_score_cat$bss_cascade_over_baseline <- 100*(1 - brier_score_cat$brier_cascade/0.25)
brier_score_cat$bss_cascade_over_historical <- 100*(1 - brier_score_cat$brier_cascade/brier_score_cat$brier_hist)

brier_score_cat

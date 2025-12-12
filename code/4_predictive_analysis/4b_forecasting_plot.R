# Code to analyze predictions

library(MCMCvis)
library(dplyr)
library(scoringRules)
library(ggplot2)
library(tidyr)

full_data <- readRDS("data/3_model_in/model_full_data_2024.rds")

blinded_yn <- 0
final_year <- 2024

hist_forecast_predictions_irr <- c()
hist_forecast_predictions_dis <- c()
for (nn in 21:1){
  
  # Read in time machine model fit
  yrs_back <- nn
  print(yrs_back)
  
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
  
  median_pisi_w <- median(pisi_w_pred)
  upper97_5_pisi_w <- quantile(pisi_w_pred,.975)
  lower2_5_pisi_w <- quantile(pisi_w_pred,.025)
  
  median_pisi_e <- median(pisi_e_pred)
  upper97_5_pisi_e <- quantile(pisi_e_pred,.975)
  lower2_5_pisi_e <- quantile(pisi_e_pred,.025)

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
  
  dis_w_score_cascade <- crps_sample(dis_act_w,dis_w_pred)
  dis_e_score_cascade <- crps_sample(dis_act_e,dis_e_pred)
  
  dis_w_score_historical <- crps_sample(dis_act_w,dis_historical_baseline_w)
  dis_e_score_historical <- crps_sample(dis_act_e,dis_historical_baseline_e)
  
  #Get the chance of a large outbreak (>=100 individuals)
  len_chains_w <- length(dis_w_pred)
  real_ob_size_w <- exp(dis_w_pred)-1
  
  zero_inf_w <- sum(real_ob_size_w == 0)/len_chains_w
  zero_to_hundred_inf_w <- sum(real_ob_size_w > 0 & real_ob_size_w <= 100)/len_chains_w
  hundred_to_fivehundred_inf_w <- sum(real_ob_size_w > 100 & real_ob_size_w <= 500)/len_chains_w
  over_fivehundred_inf_w <- sum(real_ob_size_w > 500)/len_chains_w
  
  len_chains_e <- length(dis_e_pred)
  real_ob_size_e <- exp(dis_e_pred)-1
  
  zero_inf_e <- sum(real_ob_size_e == 0)/len_chains_e
  zero_to_hundred_inf_e <- sum(real_ob_size_e > 0 & real_ob_size_e <= 100)/len_chains_e
  hundred_to_fivehundred_inf_e <- sum(real_ob_size_e > 100 & real_ob_size_e <= 500)/len_chains_e
  over_fivehundred_inf_e <- sum(real_ob_size_e > 500)/len_chains_e
  
  hist_forecast_predictions_irr <- rbind(hist_forecast_predictions_irr,
                                     c(median_pisi_w,
                                       upper97_5_pisi_w,
                                       lower2_5_pisi_w,
                                       median_pisi_e,
                                       upper97_5_pisi_e,
                                       lower2_5_pisi_e))
  hist_forecast_predictions_dis <- rbind(hist_forecast_predictions_dis,
                                         c(zero_inf_w,
                                           zero_to_hundred_inf_w,
                                           hundred_to_fivehundred_inf_w,
                                           over_fivehundred_inf_w,
                                           zero_inf_e,
                                           zero_to_hundred_inf_e,
                                           hundred_to_fivehundred_inf_e,
                                           over_fivehundred_inf_e))
}


### West plot ####
west_plt_dt <- hist_forecast_predictions_dis[,1:4] %>% as.data.frame()
colnames(west_plt_dt) <- c("No outbreak","0-100","101-500",">500")
west_plt_dt$years_back <- 21:1
west_plt_dt_long <- west_plt_dt %>% pivot_longer(cols = 1:4)

real_obs_w <- exp(full_data$west_disease_ts[length(full_data$west_disease_ts)-21:1])-1
real_obs_w_f <- rep(NA,length(real_obs_w))
real_obs_w_f[real_obs_w == 0] <- "No outbreak"
real_obs_w_f[real_obs_w > 0 & real_obs_w <= 100] <- "0-100"
real_obs_w_f[real_obs_w > 100 & real_obs_w <= 500] <- "101-500"
real_obs_w_f[real_obs_w > 500] <- ">500"

real_obs_w_f <- real_obs_w_f %>% 
  as.data.frame() %>% 
  mutate(years_back = 21:1) %>%
  arrange(years_back)
colnames(real_obs_w_f)[1] <- "real_ob_size"

west_plt_dt_long$real_ob_size <- NA

for (mm in 1:nrow(west_plt_dt_long)){
  real_ob_cat <- real_obs_w_f$real_ob_size[west_plt_dt_long$years_back[mm]]
  if (real_ob_cat == west_plt_dt_long$name[mm]){
    west_plt_dt_long$real_ob_size[mm] <- "black"
  } else {
    west_plt_dt_long$real_ob_size[mm] <- "none"
  }
}

west_plt_dt_long$transp_test <- .75
west_plt_dt_long$transp_test[west_plt_dt_long$real_ob_size == "black"] <- 1
west_plt_dt_long$year_labels = final_year - west_plt_dt_long$years_back + 1

### East plot ####
east_plt_dt <- hist_forecast_predictions_dis[,5:8] %>% as.data.frame()
colnames(east_plt_dt) <- c("No outbreak","0-100","101-500",">500")
east_plt_dt$years_back <- 21:1
east_plt_dt_long <- east_plt_dt %>% pivot_longer(cols = 1:4)

real_obs_e <- exp(full_data$east_disease_ts[length(full_data$east_disease_ts)-21:1])-1
real_obs_e_f <- rep(NA,length(real_obs_e))
real_obs_e_f[real_obs_e == 0] <- "No outbreak"
real_obs_e_f[real_obs_e > 0 & real_obs_e <= 100] <- "0-100"
real_obs_e_f[real_obs_e > 100 & real_obs_e <= 500] <- "101-500"
real_obs_e_f[real_obs_e > 500] <- ">500"

real_obs_e_f <- real_obs_e_f %>% 
  as.data.frame() %>% 
  mutate(years_back = 21:1) %>%
  arrange(years_back)
colnames(real_obs_e_f)[1] <- "real_ob_size"

east_plt_dt_long$real_ob_size <- NA

for (mm in 1:nrow(east_plt_dt_long)){
  real_ob_cat <- real_obs_e_f$real_ob_size[east_plt_dt_long$years_back[mm]]
  if (real_ob_cat == east_plt_dt_long$name[mm]){
    east_plt_dt_long$real_ob_size[mm] <- "black"
  } else {
    east_plt_dt_long$real_ob_size[mm] <- "none"
  }
}

east_plt_dt_long$transp_test <- .75
east_plt_dt_long$transp_test[east_plt_dt_long$real_ob_size == "black"] <- 1

east_plt_dt_long$year_labels = final_year - east_plt_dt_long$years_back + 1


for(ii in 1:2){
  if (ii == 1){
    data_for_plot <- west_plt_dt_long
  } else {
    data_for_plot <- east_plt_dt_long
  }
  
  p <- ggplot(data=data_for_plot,aes(x=(year_labels),
                                     y=value, 
                                     label = paste(round(value,2)*100,"%",sep=""),
                                     fill = factor(name,levels=rev(c("No outbreak","0-100","101-500",">500"))),
                                     colour=real_ob_size,
                                     alpha=I(transp_test)))+
    geom_hline(yintercept=c(0,.25,.5,.75,1),lty="dashed",alpha=1) +
    geom_bar(stat="identity",linewidth=.75) + 
    geom_text(size = 5, position = position_stack(vjust = 0.5)) +
    scale_x_continuous(breaks=2004:2024) + 
    ylab("") +
    xlab("") +
    theme_bw() +
    scale_fill_manual("legend", 
                      values = rev(c("No outbreak" = "darkolivegreen3",
                                     "0-100" = "gold2",
                                     "101-500" = "darkorange2",
                                     ">500" = "firebrick4"))) +
    scale_color_manual("legend", 
                       values = c("black","transparent")) +
    theme(#legend.position = "none",
          axis.text.y=element_text(size=18),
          axis.text.x=element_text(size=18),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
    ) +
    coord_flip()
    
  plot(p)
}

p + scale_color_manual("legend", 
                       values = c("black","transparent")) +
  theme(#legend.position = "none",
        legend.key.width = unit(1, "cm"),  # Adjust width
        legend.key.height = unit(1, "cm"),
        legend.text = element_text(size = 14),
        axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

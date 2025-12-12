# Code to take full, current dataset and trim to certain years.
# For some metrics (including irruptions), data needs to be read in
library(readr)

tm_total <- 21 # number of years to go back in time
cone_blinded <- 0 # 0 for no, 1 for yes

for (num_years_time_machine in 1:tm_total){
  full_data <- readRDS("data/3_model_in/model_full_data_2024.rds")
  # Data that remains unchanged:
  
  # Ncls - Number of cells
  # Nsp - Number of species
  # N_edges - Number of cell connections
  # node1/node2 - Connections
  # max_wf/max_ef - relative masting potential
  # dt_miss_ind - missing dt
  # fc - Forest cover by cell
  # er_identity - Which region each cell belongs to
  # Ncells_wf - number of western cells
  # Ncells_nf - number of eastern cells
  # N_dis_unobserved - Number of unobserved disease years
  
  # Data that does change
  # Nyrs - Minus number of years in time machine
  full_data$Nyrs <- full_data$Nyrs - num_years_time_machine
  
  # dt_years - Remove last X number of years
  full_data$dt_years <- full_data$dt_years[1:(length(full_data$dt_years)-num_years_time_machine)] 
  
  # mast_years - Remove last X number of years
  full_data$mast_years <- full_data$mast_years[1:(length(full_data$mast_years)-num_years_time_machine)] 
  
  # irr_years - Remove last X number of years
  full_data$irr_years <- full_data$irr_years[1:(length(full_data$irr_years)-num_years_time_machine)] 
  
  # dis_years - Remove last X number of years
  full_data$dis_years <- full_data$dis_years[1:(length(full_data$dis_years)-num_years_time_machine)] 
  
  # N_dt_yrs_obs_else - Minus number of years
  full_data$N_dt_yrs_obs_else <- full_data$N_dt_yrs_obs_else-num_years_time_machine
  
  # dt_obs_in - Remove last X number of years
  full_data$dt_obs_ind <- full_data$dt_obs_ind[1:(length(full_data$dt_obs_ind)-num_years_time_machine)] 
  
  # N_obs_mast - Remove last X number of years
  full_data$N_obs_mast <- full_data$N_obs_mast[1:(length(full_data$N_obs_mast)-num_years_time_machine)] 
  
  # N_miss_mast - Remove last X number of years
  full_data$N_miss_mast <- full_data$N_miss_mast[1:(length(full_data$N_miss_mast)-num_years_time_machine)] 
  
  # yy_obs_mast - Remove last X number of years (rows!)
  full_data$yy_obs_mast <- full_data$yy_obs_mast[1:(nrow(full_data$yy_obs_mast)-num_years_time_machine),]
  
  if(cone_blinded == 1){
    full_data$yy_obs_mast[nrow(full_data$yy_obs_mast),] <- 0.00001
    full_data$N_miss_mast[length(full_data$N_miss_mast)] <- full_data$Ncls
    full_data$N_obs_mast[length(full_data$N_obs_mast)] <- 0
  }
  
  #For each column with at least one time series, recalculate the observed masting values relative to the mean
  for (each_col in 1:ncol(full_data$yy_obs_mast)){
    col_in <- full_data$yy_obs_mast[,each_col]
    if(length(which(col_in>0.00001))>0){
      col_in_known <- col_in[which(col_in>0.00001)]
      col_in_known <- 0.001 + (col_in_known - 0.001)/mean(col_in_known)
      full_data$yy_obs_mast[which(col_in>0.00001),each_col] <- col_in_known
    }
  }
  
  # mean_mast_val,sd_mean_mast - recalculate
  # Get the mean of sd of masting data using bootstrapping
  known_mast_dist <- as.numeric(full_data$yy_obs_mast[which(full_data$yy_obs_mast!= 0.00001)])
  mean_mast_val <- mean(known_mast_dist)
  
  total_mast_vals <- length(as.numeric(full_data$yy_obs_mast))
  
  bs_means <- c()
  for (n in 1:10000){
    bs_sample <- sample(known_mast_dist,total_mast_vals,replace=TRUE)
    bs_means <- c(bs_means,mean(bs_sample))
  }
  sd_mean_mast <- sd(bs_means)
  
  full_data$mean_mast_val <- mean_mast_val
  full_data$sd_mean_mast <- sd_mean_mast
  
  # dt_obs - Remove last X number of years (rows!)
  full_data$dt_obs <- full_data$dt_obs[1:(nrow(full_data$dt_obs)-num_years_time_machine),]
  
  # N_irr_obs_year - Minus X years
  full_data$N_irr_obs_year <- full_data$N_irr_obs_year-num_years_time_machine
  
  # irr_obs_yrs - Remove last X number of years
  full_data$irr_obs_yrs <- full_data$irr_obs_yrs[1:(length(full_data$irr_obs_yrs)-num_years_time_machine)] 
  
  # irr_missing_yrs - Minus one
  full_data$irr_missing_yrs <- full_data$irr_missing_yrs - num_years_time_machine
  
  # irr_ind_e - Read in data for contemporaneously estimated irruptions, add new last year with 0.00001
  contemp_irr_data_e <- read_csv(paste("data/2_formatted/Irruptions/time_machine/irruptions_detrended_east_tm_",num_years_time_machine,".csv",sep=""))
  full_data$irr_ind_e <- contemp_irr_data_e[,2:9]
  full_data$irr_ind_e[nrow(full_data$irr_ind_e)+1,] <- 0.00001 
  
  # irr_ind_w - Remove last X number of years (rows!), add new last year with 0.00001
  contemp_irr_data_w <- read_csv(paste("data/2_formatted/Irruptions/time_machine/irruptions_detrended_west_tm_",num_years_time_machine,".csv",sep=""))
  full_data$irr_ind_w <- contemp_irr_data_w[,2:9]
  full_data$irr_ind_w[nrow(full_data$irr_ind_w)+1,] <- 0.00001 
  
  # ii_obs_mast - Remove last X number of years (rows!)
  full_data$ii_obs_mast <- full_data$ii_obs_mast[1:(nrow(full_data$ii_obs_mast)-num_years_time_machine),]
  
  # ii_miss_mast - Remove last X number of years (rows!)
  full_data$ii_miss_mast <- full_data$ii_miss_mast[1:(nrow(full_data$ii_miss_mast)-num_years_time_machine),]
  
  if(cone_blinded == 1){
    full_data$ii_obs_mast[nrow(full_data$ii_obs_mast),] <- 0
    full_data$ii_miss_mast[nrow(full_data$ii_miss_mast),] <- 1:full_data$Ncls
  }
  
  #outbreak_yn_west - Remove last X number of years, change last year to 0
  full_data$outbreak_yn_west <- full_data$outbreak_yn_west[1:(length(full_data$outbreak_yn_west)-num_years_time_machine)]
  full_data$outbreak_yn_west[length(full_data$outbreak_yn_west)] <- 0
  
  #outbreak_yn_east - Remove last X number of years, change last year to 0
  full_data$outbreak_yn_east <- full_data$outbreak_yn_east[1:(length(full_data$outbreak_yn_east)-num_years_time_machine)]
  full_data$outbreak_yn_east[length(full_data$outbreak_yn_east)] <- 0
  
  # west_disease_ts - Remove last X number of years, change last year to 0
  full_data$west_disease_ts <- full_data$west_disease_ts[1:(length(full_data$west_disease_ts)-num_years_time_machine)]
  full_data$west_disease_ts[length(full_data$west_disease_ts)] <- 0
  
  # east_disease_ts - Remove last X number of years, change last year to 0
  full_data$east_disease_ts <- full_data$east_disease_ts[1:(length(full_data$east_disease_ts)-num_years_time_machine)]
  full_data$east_disease_ts[length(full_data$east_disease_ts)] <- 0
  
  # disease_ts_observed - Remove last X number of years
  full_data$disease_ts_observed <- full_data$disease_ts_observed[1:(length(full_data$disease_ts_observed)-num_years_time_machine)] 
  
  # disease_ts_unobserved - Minus X number of years from last entry
  full_data$disease_ts_unobserved[length(full_data$disease_ts_unobserved)] <- full_data$disease_ts_unobserved[length(full_data$disease_ts_unobserved)] - num_years_time_machine
  
  #dis_obs_vect - Remove last X number of years, change last year to 0
  full_data$dis_obs_vect <- full_data$dis_obs_vect[1:(length(full_data$dis_obs_vect)-num_years_time_machine)]
  full_data$dis_obs_vect[length(full_data$dis_obs_vect)] <- 0
  
  if (cone_blinded ==0){
    saveRDS(full_data,paste("data/3_model_in/time_machine/not_blinded/time_machine_minus_",num_years_time_machine,"_yrs.rds",sep=""))
  } else {
    saveRDS(full_data,paste("data/3_model_in/time_machine/blinded/time_machine_minus_blinded_",num_years_time_machine,"_yrs.rds",sep=""))
  }
  
}

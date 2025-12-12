library(terra)

#For each two year pair, calculate the difference in July temps (year 2 - year 1)
for (first_year in 1980:2022){
  second_year = first_year + 1
  first_july_temp <- rast(paste("data/1_raw/Daymet/tmax_july_",first_year,"_ncss.nc",sep=""))
  second_july_temp <- rast(paste("data/1_raw/Daymet/tmax_july_",second_year,"_ncss.nc",sep=""))
  dt <- second_july_temp - first_july_temp
  
  plot(dt,main=paste(first_year,"-",second_year))
  writeCDF(dt,paste("data/2_formatted/DeltaT/deltaT_rasters/deltaT_",first_year,"_",second_year,".nc",sep=""))
  
}

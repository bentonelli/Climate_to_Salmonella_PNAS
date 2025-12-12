library(daymetr)
library(terra)

for (each_year in 1985:2023) {
  download_daymet_ncss(location = c(80, -172, 0, -45),
                       start = each_year,
                       end = each_year,
                       frequency = "monthly",
                       param = c("tmax"),
                       path = "data/1_raw/Daymet/",
                       silent = FALSE)
  
  year_look <- rast(paste("data/1_raw/Daymet/tmax_monavg_",each_year,"_ncss.nc",sep=""))
  plot(year_look$tmax_7,main=each_year)
  year_save <- year_look$tmax_7
  writeCDF(year_save,paste("data/1_raw/Daymet/tmax_july_",each_year,"_ncss.nc",sep=""))
}


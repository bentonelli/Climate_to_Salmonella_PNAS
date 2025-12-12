# Daymet delta-t rasters to cell format

library(terra)
library(dggridR)
library(dplyr)

#Read in delta t for given year
pdf("delta_t_raster_to_cell.pdf")
for (nn in 1980:2022){
  print(nn)
  start_year <- nn
  end_year <- nn+1
  
  delta_t_raster <- rast(paste("data/2_formatted/DeltaT/deltaT_rasters/deltaT_",start_year,"_",end_year,".nc",sep=""))
  delta_t_raster <- terra::project(delta_t_raster,"+proj=longlat +datum=WGS84")
  plot(delta_t_raster,xlim=c(-180,-30),main=nn)
  
  #Take sample
  reg_sample <- spatSample(delta_t_raster, size=10000000, na.rm=TRUE,
                           method="regular", as.points=FALSE,
                           xy= TRUE,ext=c(-179,-30,10,89))
  
  #Group by cell, average
  hexgrid6 <- dggridR::dgconstruct(res = 6)
  reg_sample$cell <- dgGEO_to_SEQNUM(hexgrid6, reg_sample$x,reg_sample$y)$seqnum
  colnames(reg_sample)[3] <- "dt"
  
  by_cell_mean_dt <- reg_sample %>% 
    group_by(cell) %>% 
    summarise(mean_dt = mean(dt),
              count = n())
  
  #Get cell centers
  #Create a maps to look at where data is coming from:
  by_cell_mean_dt$center_lon   <- dgSEQNUM_to_GEO(hexgrid6,by_cell_mean_dt$cell)$lon_deg
  by_cell_mean_dt$center_lat   <- dgSEQNUM_to_GEO(hexgrid6,by_cell_mean_dt$cell)$lat_deg
  
  saveRDS(by_cell_mean_dt,paste("data/2_formatted/DeltaT/deltaT_cells/deltaT_",start_year,"_",end_year,".rds",sep=""))
  grid  <- dgcellstogrid(hexgrid6,by_cell_mean_dt$cell,return_sf = FALSE)
  grid  <- merge(grid,by_cell_mean_dt,by.x="seqnum",by.y="cell")
  
  grid <- filter(grid,center_lon>-170)
  
  countries <- map_data("world")
  states <- map_data("state")
  countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                              "Belize","Honduras","El Salvador","Nicaragua",
                                              "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                              "Dominican Republic","Puerto Rico","Colombia",
                                              "Venezuela","Brazil","Guyana","Suriname","French Guiana"))
  
  p <- ggplot() + coord_map("mollweide",xlim=c(-170,-60),ylim=c(25,85)) + ggtitle(nn) + 
    geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),axis.text.y = element_blank(),
          axis.title.x=element_blank(),axis.title.y=element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "null"),
          panel.spacing = unit(c(0, 0, 0, 0), "null")) +
    theme(panel.background = element_rect(fill = alpha("lightskyblue2",.5))) +
    geom_path   (data=grid,      aes(x=x, y=y, group=seqnum), alpha=0.4, color="white") +
    geom_polygon(data=grid,      aes(x=x, y=y, group=seqnum,fill=mean_dt), alpha=0.4)    +
    #geom_text(data = cellcenters, aes(x=lon_deg,y=lat_deg,label=cell)) +
    scale_fill_gradient2(low="skyblue3", mid="grey",high="firebrick4",
                         name = "",midpoint=0)
  print(p)
  
}
dev.off()

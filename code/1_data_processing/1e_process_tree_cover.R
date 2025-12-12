# Script to get tree cover by cell

library(terra)
library(dplyr)
library(dggridR)
library(readr)

# Get forest cover
fc_2000 <- rast("data/1_raw/Tree_cover/FC_2000_1km.nc")
fc_2000 <- terra::project(fc_2000,"+proj=longlat +datum=WGS84")
fc_2000_zeroed <- fc_2000

#Take sample
reg_sample <- spatSample(fc_2000_zeroed, 10000000, method="regular", as.points=FALSE,xy= TRUE,ext=c(-179,-30,20,89))


#Group by cell, average
hexgrid6 <- dggridR::dgconstruct(res = 6)
reg_sample$cell <- dgGEO_to_SEQNUM(hexgrid6, reg_sample$x,reg_sample$y)$seqnum

reg_sample$FC_2000_1km[which(is.na(reg_sample$FC_2000_1km))] <- 0
avg_fc_by_cell <- reg_sample %>% group_by(cell) %>% summarise(mean_fc = mean(FC_2000_1km))

table(reg_sample$cell)

saveRDS(avg_fc_by_cell,"data/2_formatted/Tree_cover/FC_by_cell.rds")


### Plot to check out ####

grid  <- dgcellstogrid(hexgrid6,avg_fc_by_cell$cell,return_sf = FALSE)
grid  <- merge(grid,avg_fc_by_cell,by.x="seqnum",by.y="cell")

grid <- filter(grid,mean_fc >0)

countries <- map_data("world")
states <- map_data("state")
countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                            "Belize","Honduras","El Salvador","Nicaragua",
                                            "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                            "Dominican Republic","Puerto Rico","Colombia",
                                            "Venezuela","Brazil","Guyana","Suriname","French Guiana"))

p <- ggplot() + coord_map("mollweide",xlim=c(-180,-30),ylim=c(15,85)) + ggtitle("Tree Cover") +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = alpha("lightskyblue2",.5))) +
  geom_path   (data=grid,      aes(x=x, y=y, group=seqnum), alpha=0.4, color="white") +
  geom_polygon(data=grid,      aes(x=x, y=y, group=seqnum,fill=mean_fc), alpha=0.4)    +
  #geom_text(data = cellcenters, aes(x=lon_deg,y=lat_deg,label=cell)) +
  scale_fill_gradient2(low="grey",high="darkgreen",
                       name = "",lim=c(0,1))
print(p)

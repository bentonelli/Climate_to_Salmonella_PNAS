# Get ecoregion by cell

library(terra)
library(dplyr)
library(dggridR)
library(readr)

#Get ecoregions
us_epa_ecoregions <- vect("data/1_raw/Ecoregions/NA_CEC_Eco_Level1.shp")
us_epa_ecoregions <- terra::project(us_epa_ecoregions,"epsg:4326")
us_epa_ecoregions_trim <- simplifyGeom(us_epa_ecoregions,tolerance=.01)

spat_sample_er <- spatSample(us_epa_ecoregions_trim, 100000, method="regular")
coords_er <- geom(spat_sample_er) %>% as.data.frame()
coords_er$ecoregion <- spat_sample_er$NA_L1NAME

#Group by cell, average
hexgrid6 <- dggridR::dgconstruct(res = 6)
coords_er$cell <- dgGEO_to_SEQNUM(hexgrid6, coords_er$x,coords_er$y)$seqnum

#Get top ecoregion by cell, take out water first
avg_fc_by_cell <- coords_er %>% 
  filter(ecoregion != "WATER") %>%
  #Filter >95% Water cells
  filter(!(cell %in% c(3643,6586,6589,3619,3592,
                       7,84,252,279,307,736,734,3670,
                       3671,344,427,181,121,768,767,719,785))) %>% 
  group_by(cell,ecoregion) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
head(avg_fc_by_cell)

#Save
saveRDS(avg_fc_by_cell,"data/2_formatted/Ecoregions/ER_by_cell.rds")
avg_fc_by_cell <- readRDS("data/2_formatted/Ecoregions/ER_by_cell.rds")

grid  <- dgcellstogrid(hexgrid6,avg_fc_by_cell$cell,return_sf = FALSE)
grid  <- merge(grid,avg_fc_by_cell,by.x="seqnum",by.y="cell")

cellcenters   <- dgSEQNUM_to_GEO(hexgrid6,avg_fc_by_cell$cell)
cellcenters <- as.data.frame(cellcenters)
cellcenters$cell <- avg_fc_by_cell$cell
countries <- map_data("world")
states <- map_data("state")
countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                            "Belize","Honduras","El Salvador","Nicaragua",
                                            "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                            "Dominican Republic","Puerto Rico","Colombia",
                                            "Venezuela","Brazil","Guyana","Suriname","French Guiana"))
#Full map
p <- ggplot() + coord_map("mollweide",xlim=c(-180,-60),ylim=c(15,75)) + 
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = alpha("lightskyblue2",.5))) +
  geom_polygon(data=grid,      aes(x=x, y=y, group=seqnum), alpha=0.4)    +
  geom_path   (data=grid,      aes(x=x, y=y, group=seqnum), alpha=0.4, color="white") +
  geom_text(data = cellcenters, aes(x=lon_deg,y=lat_deg,label=cell)) +
  scale_fill_gradient(low="red", high="blue")
p

### Specific ecoregion ####
er_cells_in <- filter(avg_fc_by_cell,ecoregion == "NORTHERN FORESTS")

grid  <- dgcellstogrid(hexgrid6,er_cells_in$cell,return_sf = FALSE)
grid  <- merge(grid,er_cells_in,by.x="seqnum",by.y="cell")

cellcenters   <- dgSEQNUM_to_GEO(hexgrid6,er_cells_in$cell)
cellcenters <- as.data.frame(cellcenters)
cellcenters$cell <- er_cells_in$cell

p <- ggplot() + coord_map("mollweide",xlim=c(-180,-60),ylim=c(15,75)) + 
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = alpha("lightskyblue2",.5))) +
  geom_polygon(data=grid,      aes(x=x, y=y, group=seqnum), alpha=0.4)    +
  geom_path   (data=grid,      aes(x=x, y=y, group=seqnum), alpha=0.4, color="white") +
  geom_text(data = cellcenters, aes(x=lon_deg,y=lat_deg,label=cell)) +
  scale_fill_gradient(low="red", high="blue")
p

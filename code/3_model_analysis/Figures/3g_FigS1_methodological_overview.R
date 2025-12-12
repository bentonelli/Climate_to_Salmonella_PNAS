#Methodlogical plots
library(dplyr)
library(dggridR)
library(rasterVis)
library(ggplot2)
library(terra)
library(MCMCvis)
library(sf)

### Working backwards, starting with disease data ####
dis_loc_data <- readRDS("data/2_formatted/Disease/disease_locations.rds")
west_salmonellosis_species <- readRDS("data/2_formatted/Disease/west_salmonellosis_ts.rds")
east_salmonellosis_species <- readRDS("data/2_formatted/Disease/east_salmonellosis_ts.rds")
us_epa_ecoregions <- vect("data/1_raw/Ecoregions/NA_CEC_Eco_Level1.shp")

#try to merge ecoregions
names_superregions <- us_epa_ecoregions$NA_L1NAME

west_super_region_list <- c("MARINE WEST COAST FOREST","MEDITERRANEAN CALIFORNIA",
                            "NORTH AMERICAN DESERTS","NORTHWESTERN FORESTED MOUNTAINS",
                            "SOUTHERN SEMIARID HIGHLANDS","TEMPERATE SIERRAS")

east_super_region_list <- c("EASTERN TEMPERATE FORESTS","GREAT PLAINS",
                            "NORTHERN FORESTS","TROPICAL WET FORESTS")

names_l1 <- us_epa_ecoregions$NA_L1NAME
names_superregions[which(names_l1 %in% west_super_region_list)] <- "WEST"
names_superregions[which(names_l1 %in% east_super_region_list)] <- "EAST"
names_superregions[which(!names_l1 %in% c(east_super_region_list,west_super_region_list))] <- "NEITHER"

us_epa_ecoregions$NA_SR <- names_superregions

us_epa_ecoregions <- terra::project(us_epa_ecoregions,"epsg:4326")

ag_map <- aggregate(us_epa_ecoregions,by="NA_SR")
us_epa_ecoregions_trim <- simplifyGeom(ag_map,tolerance=.1) # speed things up
er_trim_crop <- crop(us_epa_ecoregions_trim,c(-170,-50,22,90))
sf_regions <- st_as_sf(er_trim_crop)

#Remove two canadian points, limit to wint_season 2020-2021
dis_loc_data <- dis_loc_data[complete.cases(dis_loc_data$count_lon),] %>% filter(wint_season == 2021)

dis_loc_data$lon_jitter <- dis_loc_data$count_lon + rnorm(length(dis_loc_data$count_lon),0,.1)
dis_loc_data$lat_jitter <- dis_loc_data$count_lat + rnorm(length(dis_loc_data$count_lat),0,.1)

projection_code <- 4326

d_points <- data.frame(long = dis_loc_data$lon_jitter, 
                       lat  = dis_loc_data$lat_jitter) %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  st_transform(crs = 4326)

regional_sf <- sf_regions %>%
  st_transform(crs = projection_code)

countries <- map_data("world")
states <- map_data("state")
countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                            "Belize","Honduras","El Salvador","Nicaragua",
                                            "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                            "Dominican Republic","Puerto Rico","Colombia",
                                            "Venezuela","Brazil","Guyana","Suriname","French Guiana"))
lakes <- map_data("lakes") %>% filter(long < -40)


p1 <- ggplot() + 
  geom_sf(data=regional_sf,aes(fill=NA_SR),alpha=.15,col="grey40",lwd=.25) + 
  scale_fill_manual(values = c("#E98A15","grey40","#648DE5")) +
  geom_sf(data=d_points,aes(col = dis_loc_data$ew_superregion),alpha=.5,size=4) +
  #geom_point(data=d_points,cex=2,alpha=.5,pch=19) + 
  scale_color_manual(values = c("#E98A15","#648DE5","transparent")) +
  theme_classic() + 
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank()) +
  theme(legend.position="none") + 
  xlab("") + ylab("") + 
  coord_sf(xlim = c(-170, -50), ylim = c(25, 75), expand = FALSE)

p1
### Irruption data ####
cbc_data <- readRDS("data/2_formatted/Irruptions/cbc_core_w_cell_filt.rds")

spec_name_list <- c("Pine Siskin","Pine Siskin (green morph)","Pine Siskin (Northern)")

targ_yr_irr <- cbc_data %>% filter(wint_season == 2021) %>% dplyr::select(c(cell,com_name,East_region,West_region))

all_cells <- targ_yr_irr %>% dplyr::select(c(cell,East_region,West_region)) %>% unique()

pisi_cells <- targ_yr_irr$cell[which(targ_yr_irr$com_name %in% spec_name_list)] %>% unique()

all_cells$pisi <- 0

all_cells$pisi[which(all_cells$cell %in% pisi_cells)] <- 1

hexgrid8 <- dggridR::dgconstruct(res = 8) 

grid  <- dgcellstogrid(hexgrid8,all_cells$cell)
grid  <- merge(grid,all_cells,by.x="seqnum",by.y="cell")

cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid8,unique(all_cells$cell))
cellcenters$cell <- unique(all_cells$cell)
cellcenters <- cellcenters %>% as.data.frame()

countries <- map_data("world")
states <- map_data("state")
countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                            "Belize","Honduras","El Salvador","Nicaragua",
                                            "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                            "Dominican Republic","Puerto Rico","Colombia",
                                            "Venezuela","Brazil","Guyana","Suriname","French Guiana"))

col_options <- c("#648DE5","#E98A15")
p2 <- ggplot() + 
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),alpha=1,fill="white",color="darkgrey") +
  #geom_polygon(data=states, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  geom_sf(data=grid,fill=alpha(col_options[all_cells$East_region+1],(grid$pisi+.25)/2))    +
  coord_sf(xlim = c(-170, -50), ylim = c(25, 75), expand = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = alpha("white",.5)))

p2
#geom_polygon(data=grid_e,      aes(x=lon_deg, y=lat, group=group), alpha=0.4,fill=alpha("forestgreen",.4))    +
#geom_path   (data=grid_e,      aes(x=lon_deg, y=lat, group=group), alpha=0.4, color="white") +
#geom_polygon(data=grid_w,      aes(x=lon_deg, y=lat, group=group), alpha=0.4,fill=alpha("dodgerblue4",.4))    +
#geom_path   (data=grid_w,      aes(x=lon_deg, y=lat, group=group), alpha=0.4, color="white") +
#geom_text(data = cellcenters, aes(x=lon_deg,y=lat_deg,label=cell)) 
#geom_point(data=ts_locations_else,aes(x=Longitude,y=Latitude),pch=16,alpha=.8) 

### Masting predictions ####
mdl_in <- readRDS("data/4_model_out/mdl_fit.rds")
mdl_data_in <- readRDS("data/3_model_in/model_full_data_2024.rds")

mast_pred <- MCMCsummary(mdl_in,params="yy_else")
er_id <- mdl_data_in$er_identity

yy_ind <- rep(1:120,42)
yr_2020_2021<- mast_pred[(120*38+1):(120*39),]

hexgrid6 <- dggridR::dgconstruct(res = 6) 

grid2  <- dgcellstogrid(hexgrid6,mdl_data_in$cell)
grid2$region <- er_id[1,]
grid2$mean_pred <- yr_2020_2021$mean
#grid2  <- merge(grid,mdl_data_in,by.x="seqnum",by.y="cell")

cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6,unique(mdl_data_in$cell))
cellcenters$cell <- mdl_data_in$cell
cellcenters <- cellcenters %>% as.data.frame()

countries <- map_data("world")
states <- map_data("state")
countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                            "Belize","Honduras","El Salvador","Nicaragua",
                                            "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                            "Dominican Republic","Puerto Rico","Colombia",
                                            "Venezuela","Brazil","Guyana","Suriname","French Guiana"))

col_options <- c("#E98A15","#648DE5")

size_label <- 100*round(grid2$mean_pred,2)
size_label <- paste(size_label,"%",sep="")
#grid2$mean_pred/max(grid2$mean_pred))
#fill=alpha(col_options[all_cells$East_region+1],(grid$pisi+.25)/2)
p3 <- ggplot() + 
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  geom_polygon(data=lakes, aes(x=long, y=lat, group=group),alpha=1,fill="white",color="darkgrey") +
  #geom_polygon(data=states, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  geom_sf(data=grid2,col=alpha(col_options[grid2$region+1],1),fill=alpha(col_options[grid2$region+1],grid2$mean_pred/max(grid2$mean_pred)))+
  coord_sf(xlim = c(-170, -50), ylim = c(25, 75), expand = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = alpha("white",.5))) +
  #geom_polygon(data=grid_e,      aes(x=lon_deg, y=lat, group=group), alpha=0.4,fill=alpha("forestgreen",.4))    +
  #geom_path   (data=grid_e,      aes(x=lon_deg, y=lat, group=group), alpha=0.4, color="white") +
  #geom_polygon(data=grid_w,      aes(x=lon_deg, y=lat, group=group), alpha=0.4,fill=alpha("dodgerblue4",.4))    +
  #geom_path   (data=grid_w,      aes(x=lon_deg, y=lat, group=group), alpha=0.4, color="white") +
  geom_text(data = cellcenters, aes(x=lon_deg,y=lat_deg,label=size_label),cex=3) 
#geom_point(data=ts_locations_else,aes(x=Longitude,y=Latitude),pch=16,alpha=.8) 
p3



#Need to map these on to cells

### Tree Cover ####
tree_cover_map <- terra::rast("data/1_raw/Tree_cover/FC_2000_1km.nc")
tree_cover_map[is.na(tree_cover_map[])] <- 0 
tree_cover_map <- terra::project(tree_cover_map,"epsg:4326")

gplot(tree_cover_map) + 
  geom_tile(aes(fill = value)) +
  #facet_wrap(~ variable) +
  scale_fill_gradient(low="grey90",high="forestgreen",na.value = "grey90") +
  labs(fill = "") +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="transparent", color="grey20") +
  coord_equal() + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        panel.background = element_rect(fill = "grey90", color = NA), # Set panel background
        plot.background = element_rect(fill = "grey90", color = NA)) +
  coord_sf(xlim = c(-170, -50), ylim = c(20, 85), expand = FALSE)

### Delta-T ####

dt_1 <- terra::rast("data/1_raw/Daymet/tmax_july_2018_ncss.nc")
dt_2 <- terra::rast("data/1_raw/Daymet/tmax_july_2019_ncss.nc")
dt_comb <- dt_2 - dt_1 
dt_comb <- terra::project(dt_comb,"epsg:4326")

#plot(dt_comb)
cols_in <- colorRampPalette(c('skyblue3', 'grey90', 'firebrick3'))(100)

gplot(dt_comb) + 
  geom_tile(aes(fill = value)) +
  #facet_wrap(~ variable) +
  # From: https://stackoverflow.com/questions/75685834/how-do-i-set-a-gradient-color-scale-with-a-fixed-midpoint-with-ggplot
  scale_fill_gradientn(
    colours = cols_in,
    breaks=c(-12,-8,-4,0,4,8,12),
    limits=c(-12, 12),
    na.value = "grey90",
    name = "Delta-T"
  ) + 
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="transparent", color="grey20") +
  coord_equal() + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        panel.background = element_rect(fill = "grey90", color = NA), # Set panel background
        plot.background = element_rect(fill = "grey90", color = NA)) +
  coord_sf(xlim = c(-170, -50), ylim = c(20, 85), expand = FALSE)


### Now time series, starting with disease ####
par(las=1,mar=c(4,8,4,4))

west_disease_ts <- readRDS("data/2_formatted/Disease/west_salmonellosis_ts.rds") %>% filter(wint_season <= 2024 & wint_season >= 1988)
east_disease_ts <- readRDS("data/2_formatted/Disease/east_salmonellosis_ts.rds") %>% filter(wint_season <= 2024 & wint_season >= 1988)

plot(west_disease_ts$wint_season,west_disease_ts$total_cases,type="l",col="#648DE5",xlim=c(1983,2024),lwd=6,ylab="",xlab="",cex.lab=1.6,cex.axis=1.6)
points(east_disease_ts$wint_season,east_disease_ts$total_cases,type="l",col="#E98A15",lwd=6)
title(ylab="", mgp=c(4.25,1,0),cex.lab=1.6)
abline(v=2021,lwd=3,lty=2)

### Irruption time series ####
irr_e <- mdl_data_in$irr_ind_e[,1]
irr_w <- mdl_data_in$irr_ind_w[,1]

plot(mdl_data_in$irr_years[1:(length(mdl_data_in$irr_years)-1)],irr_e[1:(length(irr_e)-1)],type="l",col="#E98A15",lwd=6,ylim=c(-3,3),xlim=c(1983,2024),ylab="",xlab="",cex.lab=1.6,cex.axis=1.6)
points(mdl_data_in$irr_years[1:(length(mdl_data_in$irr_years)-1)],irr_w[1:(length(irr_w)-1)],type="l",col="#648DE5",lwd=6,ylab="",xlab="",cex.lab=1.6,cex.axis=1.6)
abline(v=2021,lwd=3,lty=2)

### Lastly, plot regional masting estimates ####
wrm <- MCMCsummary(mdl_in,params="WRM",probs = c(.055,.5,.945))
erm <- MCMCsummary(mdl_in,params="ERM",probs = c(.055,.5,.945))

plot(mdl_data_in$mast_years,erm$`50%`,type="l",col=alpha("#E98A15",.75),ylab="",xlab="",lwd=6,ylim=c(-1,2.25),xlim=c(1983,2024),cex.lab=1.6,cex.axis=1.6)
polygon(c(mdl_data_in$mast_years,rev(mdl_data_in$mast_years)),c(erm$`5.5%`,rev(erm$`94.5%`)),col=alpha("#E98A15",.2),border=NA)

points(mdl_data_in$mast_years,wrm$`50%`,type="l",col=alpha("#648DE5",.75),lwd=6)
polygon(c(mdl_data_in$mast_years,rev(mdl_data_in$mast_years)),c(wrm$`5.5%`,rev(wrm$`94.5%`)),col=alpha("#648DE5",.2),border=NA)
abline(v=2020,lwd=3,lty=2)

# Derive Irruption Indices for each target bird, and West + East
# E/W defined here by ecoregions

library(dplyr)
library(ggplot2)
library(dggridR)
library(readr)
library(lubridate)
library(zoo)
library(terra)

wint_season_range <- 1983:2024
wint_season_last_known <- 2024

#Read in CBC data
cbc_core <- read_csv("data/1_raw/CBC/US_CA_CBC_data.csv")
cbc_core$date <- mdy_hms(cbc_core$cnt_dt)

# Add year column and winter season column
cbc_core$year <- year(cbc_core$date)
cbc_core$wint_season <- year(cbc_core$date)
cbc_core$wint_season[month(cbc_core$date)==12] <- cbc_core$wint_season[month(cbc_core$date)==12] + 1

#Look at distribution of records by season
table(cbc_core$wint_season)

hexgrid8 <- dggridR::dgconstruct(res = 8) 

#Get ecoregions info, simplify
us_epa_ecoregions <- vect("data/1_raw/Ecoregions/NA_CEC_Eco_Level1.shp")
us_epa_ecoregions <- terra::project(us_epa_ecoregions,"epsg:4326")
us_epa_ecoregions_trim <- simplifyGeom(us_epa_ecoregions,tolerance=.01)


#### Assign ecoregions to each cell based on regular sampling ####
# This only needs to be done once

# spat_sample_er <- spatSample(us_epa_ecoregions_trim, 1000000, method="regular")
#  
# coords_er <- geom(spat_sample_er) %>% as.data.frame()
# coords_er$ecoregion <- spat_sample_er$NA_L1NAME
# coords_er$cell <- dgGEO_to_SEQNUM(hexgrid8, coords_er$x,coords_er$y)$seqnum
# avg_fc_by_cell <- coords_er %>%
#   filter(ecoregion != "WATER") %>%
#   group_by(cell) %>%
#   count(ecoregion) %>%
#   slice(which.max(n))
# 
# head(avg_fc_by_cell)
# 
# saveRDS(avg_fc_by_cell,"data/2_formatted/Ecoregions/ER_by_cell8.rds")

### Load ecoregion information####
avg_fc_by_cell <- readRDS("data/2_formatted/Ecoregions/ER_by_cell8.rds")

cbc_core$cell <- dgGEO_to_SEQNUM(hexgrid8, cbc_core$longitude,cbc_core$latitude)$seqnum
cbc_core_w_cell <- merge(cbc_core,avg_fc_by_cell,by="cell",all.x=TRUE)

grid  <- dgcellstogrid(hexgrid8,cbc_core_w_cell$cell)
grid  <- merge(grid,cbc_core_w_cell,by.x="seqnum",by.y="cell")

#Eastern Range
east_er <- c("HUDSON PLAIN",
             "NORTHERN FORESTS",
             "EASTERN TEMPERATE FORESTS",
             "GREAT PLAINS",
             "TROPICAL WET FORESTS")

#Western Range
west_er <- c("NORTHWESTERN FORESTED MOUNTAINS",
             "MARINE WEST COAST FOREST",
             "NORTH AMERICAN DESERTS",
             "MEDITERRANEAN CALIFORNIA",
             "SOUTHERN SEMI-ARID HIGHLANDS")

#Filter data to areas within study area. Then add E/W designation
cbc_core_w_cell_filt <- cbc_core_w_cell %>% 
  filter(ecoregion %in% c(east_er,west_er))

cbc_core_w_cell_filt$East_region <- 0
cbc_core_w_cell_filt$East_region[which(cbc_core_w_cell_filt$ecoregion %in% east_er)] <- 1

cbc_core_w_cell_filt$West_region <- 0
cbc_core_w_cell_filt$West_region[which(cbc_core_w_cell_filt$ecoregion %in% west_er)] <- 1

#Save all species names, including subspecies reported in the CBC database
spec_name_bbl <- list(c("Pine Siskin","Pine Siskin (green morph)","Pine Siskin (Northern)"),
                      c("Red-breasted Nuthatch"),
                      c("Evening Grosbeak"),
                      c("Pine Grosbeak"),
                      c("Purple Finch","Purple Finch (Eastern)","Purple Finch (Western)"),
                      c("redpoll sp.","Common Redpoll","Common Redpoll (flammea)",
                        "Common x Hoary Redpoll (hybrid)","Common/Hoary Redpoll","Hoary Redpoll"),
                      c("Red Crossbill","Red Crossbill (Appalachian or type 1)",
                        "Red Crossbill (Douglas-fir or type 4)","Red Crossbill (Lodgepole Pine or type 5)",
                        "Red Crossbill (Ponderosa Pine or type 2)","Red Crossbill (Sierra Madre or type 6)",
                        "Red Crossbill (Sitka Spruce or type 10)","Red Crossbill (Western Hemlock or type 3)"),
                      c("White-winged Crossbill","White-winged Crossbill (leucoptera)"))


cbc_core_w_cell_filt <- filter(cbc_core_w_cell_filt,wint_season <= wint_season_last_known)

#Save intermediate data
#saveRDS(cbc_core_w_cell_filt,"data/2_formatted/Irruptions/cbc_core_w_cell_filt.rds")
cbc_core <- cbc_core_w_cell <- NULL

#### Eastern region irruption data ####

all_spec_frame <- data.frame(wint_season=wint_season_range)
total_records <- c()
for (spec_name_list in spec_name_bbl){
  #spec_name_list <- spec_name_bbl[[1]]
  
  #Get how many individuals counted of each species
  tot_spec_count <- cbc_core_w_cell_filt %>% 
    filter(com_name %in% spec_name_list) %>% 
    summarise(total_records = sum(how_many))
  total_records <- c(total_records,tot_spec_count)
  
  spec_name <- spec_name_list[1]
  print(spec_name_list[1])
  
  # Get total cells counted
  spec_irr_total_cells <- cbc_core_w_cell_filt %>% 
    filter(East_region == 1) %>%
    group_by(wint_season) %>%
    summarise(total_cell=length(unique(cell)))
  
  # Get number of cells with target bird reported for each year
  spec_irr_count <- cbc_core_w_cell_filt %>% 
    filter(East_region == 1) %>%
    group_by(wint_season) %>%
    filter(com_name %in% spec_name_list) %>% 
    summarise(unq_cells_seen = length(unique(cell)))
  
  #Calculate percent of areas that have target bird reported
  spec_irr_ind <- merge(spec_irr_total_cells,spec_irr_count,by="wint_season")
  spec_irr_ind$perc_cells <- spec_irr_ind$unq_cells_seen/spec_irr_ind$total_cell
  
  if (nrow(spec_irr_ind)>10){
    
    #Detrend percent of cells observed time series
    dt_perc <- smooth.spline(spec_irr_ind$wint_season, spec_irr_ind$perc_cells,df=8)
    spec_irr_ind$perc_cell_trend <- dt_perc$y
    spec_irr_ind$perc_cell_dt <- spec_irr_ind$perc_cells - spec_irr_ind$perc_cell_trend 
    spec_irr_ind <- filter(spec_irr_ind,wint_season %in% wint_season_range)
    spec_irr_ind$norm_irr_ind <- (spec_irr_ind$perc_cell_dt)/sd(spec_irr_ind$perc_cell_dt)

    #Add to df
    spec_to_add <- dplyr::select(spec_irr_ind,c(wint_season,norm_irr_ind))
    colnames(spec_to_add)[2] <- spec_name
    all_spec_frame <-merge(all_spec_frame,spec_to_add,by="wint_season",all=TRUE)
  }
}

all_spec_frame <- filter(all_spec_frame,wint_season %in% wint_season_range)

#Look at species correlations, just for fun!
corr_all_spec <- cor(all_spec_frame[,2:ncol(all_spec_frame)],use ="complete.obs")
corrplot::corrplot(corr_all_spec,order = "AOE")

write_csv(all_spec_frame,"data/2_formatted/Irruptions/irruptions_detrended_east.csv")

#### Western region irruption data ####

total_records <- c()
all_spec_frame <- data.frame(wint_season=wint_season_range)
for (spec_name_list in spec_name_bbl){
  #spec_name_list <- spec_name_bbl[[1]]
  
  tot_spec_count <- cbc_core_w_cell_filt %>% 
    filter(com_name %in% spec_name_list) %>% 
    summarise(total_records = sum(how_many))
  total_records <- c(total_records,tot_spec_count)
  
  spec_name <- spec_name_list[1]
  print(spec_name_list[1])
  
  # Get total cells counted
  
  spec_irr_total_cells <- cbc_core_w_cell_filt %>% 
    filter(East_region == 0) %>%
    group_by(wint_season) %>%
    summarise(total_cell=length(unique(cell)))
  
  spec_irr_count <- cbc_core_w_cell_filt %>% 
    filter(East_region == 0) %>%
    group_by(wint_season) %>%
    filter(com_name %in% spec_name_list) %>% 
    summarise(unq_cells_seen = length(unique(cell)))
  
  spec_irr_ind <- merge(spec_irr_total_cells,spec_irr_count,by="wint_season")
  spec_irr_ind$perc_cells <- spec_irr_ind$unq_cells_seen/spec_irr_ind$total_cell
  
  if (nrow(spec_irr_ind)>10){
    
    #Detrend percent of cells observed time series
    dt_perc <- smooth.spline(spec_irr_ind$wint_season, spec_irr_ind$perc_cells,df=8)
    spec_irr_ind$perc_cell_trend <- dt_perc$y
    spec_irr_ind$perc_cell_dt <- spec_irr_ind$perc_cells - spec_irr_ind$perc_cell_trend 
    spec_irr_ind <- filter(spec_irr_ind,wint_season %in% wint_season_range)
    spec_irr_ind$norm_irr_ind <- (spec_irr_ind$perc_cell_dt)/sd(spec_irr_ind$perc_cell_dt)
    
    #Add to df
    spec_to_add <- dplyr::select(spec_irr_ind,c(wint_season,norm_irr_ind))
    colnames(spec_to_add)[2] <- spec_name
    all_spec_frame <-merge(all_spec_frame,spec_to_add,by="wint_season",all=TRUE) 
  }
}

all_spec_frame <- filter(all_spec_frame,wint_season %in% wint_season_range)

#total records for each species, winters 1983-2024
print(total_records)

corr_all_spec <- cor(all_spec_frame[,2:ncol(all_spec_frame)],use ="complete.obs")
corrplot::corrplot(corr_all_spec,order = "AOE")

write_csv(all_spec_frame,"data/2_formatted/Irruptions/irruptions_detrended_west.csv")


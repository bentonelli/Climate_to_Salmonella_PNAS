
library(dplyr)
library(readr)
library(ggplot2)
library(MCMCvis)
library(terra)
library(sf)
library(dggridR)
library(reshape2)
library(maps)
library(mapproj)

#Set the full range of years of KNOWN data for the analysis
#End years for delta-t
delta_t_year_targets <- 1981:2022

#Mast years
#masting_year_targets <- 1982:2022
masting_year_targets <- delta_t_year_targets+1

#Irruption winter end years
#irruption_year_targets <- 1983:2023
irruption_year_targets <- masting_year_targets+1

#Disease winter end years
disease_year_targets <- 1988:2024

#Number of prediction years to add on to the end of the model
pred_yr_num <- 1

#WHISPers disease dates including unknowns
whispers_all_yrs <- 1983:max(irruption_year_targets+pred_yr_num)


### Masting Data ####

#Read in MASTREE data
mastree <- read_csv("data/1_raw/MASTREE/MASTREEplus_2024-06-26_V2.csv") %>% 
  filter(Year %in% masting_year_targets)

#Get ecoregions, simplify geometry (basically, round the edges of the ecoregions)
us_epa_ecoregions <- vect("data/1_raw/Ecoregions/NA_CEC_Eco_Level1.shp")
us_epa_ecoregions <- terra::project(us_epa_ecoregions,"epsg:4326")

#This speeds up processing
us_epa_ecoregions_trim <- simplifyGeom(us_epa_ecoregions,tolerance=.01)

# Filter to specific ecoregions
er_targets <- us_epa_ecoregions_trim[us_epa_ecoregions_trim$NA_L1NAME %in% 
                                       c("NORTHERN FORESTS",
                                         "NORTHWESTERN FORESTED MOUNTAINS",
                                         "MARINE WEST COAST FOREST")]

# Define two-year masting genera we want data for.
genus_patterns <- c("Abies","Picea","Tsuga")

# Filter to those groups, and only cone counts
mt_na <- mastree %>% 
  filter(grepl(paste(genus_patterns, collapse="|"), Species) & VarType == "C")

#Get locations of each times series, filter dataset to those that fall within the above ecoregions
tree_locs <- cbind(mt_na$Longitude,mt_na$Latitude) %>% 
  as.data.frame()
inter_sects <- st_intersects(st_as_sf(tree_locs, 
                                      coords = c(1:2),
                                      crs="epsg:4326"),
                             st_as_sf(er_targets))
mt_er <- mt_na[!is.na(as.numeric(inter_sects)),]

#Get unique locations
locations_unique_else <- unique(dplyr::select(mt_er,c("Longitude","Latitude")))

#Plot all unique masting data locations
countries <- map_data("world")
states <- map_data("state")
countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                            "Belize","Honduras","El Salvador","Nicaragua",
                                            "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                            "Dominican Republic","Puerto Rico","Colombia",
                                            "Venezuela","Brazil","Guyana","Suriname","French Guiana"))

p2 <- ggplot() + coord_map("mollweide",xlim=c(-140,-60),ylim=c(25,75)) +
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = ggplot2::alpha("white",.5))) +
  geom_point(data=locations_unique_else,aes(x=Longitude, y=Latitude), alpha=0.7,size=3,col="orange",pch=17)
#print(p2)

#Get the percent of the maximum value for each time series
#mt_er$perc_max <- mt_er$Value/mt_er$Max_value

# For each data point, get a cell ID
hexgrid6 <- dggridR::dgconstruct(res = 6)
mt_er$cell <- dgGEO_to_SEQNUM(hexgrid6, mt_er$Longitude,mt_er$Latitude)$seqnum

#Get cells in three target ecoregions
#Filter to where at least 25% of the cell is comprised of one of the target ecoregions
er_cells <- readRDS("data/2_formatted/Ecoregions/ER_by_cell.rds")
inc_er_cells <- er_cells %>% 
  filter(freq >= 0.25) %>% 
  filter(ecoregion %in% c("NORTHERN FORESTS",
                          "NORTHWESTERN FORESTED MOUNTAINS",
                          "MARINE WEST COAST FOREST")) %>% 
  dplyr::select(cell) %>% 
  unique() %>% 
  unlist() %>% 
  as.numeric()

#Of the cells included, pick most-represented ecoregion for indexing later on.
er_cells_ind <- er_cells %>% 
  filter(freq >= 0.25) %>% 
  filter(ecoregion %in% c("NORTHERN FORESTS",
                          "NORTHWESTERN FORESTED MOUNTAINS",
                          "MARINE WEST COAST FOREST")) %>% 
  dplyr::select(cell,freq)

mt_na_filt <- mt_er 

table(mt_na_filt$cell)

#Save a list of cells with data
cells_else <- unique(mt_na_filt$cell)

### Merge time series data from same cells to an average percent of max value ####

# This code deals with the possibility of a multiple time series coming from a single
# cell in a given year. This code calculates the average % of max metrics from all time series

#Keep only useful columns
masting_dt_else <- mt_na_filt %>% 
  dplyr::select(c(Year,cell,Alpha_Number,Site_number,Segment,Variable_number,Species,Longitude,Latitude,Value,Max_value))

#Get unique locations
unique_loc_else <- mt_na_filt %>% dplyr::select(c(Alpha_Number,Site_number,Variable_number,Species,Longitude,Latitude)) %>% unique()

unique_loc_else$unq_code <- 1:nrow(unique_loc_else)

mt_er_else <- merge(masting_dt_else,unique_loc_else,by=c("Alpha_Number","Site_number","Variable_number","Species","Longitude","Latitude"))

mt_er_else$delta_t <- NA

code_count_else <- mt_er_else %>% group_by(unq_code) %>% summarise(code_count = n())
mt_er_else <- merge(mt_er_else,code_count_else)

#mt_er_else$perc_max_trans = mt_er_else$perc_max + 0.001 # Add small number to avoid zeros, changed from 0.001 in an attempt to model as log(yy)

# Get average masting by unique code (unique time series)
mean_mast_vals <- mt_er_else %>% 
  group_by(unq_code) %>% 
  summarise(mean_mast_val = mean(Value))

mt_er_else <- merge(mt_er_else,mean_mast_vals,by="unq_code")

mt_er_else$perc_mean <- mt_er_else$Value/mt_er_else$mean_mast_val + .001 ## Add small number to avoid zeros

print("Unique time series")
length(unique(mt_er_else$unq_code))

print("Species:")
unique(mt_er_else$Species)

print("Total data points:")
nrow(mt_er_else)
#Average masting when more than one measurement per cell/year
by_cl_yr_avg_else <- mt_er_else %>% 
  group_by(cell, Year) %>% 
  summarise(perc_mean_avg = mean(perc_mean))

print("Compiled data points:")
nrow(by_cl_yr_avg_else)

print("Total cells represented:")
length(unique(by_cl_yr_avg_else$cell))

#### Create data for model, fill with NAs ####
# To estimate masting across the whole ecoregion, add cells from ecoregion not represented in the dataset
add_cells <- unique(c(er_cells_ind$cell,unique(by_cl_yr_avg_else$cell)))
add_cells <- add_cells[order(add_cells)]

#Don't add cells that will already will be estimated
add_cells_else <- add_cells[!c(add_cells %in% c(by_cl_yr_avg_else$cell))]

add_rows_else <- data.frame(cell = add_cells_else,
                            year = rep(2018,length(add_cells_else)),
                            perc_mean_avg = rep(NA,length(add_cells_else)))


dt_for_model_else <- by_cl_yr_avg_else
colnames(dt_for_model_else)[2] <- "year"

dt_for_model_else <- dt_for_model_else %>% arrange(cell,year)
dt_for_model_else <- rbind(dt_for_model_else,add_rows_else) %>% arrange(cell,year)

cells <- unique(dt_for_model_else$cell)

# Set up 0/1 matrix so that Western Forests are 1st column and Northern Forests are 2nd column
identity_cells <- c()
identity_matrix <- c()
for (each_inc_cell in cells){
  cl_to_add <- filter(er_cells,cell == each_inc_cell) %>% 
    filter(ecoregion %in% c("NORTHERN FORESTS",
                            "NORTHWESTERN FORESTED MOUNTAINS",
                            "MARINE WEST COAST FOREST"))
  new_row <- cl_to_add[which(cl_to_add$freq == max(cl_to_add$freq)),]
  identity_cells <- rbind(identity_cells,new_row)
  if (new_row$ecoregion[1] == "NORTHERN FORESTS"){
    identity_matrix <- rbind(identity_matrix,c(new_row$cell[1],0,1))
  } else if (new_row$ecoregion[1] %in% c("NORTHWESTERN FORESTED MOUNTAINS","MARINE WEST COAST FOREST")){
    identity_matrix <- rbind(identity_matrix,c(new_row$cell[1],1,0))
  } else {
    identity_matrix <- rbind(identity_matrix,c(new_row$cell[1],0,0))
  }
}

#Create a maps to look at where data is coming from:
lat_centers   <- dgSEQNUM_to_GEO(hexgrid6,cells)$lat_deg
lon_centers   <- dgSEQNUM_to_GEO(hexgrid6,cells)$lon_deg

#Get unique locations for timeseries
ts_locations_else <- mt_er_else %>% dplyr::select(c(Longitude,Latitude)) %>% unique()

grid_e  <- dgcellstogrid(hexgrid6,cells[which(identity_matrix[,3]==1)])
grid_e  <- merge(grid_e,dt_for_model_else,by.x="seqnum",by.y="cell")

grid_w  <- dgcellstogrid(hexgrid6,cells[which(identity_matrix[,2]==1)])
grid_w  <- merge(grid_w,dt_for_model_else,by.x="seqnum",by.y="cell")

cellcenters_e   <- dgSEQNUM_to_GEO(hexgrid6,grid_e$seqnum)
cellcenters_e <- as.data.frame(cellcenters_e)

cellcenters_w   <- dgSEQNUM_to_GEO(hexgrid6,grid_w$seqnum)
cellcenters_w <- as.data.frame(cellcenters_w)

grid_e <- cbind(grid_e,cellcenters_e)
grid_w <- cbind(grid_w,cellcenters_w)

countries <- map_data("world")
states <- map_data("state")
countries <- filter(countries,region %in% c("Canada","USA","Mexico","Guatemala","Bahamas",
                                            "Belize","Honduras","El Salvador","Nicaragua",
                                            "Costa Rica","Panama","Cuba","Jamaica","Haiti",
                                            "Dominican Republic","Puerto Rico","Colombia",
                                            "Venezuela","Brazil","Guyana","Suriname","French Guiana"))

p <- ggplot() + #coord_map("mollweide",xlim=c(-150,-60),ylim=c(25,70)) + 
  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  geom_polygon(data=states, aes(x=long, y=lat, group=group), fill="grey90", color="darkgrey") +
  geom_sf(data=grid_e,fill=alpha("tan4",.4))    +
  geom_sf(data=grid_w,fill=alpha("dodgerblue3",.4))    +
  coord_sf(xlim = c(-170, -50), ylim = c(25, 75), expand = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.spacing = unit(c(0, 0, 0, 0), "null")) +
  theme(panel.background = element_rect(fill = alpha("lightskyblue2",.5))) +
  #geom_polygon(data=grid_e,      aes(x=lon_deg, y=lat, group=group), alpha=0.4,fill=alpha("forestgreen",.4))    +
  #geom_path   (data=grid_e,      aes(x=lon_deg, y=lat, group=group), alpha=0.4, color="white") +
  #geom_polygon(data=grid_w,      aes(x=lon_deg, y=lat, group=group), alpha=0.4,fill=alpha("dodgerblue4",.4))    +
  #geom_path   (data=grid_w,      aes(x=lon_deg, y=lat, group=group), alpha=0.4, color="white") +
  #geom_text(data = cellcenters, aes(x=lon_deg,y=lat_deg,label=cells)) +
  geom_point(data=ts_locations_else,aes(x=Longitude,y=Latitude),pch=16,alpha=.8) 
p

#table(dt_for_model_else$cell,dt_for_model_else$year)

obs_data_func <- function(dt_for_model,yrs_target){
  #nyr <- length(unique(c(dt_for_model$year,dt_for_model$year)))
  nyr <- length(yrs_target)
  ncell <- length(unique(dt_for_model$cell))
  #years <- unique(dt_for_model$year) 
  years <- yrs_target[order(yrs_target)]
  cells <- unique(dt_for_model$cell)
  
  #create and fill sds, obs, data for PPC, and temp
  #sigma_y_in <- matrix(nrow = nyr, ncol = ncell)
  y_obs_in <- matrix(nrow = nyr, ncol = ncell)
  y_PPC <- rep(NA, nyr * ncell)
  
  #number of observation and NAs for each year
  len_y_obs_in <- rep(NA, nyr)
  len_y_mis_in <- rep(NA, nyr)
  
  #indices for observed and missing
  ii_obs_in <- matrix(NA, nrow = nyr, ncol = ncell)
  ii_mis_in <- matrix(NA, nrow = nyr, ncol = ncell)
  
  counter <- 1
  for (i in 1:nyr){
    temp_yr <- dplyr::filter(dt_for_model, year == years[i])
    
    cls_missing <- cells[!(cells %in% temp_yr$cell)]
    na_data_to_add <- data.frame(cell = cls_missing,
                                 year = rep(years[i],length(cls_missing)),
                                 perc_mean_avg = rep(NA,length(cls_missing))) # Following example
    temp_yr <- rbind(temp_yr,na_data_to_add) %>% arrange(cell)
    
    
    #which are not NA
    no_na <- temp_yr$perc_mean_avg[which(!is.na(temp_yr$perc_mean_avg))]
    
    #pad end with NAs
    if (length(no_na) < ncell){
      num_na <- ncell - length(no_na)
      
      #add NAs to end
      t_y_obs_in <- c(no_na, rep(NA, num_na))
      t_obs_in <- c(which(!is.na(temp_yr$perc_mean_avg)), rep(NA, num_na)) 
      t_mis_in <- c(which(is.na(temp_yr$perc_mean_avg)), rep(NA, length(no_na)))
      
      #fill objects
      ii_obs_in[i,] <- t_obs_in
      ii_mis_in[i,] <- t_mis_in
      y_obs_in[i,] <- t_y_obs_in
    } else {
      #no NAs to end (no missing values)
      y_obs_in[i,] <- no_na
      ii_obs_in[i,] <- 1:length(no_na)
      ii_mis_in[i,] <- rep(NA, length(no_na))
    }
    
    #length of data/miss for each year
    len_y_obs_in[i] <- length(no_na)
    len_y_mis_in[i] <- ncell - length(no_na)
  }
  
  #Convert masting data to long format
  yy_data_long <- matrix(NA, nrow = nyr, ncol = ncell)
  for (nn in 1:nrow(dt_for_model)){
    yr_ind <- which(years == as.numeric(dt_for_model[nn,2]))
    cl_ind <- which(cells == as.numeric(dt_for_model[nn,1]))
    yy_data_long[yr_ind,cl_ind] <- as.numeric(dt_for_model[nn,3])
  }
  
  #fill 0 where NA in y_obs - Stan does not like NA and zeros are not being used to estimate any param (y_obs is used to fill y)
  yy_data_long[which(is.na(yy_data_long))] <- 0.00001
  
  ii_obs_in[which(is.na(ii_obs_in))] <- 0
  ii_mis_in[which(is.na(ii_mis_in))] <- 0
  
  list_out <- list(yy_data_long,ii_obs_in,ii_mis_in,len_y_obs_in,len_y_mis_in,nyr,ncell,years,cells)
  names(list_out) <- c("yy_data_long","ii_obs_in","ii_mis_in","len_y_obs_in",'len_y_mis_in',"nyr","ncell","years","cells")
  return(list_out)
}

mast_years_targ <- min(masting_year_targets):(max(masting_year_targets)+pred_yr_num)
obs_out_else <- obs_data_func(dt_for_model_else,mast_years_targ)

### Get the mean and sd of masting data using bootstrapping ####
known_mast_dist <- as.numeric(obs_out_else$yy_data_long[which(obs_out_else$yy_data_long != 0.00001)])
mean_mast_val <- mean(known_mast_dist)
total_mast_vals <- length(as.numeric(obs_out_else$yy_data_long))

bs_means <- c()
for (n in 1:10000){
  bs_sample <- sample(known_mast_dist,total_mast_vals,replace=TRUE)
  bs_means <- c(bs_means,mean(bs_sample))
}
sd_mean_mast <- sd(bs_means)

#### ####
# create adjacency matrix -------------------------------------------------

cells <- unique(dt_for_model_else$cell)

#make hexgrid
hexgrid6 <- dggridR::dgconstruct(res = 6)

#get hexgrid cell centers
cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6,cells)

#add lat col to df
#f_out$lat <- cellcenters$lat_deg

ncell <- length(cells)
#create adjacency matrix - 1 if adjacent to cell, 0 if not
adjacency_matrix <- matrix(data = NA, nrow = ncell, ncol = ncell)

for (i in 1:ncell){
  #i <- 1
  for (j in i:ncell){
    #j <- 69
    dists <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                              c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
    adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
  }
}

#indices for 1s
ninds <- which(adjacency_matrix == 1, arr.ind = TRUE)


#if a cell doesn't border any other cells, drop it and redefine objects
DROP <- FALSE

s_cols <- apply(adjacency_matrix, 2, function(x) sum(x, na.rm = TRUE))
s_rows <- apply(adjacency_matrix, 1, function(x) sum(x, na.rm = TRUE))
to.rm.ind <- which((s_cols + s_rows) == 0)

if (length(to.rm.ind) > 0){
  DROP <- cells[to.rm.ind]
  
  cells_r <- cells[-to.rm.ind]
  ncell <- length(cells_r)
  #f_out <- dplyr::filter(f_out, cell %in% cells_r)
  cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, cells_r)
  
  #create adjacency matrix - 1 if adjacent to cell, 0 if not
  adjacency_matrix <- matrix(data = NA, nrow = ncell, ncol = ncell)
  
  for (i in 1:ncell){
    for (j in i:ncell){
      dists <- geosphere::distm(c(cellcenters$lon_deg[i], cellcenters$lat_deg[i]),
                                c(cellcenters$lon_deg[j], cellcenters$lat_deg[j]))
      adjacency_matrix[i,j] <- as.numeric((dists/1000) > 0 & (dists/1000) < 311)
    }
  }
  
  #indices for 1s
  ninds <- which(adjacency_matrix == 1, arr.ind = TRUE)
}

#### Get Delta-T data ####
# First argument is how long before masting year cone development starts
# 2 for non-pines. Second argument is know DT years. 
dt_calculator <- function(start_year_fact,years,cells){
  dt_long = matrix(NA, nrow = length(years), ncol = length(cells))
  counter = 1
  for (each_year in years){
    start_year <- each_year - start_year_fact
    end_year <- start_year + 1
    
    f_name_yr <- paste("data/2_formatted/DeltaT/deltaT_cells/deltaT_",start_year,"_",end_year,".rds",sep="")
    if(file.exists(f_name_yr) & end_year <= max(years)){
      dt_by_year <- readRDS(f_name_yr) %>% 
        filter(cell %in% cells) %>%
        dplyr::select(c(cell,mean_dt,count))
      dt_long[counter,] <- dt_by_year$mean_dt 
    }
    counter <- counter + 1
  }
  return(dt_long)
}

dt_long_else <- dt_calculator(2,mast_years_targ,cells)

yrs_dt_else_obs <- as.numeric(complete.cases(dt_long_else))

#### Tree cover data ####
fc_by_cell <- readRDS("data/2_formatted/Tree_cover/FC_by_cell.rds")

fc_by_cell_filt <- fc_by_cell[which(fc_by_cell$cell %in% cells),]

# Read in irruption data
### Irruptions ####
irr_dt_e <- read_csv("data/2_formatted/Irruptions/irruptions_detrended_east.csv") %>% 
  filter(wint_season %in% irruption_year_targets) %>%
  dplyr::select(c(wint_season,`Pine Siskin`,`Evening Grosbeak`,`redpoll sp.`, # First four are likely disease culprits
                  `Purple Finch`,`Pine Grosbeak`,`Red Crossbill`,
                  `White-winged Crossbill`,`Red-breasted Nuthatch`))

irr_dt_w <- read_csv("data/2_formatted/Irruptions/irruptions_detrended_west.csv") %>% 
  filter(wint_season %in% irruption_year_targets) %>%
  dplyr::select(c(wint_season,`Pine Siskin`,`Evening Grosbeak`,`redpoll sp.`, # First four are likely disease culprits
                  `Purple Finch`,`Pine Grosbeak`,`Red Crossbill`,
                  `White-winged Crossbill`,`Red-breasted Nuthatch`))

#Add NA years on to the end
irr_yr_add_on_targets <- min(irruption_year_targets):max(irruption_year_targets+pred_yr_num)
add_irr_yrs <- irr_yr_add_on_targets[!c(irr_yr_add_on_targets %in% irr_dt_e$wint_season)]

for (each_yr in add_irr_yrs){
  irr_dt_na_add <- c(each_yr,rep(NA,ncol(irr_dt_e)-1))
  
  irr_dt_e <- rbind(irr_dt_e,irr_dt_na_add)
  irr_dt_w <- rbind(irr_dt_w,irr_dt_na_add)
}
irr_dt_e_trim <- irr_dt_e
irr_dt_w_trim <- irr_dt_w

irr_obs_yrs <- which(!(irr_dt_e_trim$wint_season %in% add_irr_yrs))
irr_missing_yrs <- which((irr_dt_w_trim$wint_season %in% add_irr_yrs))

N_irr_obs_year <- length(irr_obs_yrs)

#Switch to matrix, replace NAs with .0001
irr_dt_trim_e_ind <- as.matrix(irr_dt_e_trim[,2:ncol(irr_dt_e_trim)])
irr_east_nas <- as.numeric(complete.cases(irr_dt_trim_e_ind))
#Set 0s to very small number, as stan can't handle 0s
irr_dt_trim_e_ind[is.na(irr_dt_trim_e_ind)] <- .00001

irr_dt_trim_w_ind <- as.matrix(irr_dt_w_trim[,2:ncol(irr_dt_w_trim)])
irr_west_nas <- as.numeric(complete.cases(irr_dt_trim_w_ind))
#Set 0s to very small number, as stan can't handle 0s
irr_dt_trim_w_ind[is.na(irr_dt_trim_w_ind)] <- .00001

#Set 0s to very small number, as stan can't handle 0s
dt_long_else[which(is.na(dt_long_else))] <- .00001

#Get disease data, log+1 and standardize

### Disease data ####

west_disease_ts <- readRDS("data/2_formatted/Disease/west_salmonellosis_ts.rds") %>% filter(wint_season %in% disease_year_targets)

dis_yr_add_on_targets <- whispers_all_yrs[!c(whispers_all_yrs %in% west_disease_ts$wint_season)] 

#miss_ind_disease <- which(!irr_yr_target %in% west_disease_ts$wint_season)
add_df <- data.frame(wint_season = dis_yr_add_on_targets,
                     total_cases = rep(0,length(dis_yr_add_on_targets)))
west_disease_ts <- rbind(west_disease_ts,add_df) %>% arrange(wint_season)
yrs_west_disease <- west_disease_ts$wint_season
west_disease_ts <- log(west_disease_ts$total_cases + 1)

east_disease_ts <- readRDS("data/2_formatted/Disease/east_salmonellosis_ts.rds") %>% filter(wint_season %in% disease_year_targets)

dis_yr_add_on_targets <- whispers_all_yrs[!c(whispers_all_yrs %in% east_disease_ts$wint_season)] 

add_df <- data.frame(wint_season = dis_yr_add_on_targets,
                     total_cases = rep(0,length(dis_yr_add_on_targets)))
east_disease_ts <- rbind(east_disease_ts,add_df) %>% arrange(wint_season)
yrs_east_disease <- east_disease_ts$wint_season
east_disease_ts <- log(east_disease_ts$total_cases + 1)

disease_ts_unobserved <- which(yrs_west_disease %in% dis_yr_add_on_targets)
disease_ts_observed <- which(!c(yrs_west_disease %in% dis_yr_add_on_targets))

#Note, this assumes each-length disease time series from the east and west
east_disease_ts[disease_ts_unobserved] <- west_disease_ts[disease_ts_unobserved] <- 0

outbreak_yn_west <- as.numeric(west_disease_ts != 0)
outbreak_yn_east <- as.numeric(east_disease_ts != 0)

dis_obs_vect = as.numeric(!c(1:length(east_disease_ts) %in% disease_ts_unobserved))

data_mdl <- list(
  #N = ncol(dt_long_else),
  
  Nyrs = length(yrs_dt_else_obs),
  
  dt_years = whispers_all_yrs-2, 
  mast_years = whispers_all_yrs-1,
  irr_years = whispers_all_yrs,
  dis_years = whispers_all_yrs,
  cells = cells,
  
  Ncls = length(cells),
  
  Nsp = ncol(irr_dt_trim_e_ind),
  
  N_edges = nrow(ninds),
  node1 = ninds[,1],
  node2 = ninds[,2],
  
  mean_wf_val = sum(identity_matrix[,2] * fc_by_cell_filt$mean_fc),
  mean_ef_val = sum(identity_matrix[,3] * fc_by_cell_filt$mean_fc),
  
  N_dt_yrs_obs_else = length(which(yrs_dt_else_obs==1)),
  dt_obs_ind = as.array(which(yrs_dt_else_obs==1)),
  dt_miss_ind = as.array(which(yrs_dt_else_obs==0)),
  
  N_obs_mast = obs_out_else$len_y_obs_in,
  N_miss_mast = obs_out_else$len_y_mis_in,
  
  yy_obs_mast = obs_out_else$yy_data_long,
  
  mean_mast_val = mean_mast_val,
  sd_mean_mast = sd_mean_mast,
  
  dt_obs = dt_long_else,
  
  fc = fc_by_cell_filt$mean_fc,
  
  er_identity = t(identity_matrix[,2:3]),
  
  Ncells_wf = sum(identity_matrix[,2]),
  Ncells_nf = sum(identity_matrix[,3]),
  
  N_irr_obs_year = N_irr_obs_year,
  irr_obs_yrs = irr_obs_yrs,
  irr_missing_yrs = irr_missing_yrs,
  
  irr_ind_e = irr_dt_trim_e_ind,
  irr_ind_w = irr_dt_trim_w_ind,
  
  ii_obs_mast = obs_out_else$ii_obs_in,
  ii_miss_mast = obs_out_else$ii_mis_in,
  
  outbreak_yn_west = outbreak_yn_west, 
  outbreak_yn_east = outbreak_yn_east, 
  
  west_disease_ts = west_disease_ts,
  east_disease_ts = east_disease_ts,
  
  N_dis_unobserved = length(disease_ts_unobserved),
  disease_ts_observed = disease_ts_observed,
  disease_ts_unobserved = disease_ts_unobserved,
  
  dis_obs_vect = dis_obs_vect
)

saveRDS(data_mdl,paste("data/3_model_in/model_full_data_",max(disease_year_targets),".rds",sep=""))

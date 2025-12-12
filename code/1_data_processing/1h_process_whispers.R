library(readr)
library(sf)
library(usa)
library(stringr)
library(lubridate)
library(dplyr)
library(terra)
library(tibble)

disease_target_years <- 1983:2024

US_counties <- read_sf(dsn = "data/1_raw/Counties/c_05mr24.shp") %>% as.data.frame()
US_states <- usa::states

#plot(US_counties$LON,US_counties$LAT,xlim=c(-130,-65),ylim=c(24,55))
whispers_dt <- read_csv("data/1_raw/WHISPers/whispers_raw.csv")

#Change columns to date format
whispers_dt$`Event Start Date` <- as_date(whispers_dt$`Event Start Date`,format = "%m/%d/%y")
whispers_dt$`Event End Date` <- as_date(whispers_dt$`Event End Date`,format = "%m/%d/%y")

# Filter to US
whispers_dt <- whispers_dt[grepl('United States', whispers_dt$Countries),]
whispers_dt$ID <- 1:nrow(whispers_dt)


#Get Lat longs of all event locations
full_location_dt <- c()
for (each_rec in 1:nrow(whispers_dt)){
  first_row <- whispers_dt[each_rec,]
  
  
  # Because a lot of the county names don't match up exactly between WHISPers and
  # the NWS database, some  need to be fixed here:
  counties <- first_row$`Counties (or equivalent)`
  
  if(!is.na(counties)){
    counties_split <- strsplit(counties,"; ")
    counties_split_vector <- as.character(counties_split[[1]])
    
    count_rec <- c()
    state_rec <- c()
    county_rec <- c()
    count_lon <- c()
    count_lat <- c()
    
    for (each_county in 1:length(counties_split_vector)){
      split_state <- strsplit(counties_split_vector[each_county],", ")
      county_name <- split_state[[1]][1]
      county_name <- gsub("St ","St. ",county_name)
      county_name <- gsub(" County","",county_name)
      county_name <- gsub(" Parish","",county_name)
      county_name <- gsub(" City and Borough","",county_name)
      county_name <- gsub(" Borough","",county_name)
      county_name <- gsub(" Census Area","",county_name)
      county_name <- gsub(" Municipality","",county_name)
      
      if (county_name == "Lac Qui Parle"){
        county_name = "Lac qui Parle"
      }
      if (county_name == "Lake Of The Woods"){
        county_name = "Lake of the Woods"
      }
      if (county_name == "Lewis And Clark"){
        county_name = "Lewis and Clark"
      }
      
      state_code <- split_state[[1]][2]
      #Add min because there are some duplicates.
      county_ind <- min(which(US_counties$STATE == state_code & US_counties$COUNTYNAME == county_name))
      location_coords <- c(US_counties$LON[county_ind],US_counties$LAT[county_ind])
      
      #For some reason one county is missing!
      if (county_name == "Monroe" & state_code == "FL"){
        location_coords <- c(-81.02,25.48)
      }
      
      count_rec <- c(count_rec,each_rec)
      state_rec <- c(state_rec,state_code)
      county_rec <- c(county_rec,county_name)
      count_lon <- c(count_lon,location_coords[1])
      count_lat <- c(count_lat,location_coords[2])
    }
    df_ind <- data.frame(ID=count_rec,state_rec=state_rec,county_rec=county_rec,
                         count_lon=count_lon,count_lat=count_lat)
    full_location_dt <- rbind(full_location_dt,df_ind)
  } else {
    # If no counties recorded (just state)
    state_abb <- US_states$abb[which(US_states$name == first_row$`States (or equivalent)`)]
    
    df_ind <- data.frame(ID=each_rec,state_rec=state_abb,county_rec=NA,
                         count_lon=US_states$long[which(US_states$abb == state_abb)],
                         count_lat=US_states$lat[which(US_states$abb == state_abb)])
    full_location_dt <- rbind(full_location_dt,df_ind)
  }
  
}
#View(full_location_dt)
# 
# 
# all_spec_list <- c()
# for (each_rec in 1:nrow(whispers_dt)){
#   species_check <- whispers_dt$Species[each_rec]
#   species_split <- strsplit(species_check,"; ")
#   spec_vect <- as.character(species_split[[1]])
#   all_spec_list <- unique(c(all_spec_list,spec_vect))
# }
# 
# all_spec_names <- data.frame(ID=1:length(all_spec_list),spec=all_spec_list)
# write_csv(all_spec_names,"all_spec_list_1_18_24.csv")

ordered_spec_list <- read_csv("data/1_raw/Other/animal_species.csv")

passerine_entries <- filter(ordered_spec_list,Passerine==1)
whispers_dt$passerines_involved <- NA
for (each_rec in 1:nrow(whispers_dt)){
  species_check <- whispers_dt$Species[each_rec]
  species_split <- strsplit(species_check,"; ")
  spec_vect <- as.character(species_split[[1]])
  yn_passerines <- 0
  for (each_vect_entry in spec_vect){
    yn_passerines <- yn_passerines + sum(each_vect_entry %in% passerine_entries$spec) 
  }
  if(yn_passerines > 0){
    whispers_dt$passerines_involved[each_rec] <- 1
  } else {
    whispers_dt$passerines_involved[each_rec] <- 0
  }
}

whispers_dt_passerines <- whispers_dt %>% filter(passerines_involved == 1)

whispers_dt_passerines$event_mid_date <- whispers_dt_passerines$`Event Start Date` + (whispers_dt_passerines$`Event End Date` - whispers_dt_passerines$`Event Start Date`)

whispers_dt_passerines$state_rec <- state.abb[match(whispers_dt_passerines$`States (or equivalent)`,state.name)]

#Set up winter season calculations - if no mid date (1 data point), just use start date
whispers_dt_passerines$event_mid_date[is.na(whispers_dt_passerines$event_mid_date)] <- whispers_dt_passerines$`Event Start Date`[is.na(whispers_dt_passerines$event_mid_date)]
whispers_dt_passerines$month_mid <- month(whispers_dt_passerines$event_mid_date)
whispers_dt_passerines$year_mid <- year(whispers_dt_passerines$event_mid_date)

#Filter to winter records
whispers_dt_passerines <- whispers_dt_passerines %>% filter(month_mid %in% c(10,11,12,1,2,3,4))

whispers_dt_passerines$wint_season <- whispers_dt_passerines$year_mid
whispers_dt_passerines$wint_season[which(whispers_dt_passerines$month_mid %in% c(10,11,12))] <- whispers_dt_passerines$wint_season[which(whispers_dt_passerines$month_mid %in% c(10,11,12))] + 1


### Read in each record, and sort out whether it was in the Western or Eastern Region ####
us_epa_ecoregions <- vect("data/1_raw/Ecoregions/NA_CEC_Eco_Level1.shp")
us_epa_ecoregions <- terra::project(us_epa_ecoregions,"epsg:4326")
us_epa_ecoregions_trim <- simplifyGeom(us_epa_ecoregions,tolerance=.01) # speed things up
table(us_epa_ecoregions_trim$NA_L1NAME)

pts <- vect(cbind(full_location_dt$count_lon,full_location_dt$count_lat), crs="+proj=longlat")

er_out <- terra::extract(us_epa_ecoregions_trim, pts)
full_location_dt$er_extract <- er_out$NA_L1NAME

west_super_region_list <- c("MARINE WEST COAST FOREST","MEDITERRANEAN CALIFORNIA",
                            "NORTH AMERICAN DESERTS","NORTHWESTERN FORESTED MOUNTAINS",
                            "SOUTHERN SEMIARID HIGHLANDS","TEMPERATE SIERRAS")

east_super_region_list <- c("EASTERN TEMPERATE FORESTS","GREAT PLAINS",
                            "NORTHERN FORESTS","TROPICAL WET FORESTS")

full_location_dt$ew_superregion <- NA
full_location_dt$ew_superregion[which(full_location_dt$er_extract %in% west_super_region_list)] <- "West"
full_location_dt$ew_superregion[which(full_location_dt$er_extract %in% east_super_region_list)] <- "East"

#Some of the island/coastal counties have centerpoints over water, so adjust those here
full_location_dt$ew_superregion[which(is.na(full_location_dt$er_extract) & is.na(full_location_dt$ew_superregion) & full_location_dt$state_rec %in% c("WA","AK"))] <- "West"
full_location_dt$ew_superregion[which(is.na(full_location_dt$er_extract) & is.na(full_location_dt$ew_superregion) & full_location_dt$state_rec%in% c("MD","GA","TX","QC","NC","Dundas & Glengarry United"))] <- "East"

#Check everything went as planned
table(full_location_dt$er_extract,full_location_dt$ew_superregion)
table(full_location_dt$state_rec,full_location_dt$ew_superregion)
table(full_location_dt$ID,full_location_dt$ew_superregion)

#Now for each event record, sort into west and east
east_salmonellosis_events <- c()
west_salmonellosis_events <- c()
for (nn in whispers_dt_passerines$ID){
  row_ind <- which(whispers_dt_passerines$ID == nn)
  event_locs <- filter(full_location_dt,ID == nn)
  if("East" %in% event_locs$ew_superregion){
    to_add_east <- whispers_dt_passerines[row_ind,]
    to_add_east <- merge(to_add_east,event_locs,by="ID")
    east_salmonellosis_events <- rbind(east_salmonellosis_events,to_add_east)
  }
  #Note this is mutually non-exlusive
  if ("West" %in% event_locs$ew_superregion){
    to_add_west <- whispers_dt_passerines[row_ind,]
    to_add_west <- merge(to_add_west,event_locs,by="ID")
    west_salmonellosis_events <- rbind(west_salmonellosis_events,to_add_west)
  }
}
table(west_salmonellosis_events$`States (or equivalent)`)
table(east_salmonellosis_events$`States (or equivalent)`)

#Combine events for saving
full_dt_loc <- rbind(west_salmonellosis_events,east_salmonellosis_events)
full_dt_loc$PISI_involved <- grepl("Pine Siskin", full_dt_loc$Species, fixed = TRUE)

saveRDS(full_dt_loc,"data/2_formatted/Disease/disease_locations.rds")

irr_yr_template <- data.frame(wint_season = 1984:2025)
west_salmonellosis_ts <- west_salmonellosis_events %>% 
  group_by(wint_season) %>% 
  dplyr::select(ID,wint_season,`Number Affected`) %>%
  filter(wint_season %in% disease_target_years) %>% 
  unique() %>%
  summarise(total_cases = sum(`Number Affected`))

west_salmonellosis_ts <- merge(irr_yr_template,west_salmonellosis_ts,all=TRUE)
west_salmonellosis_ts$total_cases[which(is.na(west_salmonellosis_ts$total_cases))] <- 0

east_salmonellosis_ts <- east_salmonellosis_events %>% 
  group_by(wint_season) %>% 
  dplyr::select(ID,wint_season,`Number Affected`) %>%
  filter(wint_season %in% disease_target_years) %>% 
  unique() %>%
  summarise(total_cases = sum(`Number Affected`))

east_salmonellosis_ts <- merge(irr_yr_template,east_salmonellosis_ts,all=TRUE)
east_salmonellosis_ts$total_cases[which(is.na(east_salmonellosis_ts$total_cases))] <- 0

plot(west_salmonellosis_ts$wint_season,log(west_salmonellosis_ts$total_cases+1),type="l",col="dodgerblue4")
points(east_salmonellosis_ts$wint_season,log(east_salmonellosis_ts$total_cases+1),type="l",col="forestgreen")

plot(log(west_salmonellosis_ts$total_cases+1),log(east_salmonellosis_ts$total_cases+1))

saveRDS(west_salmonellosis_ts,"data/2_formatted/Disease/west_salmonellosis_ts.rds")
saveRDS(east_salmonellosis_ts,"data/2_formatted/Disease/east_salmonellosis_ts.rds")

#### Get breakdown by species ####
west_salmonellosis_species <- west_salmonellosis_events %>% 
  dplyr::select(ID,wint_season,Species,`Number Affected`) %>% unique() %>%
  filter(wint_season %in% disease_target_years)
print(nrow(west_salmonellosis_species))

all_spec <- c()
for (each_row in 1:nrow(west_salmonellosis_species)){
  spec_row <- west_salmonellosis_species$Species[each_row]
  spec_row <- str_split(spec_row,"; ")
  for (each_spec in 1:length(spec_row[[1]])){
    all_spec <- c(all_spec,spec_row[[1]][each_spec])
  }
}

spec_freq_w <- as.data.frame(table(all_spec))
spec_freq_w$prop_events <- spec_freq_w$Freq/nrow(west_salmonellosis_species)

any_irr_species <- c()
irr_species_list <- c("Pine Siskin","Red-breasted Nuthatch","Evening Grosbeak","Pine Grosbeak",
                      "Purple Finch","Red Crossbill","White-winged Crossbill","Common Redpoll")
for (each_row in 1:nrow(west_salmonellosis_species)){
  spec_row <- west_salmonellosis_species$Species[each_row]
  spec_row <- str_split(spec_row,"; ")
  any_irr_species <- c(any_irr_species,as.numeric(sum(unlist(spec_row) %in% irr_species_list) > 0))
}
west_salmonellosis_species$any_irr_spec <- any_irr_species
irr_involved_rate_w <- sum(west_salmonellosis_species$any_irr_spec)/nrow(west_salmonellosis_species) 

irr_involved_rate_w_large <- west_salmonellosis_species$any_irr_spec[west_salmonellosis_species$`Number Affected` >= 200]
sum(irr_involved_rate_w_large)/length(irr_involved_rate_w_large)

#East
east_salmonellosis_species <- east_salmonellosis_events %>% 
  dplyr::select(ID,wint_season,Species,`Number Affected`) %>% unique() %>%
  filter(wint_season %in% disease_target_years)
print(nrow(east_salmonellosis_species))

all_spec <- c()
for (each_row in 1:nrow(east_salmonellosis_species)){
  spec_row <- east_salmonellosis_species$Species[each_row]
  spec_row <- str_split(spec_row,"; ")
  for (each_spec in 1:length(spec_row[[1]])){
    all_spec <- c(all_spec,spec_row[[1]][each_spec])
  }
}

spec_freq_e <- as.data.frame(table(all_spec))
spec_freq_e$prop_events <- spec_freq_e$Freq/nrow(east_salmonellosis_species)

any_irr_species <- c()
for (each_row in 1:nrow(east_salmonellosis_species)){
  spec_row <- east_salmonellosis_species$Species[each_row]
  spec_row <- str_split(spec_row,"; ")
  any_irr_species <- c(any_irr_species,as.numeric(sum(unlist(spec_row) %in% irr_species_list) > 0))
}
east_salmonellosis_species$any_irr_spec <- any_irr_species
irr_involved_rate_e <- sum(east_salmonellosis_species$any_irr_spec)/nrow(east_salmonellosis_species) 

irr_involved_rate_e_large <- east_salmonellosis_species$any_irr_spec[east_salmonellosis_species$`Number Affected` >= 200]
sum(irr_involved_rate_e_large)/length(irr_involved_rate_e_large)
write_csv(spec_freq_w,"data/2_formatted/Disease/spec_freq_w.csv") 
write_csv(spec_freq_e,"data/2_formatted/Disease/spec_freq_e.csv") 

# Cross-regional outbreak size 
up_90th <- quantile(c(east_salmonellosis_species$`Number Affected`,west_salmonellosis_species$`Number Affected`),c(.9))

comb_large_outbreaks <- c(irr_involved_rate_e_large,irr_involved_rate_w_large)
sum(comb_large_outbreaks)/length(comb_large_outbreaks)

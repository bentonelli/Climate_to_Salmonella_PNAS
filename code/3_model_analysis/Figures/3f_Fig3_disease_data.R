#Script to take disease data and plot

library(dplyr)
library(terra)
library(sf)
library(ggplot2)
library(tidyr)
library(stringr)
library(readr)

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

#Remove two canadian points
dis_loc_data <- dis_loc_data[complete.cases(dis_loc_data$count_lon),]

dis_loc_data$lon_jitter <- dis_loc_data$count_lon + rnorm(length(dis_loc_data$count_lon),0,.1)
dis_loc_data$lat_jitter <- dis_loc_data$count_lat + rnorm(length(dis_loc_data$count_lat),0,.1)

projection_code <- 5070

d_points <- data.frame(long = dis_loc_data$lon_jitter, 
                       lat  = dis_loc_data$lat_jitter) %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  st_transform(crs = projection_code)

regional_sf <- sf_regions %>%
  st_transform(crs = projection_code)

ggplot() + 
  geom_sf(data=regional_sf,aes(fill=NA_SR),alpha=.15,col="darkgrey",lwd=.2) + 
  scale_fill_manual(values = c("#E98A15","grey","#648DE5")) +
  geom_sf(data=d_points,aes(col = dis_loc_data$ew_superregion),alpha=.5,size=3) +
  #geom_point(data=d_points,cex=2,alpha=.5,pch=19) + 
  scale_color_manual(values = c("#E98A15","#648DE5","transparent")) +
  theme_classic() + 
  theme(line = element_blank(),
        text = element_blank(),
        title = element_blank()) +
  theme(legend.position="none") + 
  xlab("") + ylab("")

### BY species ####

spec_freq_w <- read_csv("data/2_formatted/Disease/spec_freq_w.csv")
spec_freq_w$Region <- "West"
spec_freq_e <- read_csv("data/2_formatted/Disease/spec_freq_e.csv")
spec_freq_e$Region <- "East"

spec_freq_comb_irr <- rbind(spec_freq_e,spec_freq_w) %>% 
  dplyr::select(-Freq) %>% filter(all_spec %in% c("Pine Siskin",
                                                  "Red Crossbill",
                                                  "White-winged Crossbill",
                                                  "Red-breasted Nuthatch",
                                                  "Evening Grosbeak",
                                                  "Pine Grosbeak",
                                                  "Common Redpoll",
                                                  "Purple Finch"))

spec_freq_comb_irr$all_spec[which(spec_freq_comb_irr$all_spec == "Common Redpoll")] <- "Redpoll"

#Some species don't have any records, so manually input those zeros here
spec_freq_comb_irr <- rbind(spec_freq_comb_irr,data.frame(all_spec = c("Red-breasted Nuthatch",
                                                                       "Red-breasted Nuthatch",
                                                                       "White-winged Crossbill",
                                                                       "White-winged Crossbill",
                                                                       "Red Crossbill",
                                                                       "Pine Grosbeak"),
                                                          prop_events = c(0,0,0,0,0,0),
                                                          Region = c("West","East","West","East","West","West")))
spec_freq_comb_irr$variable = factor(spec_freq_comb_irr$all_spec, 
                                     levels = c("Pine Siskin","Evening Grosbeak",
                                                "Redpoll","Purple Finch",
                                                "Pine Grosbeak","Red Crossbill",
                                                "Red-breasted Nuthatch","White-winged Crossbill"), 
                                     ordered = TRUE)
spec_freq_comb_irr$Region = factor(spec_freq_comb_irr$Region, 
                                   levels = c("West","East"), 
                                   ordered = TRUE)

spec_freq_comb_irr$perc_events <- spec_freq_comb_irr$prop_events*100

pp_perc <- ggplot(data=spec_freq_comb_irr,aes(fill=Region,y=perc_events,x=variable)) + 
  geom_bar(stat="identity", colour = "black", position = position_dodge(width = .7), 
           width = 0.7,alpha=.9)+
  scale_fill_manual(values = rev(c("#E98A15","#648DE5"))) +
  scale_x_discrete(labels = c("Pine siskin","Evening grosbeak","Redpoll","Purple finch",
                              "Pine grosbeak","Red crossbill","Red-breasted nuthatch","White-winged crossbill"), name = "Group") +
  geom_text(aes(label = round(perc_events,0)),size=6, 
            position = position_dodge(width = .7), hjust = .5,vjust=-.5)+
  theme_grey(base_size=16) +
  ylim(0, 100) +
  theme(legend.position=c(.80,.85),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA),
        legend.key.size = unit(1, 'cm'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,color="black"),
        axis.text.y = element_text(color="black"),
        axis.line = element_line(color='black',size=1),
        text=element_text(size=20,  family="Arial"),
        panel.grid.major = element_line(color="grey95"),
        panel.background = element_rect(fill="white"))
pp_perc

### By Year ####
west_disease_ts <- readRDS("data/2_formatted/Disease/west_salmonellosis_ts.rds") %>% filter(wint_season <= 2024 & wint_season >= 1988)
east_disease_ts <- readRDS("data/2_formatted/Disease/east_salmonellosis_ts.rds") %>% filter(wint_season <= 2024 & wint_season >= 1988)

w_transition <- west_disease_ts
w_transition$region="West"

e_transition <- east_disease_ts
e_transition$region="East"

comb_disease_ts <- rbind(w_transition,e_transition)
comb_disease_ts$region <- factor(comb_disease_ts$region, levels = c("West","East"))
par(las=1)
par(mar=c(4,6,2,2))

y_log_scale <- c(0,10,100,1000,10000)+1
# Line plot in ggplot
ggplot(comb_disease_ts) +
  geom_line(aes(x=wint_season,y=log(total_cases+1),col=region),
            size=3,alpha=.95) +
  scale_color_manual(values = c("#648DE5","#E98A15")) +
  xlim(1988,2024) + 
  scale_y_continuous(breaks=log(y_log_scale),name = "",
                     labels = y_log_scale-1,limits = c(0,11))+
  theme(legend.position=c(.5,.9),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA),
        legend.key.size = unit(2, 'cm'),
        legend.direction="horizontal",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.x = element_text(angle = 0, vjust = 0, hjust=0,color="black"),
        axis.text.y = element_text(color="black"),
        axis.line = element_line(color='black',size=1),
        text=element_text(size=20,  family="Arial"),
        panel.grid.major = element_line(color="grey95"),
        panel.background = element_rect(fill="white"))


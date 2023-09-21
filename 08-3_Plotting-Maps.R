###################################################

# Script for creating maps for the paper
# Have several different files that need to be 
# loaded

###################################################
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(terra)
library(stringr)
library(forcats)
library(RColorBrewer)
library(geosphere)
library(forcats)
library(ggridges)
library(grid)
library(gridExtra)
library(gtable)
library(purrr)
library(matrixStats)
library(sfnetworks)
library(rcartocolor)
library(ggnewscale)
library(rgeos)
library(rnaturalearth)
library(ggpattern)
library(patchwork)

rm(list = ls())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
m_dir <- paste(working.dir, "Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")
pop_dir <-  paste(working.dir, "Output_Population", sep="/")
sim_dir <-  paste(working.dir, "Simulations", sep="/")


## Functions
setwd(working.dir)
source("X_Functions.R")



model.name <- "ningaloo"

colours <- c("#69BE28", "#005594", "#8AD2D8", "#53AF8B")
a4.width <- 160

#### READ IN DATA ####                                       
#* For mapping ####

# Set CRS for transformations
wgscrs <- "+proj=longlat +datum=WGS84"
gdacrs <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
sppcrs <- CRS("+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs")       # crs for sp objects

# Set cropping extent - larger than most zoomed out plot
e <- ext(112, 120.0, -28, -16)

# Load necessary spatial files

sf_use_s2(F)                                                                    # Switch off spatial geometry for cropping

# Australian outline and state and commonwealth marine parks
setwd(sp_dir)
aus <- st_read("cstauscd_r.mif") %>%                    # Geodata 100k coastline available: https://data.gov.au/dataset/ds-ga-a05f7892-eae3-7506-e044-00144fdd4fa6/
  dplyr::filter(FEAT_CODE %in% c("mainland", "island"))
st_crs(aus) <- gdacrs
ausc <- st_crop(aus, e)

# Commonwealth parks
aumpa  <- st_read("AustraliaNetworkMarineParks.shp")                            # All aus mpas
mpa <- st_crop(aumpa, e)                                                        # Crop to the study area
# Reorder levels so everything plots nicely
unique(mpa$ZoneName)
mpa$ZoneName <- factor(mpa$ZoneName, levels = c("Multiple Use Zone", 
                                                "Recreational Use Zone",
                                                "Habitat Protection Zone",
                                                "National Park Zone"))
npz <- mpa[mpa$ZoneName %in% "National Park Zone", ]                            # Just National Park Zones
plot(mpa$geometry)


# State parks
wampa <- st_read("WA_MPA_2020.shp")
st_crs(wampa) <- gdacrs
# Simplify names for plot legend
wampa$waname <- gsub("( \\().+(\\))", "", wampa$ZONE_TYPE)
wampa$waname <- gsub(" [1-4]", "", wampa$waname)
wampa$waname[wampa$NAME == "Hamelin Pool"]     <- "Marine Nature Reserve"
wampa$waname[wampa$NAME == "Abrolhos Islands"] <- "Fish Habitat Protection Area"
wampa$waname <- dplyr::recode(wampa$waname, 
                              "General Use" = "General Use Zone",
                              "Special Purpose Zone (Shore Based Activities)" = 
                                "Special Purpose Zone\n(Shore Based Activities)",
                              "Special Purpose Zone (Seagrass Protection) (IUCN IV)" = 
                                "Special Purpose Zone",
                              "MMA" = 'Marine Management Area' )

wampa <- wampa %>% 
  mutate(waname = ifelse(NAME %in% "Barrow Island" & TYPE %in% "Marine Park", "Sanctuary Zone", waname)) %>% 
  mutate(waname = ifelse(NAME %in% "Barrow Island" & TYPE %in% "Marine Management Area", "Marine Management Area", waname)) %>% 
  mutate(waname = ifelse(NAME %in% "Shark Bay" & ZONE_TYPE %in% "Recreation Zone (IUCN II)", "Recreation Area", waname)) %>% 
  mutate(waname = ifelse(NAME %in% "Hamelin Pool" & TYPE %in% "Marine Nature Reserve", "Sanctuary Zone", waname)) %>% 
  filter(!NAME %in% c("Miaboolya Beach"))

wampa$waname <- factor(wampa$waname, levels = c("Unassigned", 
                                                "Marine Management Area",
                                                "Conservation Area",
                                                "General Use Zone",
                                                "Recreation Area",
                                                "Sanctuary Zone",
                                                "Special Purpose Zone"))

wampa <- st_crop(wampa, e)                                                      # Crop to the study area
wasanc <- wampa[wampa$ZONE_TYPE %in% "Sanctuary Zone (IUCN IA)", ]
plot(wampa$geometry)

# Terrestrial parks
# terrnp <- st_read("Legislated_Lands_and_Waters_DBCA_011.shp") %>%  # Terrestrial reserves
#   dplyr::filter(leg_catego %in% c("Nature Reserve", "National Park"))
# terrnp <- st_crop(terrnp, e)       # Crop to the study area - using a different extent as this is on land

# Key Ecological Features
# kef <- st_read("AU_DOEE_KEF_2015.shp")
# kef <- st_crop(kef, e)                                                          # Crop
# unique(kef$NAME)
# # Simplify names for plot legend
# unique(kef$NAME)
# kef$NAME <- dplyr::recode(kef$NAME,
#                           "Ancient coastline at 125 m depth contour" = "Ancient coastline",
#                           "Continental Slope Demersal Fish Communities" = "Continental slope fish",
#                           "Canyons linking the Cuvier Abyssal Plain and the Cape Range Peninsula" = "Cuvier Abyssal Plain canyons",
#                           "Commonwealth waters adjacent to Ningaloo Reef" = "Ningaloo Reef")
# # Reorder levels so everything plots nicely
# # kef$NAME <- factor(kef$NAME, levels = c("Western rock lobster", "Western demersal fish", "Wallaby Saddle", 
# #                                         "Abrolhos Islands", "Ancient coastline", 
# #                                         "West coast canyons", "West coast lagoons"))

# Coastal waters limit
cwatr <- st_read("amb_coastal_waters_limit.shp")       # Coastal waters limit
cwatr <- st_crop(cwatr, e)

# Bathymetry data
cbathy <- raster("bath_250_good.tif")                    
bath_r <- rast(cbathy)
crs(bath_r) <- wgscrs
bath_r <- crop(bath_r, ext(112, 120.0, -28, -16))
bath_df <- as.data.frame(bath_r, xy = T, na.rm = T)                             # Dataframe - cropped and above 0 use for bath cross section
bath_r <- clamp(bath_r, upper = 0, value = F)                               # Only data below 0
bathy <- as.data.frame(bath_r, xy = T, na.rm = T)

# Fisheries layers
CEO_notices <- st_read("Fisheries_Guide_CEO_Notices_Determinations_DPIRD_060.shp") %>% 
  filter(ufi %in% c( 912, 913, 914, 916, 923, 917, 918, 920, 921, 922, 924, 925)) %>% 
  dplyr::select(-notc_detr)

Management_plans <- st_read("Fisheries_Guide_Consolidated_Management_Plans_DPIRD_062.shp") %>% 
  filter(ufi %in% c(699, 700, 563, 571, 568, 572, 554, 694, 695))

plot(Management_plans$geometry)
plot(CEO_notices$geometry, col="blue", add=T)

License_conditions <- st_read("Fisheries_Guide_Licence_Conditions_DPIRD_050.shp") %>% 
  filter(ufi %in% 13) %>% 
  dplyr::select(-object_id)

fisheries_closures <- rbind(CEO_notices, Management_plans, License_conditions) %>% 
  filter(!descript %in% c("Schedule 1 - Description of the Fishery", "Schedule - Item 1 (The Fishery)"))

fisheries_closures <- st_crop(fisheries_closures, e)

fisheries_closures$name <- dplyr::recode(fisheries_closures$name,
                                         "Onslow Prawn Managed Fishery Management Plan 1991 - Notice of Areas Closed to Fishing for Prawns" = "Onslow Prawn Fishery",
                                         "Pilbara Fish Trawl Interim Managed Fishery Management Plan 1997"  = "Pilbara Fish Trap and Trawl Fishery",
                                         "Exmouth Gulf Prawn Limited Entry Fishery Notice 1989" = "Exmouth Gulf Prawn Fishery",
                                         "Point Maud to Tantabiddi Closure" = "Point Maud to Tantabiddi Closure",
                                         "Shark Bay Prawn Managed Fishery Management Plan 1993 - Determination of Areas Closed to Fishing for Prawns" = "Shark Bay Prawn Fishery",
                                         "Gascoyne Demersal Scalefish Management Plan 2010" = "Gascoyne Demersal Scalefish Fishery")

fisheries_closures <- fisheries_closures %>% 
  mutate(name = ifelse(name %in% "Pilbara Fish Trap and Trawl Fishery" & descript %in% "Schedule 3 (Item 3 - Zone 2 - Area 3)", "Pilbara Fish Trap, Trawl and Line Fishery", name)) %>%
  mutate(pattern_dir = ifelse(name %in% c("Pilbara Fish Trap, Trawl and Line Fishery","Exmouth Gulf Prawn Fishery", "Shark Bay Prawn Fishery" ), "Left", "Right"))

## Make overview map
nmpa_fills <- scale_fill_manual(values = c("National Park Zone" = "#7bbc63",
                                           "Multiple Use Zone" = "#b9e6fb",
                                           "Recreational Use Zone" = "#ffb36b",
                                           "Habitat Protection Zone" = "#fff8a3"
), 
name = "Australian Marine Parks")

wampa_fills <- scale_fill_manual(values = c("Marine Management Area" = "#b7cfe1",
                                            "Conservation Area" = "#b3a63d",
                                            "Sanctuary Zone" = "#bfd054",
                                            "General Use Zone" = "#bddde1",
                                            "Recreation Area" = "#f4e952",
                                            "Special Purpose Zone" = "#c5bcc9",
                                            "Marine Nature Reserve" = "#bfd054"
),
name = "State Marine Parks")
# Grey area off barrow island have a SZ and a marine management area

# terr_fills <- scale_fill_manual(values = c("National Park" = "#c4cea6",          # Set the colours for terrestrial parks
#                                            "Nature Reserve" = "#e4d0bb"))
closed_fills <- scale_pattern_color_manual(values= c("Pilbara Fish Trap and Trawl Fishery" = "deeppink",
                                            "Pilbara Fish Trap, Trawl and Line Fishery" = "blue4",
                                            "Exmouth Gulf Prawn Fishery" = "darkviolet",
                                            "Point Maud to Tantabiddi Closure" = "goldenrod2",
                                            "Shark Bay Prawn Fishery" = "red",
                                            "Gascoyne Demersal Scalefish Fishery" = "chartreuse4"),
                                           name = "Fishery Closures")
closure_pattern <- scale_pattern_angle_manual(values = c(-30, 30), guide="none")

# nmpa <- mpa %>%
#   dplyr::filter(ResName %in% c("Ningaloo", "Shark Bay"))

# gmpa <- mpa %>%
#   dplyr::filter(ResName %in% "Gascoyne")

p3 <- ggplot() +
  geom_contour_filled(data = bathy, aes(x = x, y = y, z = bath_250_good,
                                        fill = after_stat(level)),
                      breaks = c(0, -30, -70, -200, - 700, -2000 , -4000,-6000)) +
  geom_contour(data = bathy, aes(x = x, y = y, z = bath_250_good),
               breaks = c(-30, -70, -200, - 700, -2000 , -4000,-6000), colour = "white", alpha = 3/5, size = 0.1) +
  scale_fill_grey(start = 1, end = 0.5, guide = "none") +
  geom_sf(data = ausc, fill = "seashell2", colour = "grey80", size = 0.1) +
  new_scale_fill() +
  geom_sf_pattern(data=fisheries_closures, aes(pattern_colour = name, pattern_angle=pattern_dir),
                  pattern= "stripe" ,pattern_size=0.1, pattern_spacing=0.01, pattern_alpha=0.5,
                  colour=NA, alpha=0)+
  closed_fills + 
  closure_pattern +
  new_scale_fill() +
  geom_sf(data = wampa, aes(fill = waname), alpha = 2/5, colour = NA) +
  wampa_fills +
  labs(fill = "State Marine Parks") +
  new_scale_fill() +
  labs(fill = "Terrestrial Managed Areas") +
  new_scale_fill() +
  geom_sf(data = mpa, aes(fill = ZoneName), alpha = 0.4, colour = NA) +
  nmpa_fills + 
  labs(fill = "Australian Marine Parks") +
  new_scale_fill() +
  geom_sf(data = cwatr, colour = "firebrick", alpha = 4/5, size = 0.4) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_legend(order = 1)) +
  annotate(geom = "text", x = c((114.1279 + 0.4), (113.6775 + 0.45)), 
           y = c(-21.9323, -22.7212), label = c("Exmouth", "Pt Cloates"),
           size = 3) +
  annotate(geom = "point", x = c(114.1279, 113.6775), 
           y = c(-21.9323, -22.7212)) +
  coord_sf(xlim = c(112, 120.00), ylim = c(-26.5, -17.75)) +
  theme_minimal() +
  theme(legend.justification = "top",
        legend.box.margin = margin(c(160,0,0,0)),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.spacing = unit(0.01, "cm"))
p3 

p3.1 <- ggplot() +
  geom_sf(data = aus, fill = "seashell2", colour = "grey90", size = 0.075) +
  #geom_sf(data = aumpa, alpha = 5/6, colour = "grey85", size = 0.02) +
  coord_sf(xlim = c(108, 125), ylim = c(-37, -13)) +
  annotate(geom = "text", x=c(110), y=c(-29.4), label = c("Indian\nOcean"), size=4)+
  annotate(geom = "text", x=c(120.75), y=c(-25.94), label = c("Western\nAustralia"), size=4)+
  annotate("rect", xmin = 113, xmax = 114.35, ymin = -22.8, ymax = -21.5,   # Change here 
           colour = "grey25", fill = "white", alpha = 1/5, size = 0.2) +
  theme_bw() +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "grey70"))+
  ylab(NULL)+
  xlab(NULL)
p3.1

p3 + inset_element(p3.1, left = -1.25, right = 2.874, top = 0.4, bottom = 0)  

setwd(fig_dir)
ggsave('broad-site-plot.png', dpi = 200, width = 10, height = 10)


#### PLOTS OF THE MODEL ####
setwd(sg_dir)
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
water <- readRDS(paste0(model.name, sep="_","water"))

setwd(sp_dir)
bathy <- raster("ga_bathy_ningaloocrop.tif")
WHA <- st_read("2013_02_WorldHeritageMarineProgramme.shp") %>% 
  st_transform(4283)%>%
  st_make_valid %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) 
BR <- st_read("Boat_Ramps.shp") %>% 
  st_transform(4283)%>%
  st_make_valid() 
BR <- BR[1:4,]

network <- st_read(paste0(model.name, sep="_","network.shapefile.shp"))

NCELL <- nrow(water)

setwd(sp_dir)

MP <- st_read("WA_MPA_2018.shp")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) %>% 
  mutate(Year.Sanct = 2005)

NTZ <- MP%>%
  filter(IUCN == "IA") %>% 
  filter(!COMMENTS == "Cloates Sanctuary Zone") %>% # Cloates is included in the Australian marine parks shapefile because I couldn't get them to line up nicely using the State MP shapefile
  rename(Name = "COMMENTS") %>% 
  filter(!Name %in% c("Previously PA_ID WA_42756 in terrestrial CAPAD", "Conservation Area")) %>% 
  dplyr::select(Name, Year.Sanct, geometry) 

AMP <- st_read("NTZ and Fished areas for status.shp") %>% 
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5)
plot(AMP)


AMP_NTZ <- AMP %>% 
  filter(ResName == "Ningaloo"|COMMENTS == "Cloates Sanctuary Zone") %>% # Select just Cloates
  mutate(Year.Sanct = ifelse(ResName %in% c("Ningaloo"), 2017, 2005)) %>% 
  rename(Name = "COMMENTS") %>% 
  mutate(Name = ifelse(is.na(Name), "Comm Cloates", Name)) %>% 
  dplyr::select(Name, Year.Sanct, geometry) 
plot(AMP_NTZ$geometry)

# Put all of the NTZs together into one object
NTZ <- rbind(NTZ, AMP_NTZ) # Put all the different NTZs together into one object
NTZ <- st_make_valid(NTZ) 

#* Create list of cells to restrict plots to shallow water (<30m)
water_points <- st_centroid_within_poly(water) 

water_bathy <- raster::extract(bathy, water_points, fun=mean, df=TRUE)

water_bathy <- water_bathy %>% 
  mutate(ID = as.factor(ID))

model_WHA <- water %>% 
  st_intersects(., WHA) %>% 
  as.data.frame()

water_WHA <-water[c(as.numeric(model_WHA$row.id)), ]

water <- water_WHA %>% 
  mutate(ID = as.factor(ID)) %>% 
  left_join(., water_bathy, by="ID") %>% 
  rename(bathy = "ga_bathy_ningaloocrop") %>% 
  #filter(bathy >= c(-30)) %>% 
  filter(!is.na(bathy))


#### MAP OF THE NINGALOO MODEL ####
water <- water %>% 
  mutate(WHA = ifelse(ID %in% c(model_WHA$row.id), "Y", "N")) %>% 
  mutate(Map_Colour = ifelse(Fished_2017=="N", "NTZ", ifelse(WHA=="Y", "WHA", "None")))


map <- water %>% 
  filter(!ID==387) %>% # Weird cell on the side that makes things look odd
  ggplot(.)+
  geom_sf(aes(fill=Map_Colour), colour="grey20", lwd=0.2)+
  scale_fill_manual(values=c("NTZ"= "#48A02D", "WHA"="#D6CF7D", "None"="skyblue1"),
                    labels = c("No-take zone", "World heritage area", "Outside world heritage\nand marine park area"),
                    name="Zone type")+
  theme_void() +
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(1.5,"line")) 
map

setwd(fig_dir)
ggsave(map, filename="Whole_map.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )


#### CATCH PER UNIT EFFORT MAP ####

Names <- c("Historical and Current Management", "No Spatial Management", 
           "Temporal and Spatial Management","Temporal Management Only" )

Effort_Scen <- list()
Spatial_Qs <- list()

setwd(sg_dir)
Effort_Scen[[1]] <- readRDS(paste0(model.name, sep="_", "fishing_test"))
Spatial_Qs[[1]] <- readRDS(paste0(model.name, sep="_", "Spatial_q_NTZ"))
Spatial_Qs[[2]] <- readRDS(paste0(model.name, sep="_", "Spatial_q_No_NTZ"))

setwd(sim_dir)
Effort_Scen[[2]] <- readRDS(paste0(model.name, sep="_", "burn_in_fishing_High_M"))
Effort_Scen[[3]] <- readRDS(paste0(model.name, sep="_", "S02_fishing"))
Effort_Scen[[4]] <- readRDS(paste0(model.name, sep="_", "S03_fishing"))

#* Convert back to effort in boat days ####

Boat_Days <- list()

Boat_Days_Scen <- array(0, dim=c(NCELL, 12, 59))

for(S in 1:4){
  
  if(S==1|S==3){
    for (YEAR in 1:59){
      Boat_Days_Scen[,,YEAR] <- Effort_Scen[[S]][,,YEAR] / Spatial_Qs[[1]][,YEAR]
    }
  } else {
    for (YEAR in 1:59){
      Boat_Days_Scen[,,YEAR] <- Effort_Scen[[S]][,,YEAR] / Spatial_Qs[[2]][,YEAR]
    }
  }
  Boat_Days_Scen[is.nan(Boat_Days_Scen)] <- 0
  
  Boat_Days[[S]] <- Boat_Days_Scen
}
colSums(Boat_Days[[2]][,,2])
Effort_Scen[[1]][1000:1834,,50]
Effort_Scen[[3]][1000:1834,,50]
#* Sum up effort for every year in each of the scenarios

Boat_Days_sum <- NULL

for(S in 1:4){
  temp <- Boat_Days[[S]]
  temp2 <- colSums(temp, dim=2) %>% 
    as.data.frame() %>% 
    mutate(Scenario = Names[S])
  
  Boat_Days_sum <- rbind(Boat_Days_sum, temp2)
}

Boat_Days_sum <- Boat_Days_sum %>% 
  rename(Effort = ".") %>% 
  mutate(Mortality = Effort*(0.00001)) %>% 
  mutate(Finite = 1-exp(-Mortality))
#* Spatial Plot of Peak Effort ####

## Extract peak effort - 
Normal_Effort <- Boat_Days[[2]]
peak <- rowSums(Normal_Effort[,,40])
peak <- peak[as.numeric(water_WHA$ID)]
#peak <- round(peak, digits=1)

water_WHA$Effort <- peak

water_WHA_2 <- water_WHA %>% 
  #filter(as.numeric(cell_area)>1000000) %>% 
  mutate(cell_area = cell_area/1000000) %>% 
  mutate(cell_effect = as.numeric(cell_area/sum(cell_area))) %>% 
  mutate(Effort = as.numeric(Effort/cell_area))

map <- ggplot()+
  geom_sf(data=water_WHA_2, aes(fill=Effort), color = NA, lwd=0)+
  scale_fill_carto_c(bquote(Fishing~effort~(Boat~days~km^-1)), palette="BluYl", direction=-1)+
  #annotate("text", x = 113.45, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
  geom_sf(data=ausc)+
  geom_sf(data=BR)+
  geom_richtext(data = BR, x=c(114.1730+0.14, 114.1400+0.15, 113.9784-0.165, 113.7665+0.165), 
                y=c(-21.83106, -21.95587, -21.91276, -23.15521), label = c("Bundegi", "Exmouth", "Tantabiddi", "Coral Bay"), size=2.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
map

setwd(fig_dir)
a4.width <- 160
ggsave(map, filename="ningaloo_spatial_Effort_Plot.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )



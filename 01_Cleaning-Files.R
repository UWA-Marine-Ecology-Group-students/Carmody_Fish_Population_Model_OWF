library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(stringr)
library(forcats)
library(RColorBrewer)
library(geosphere)
library(sfnetworks)
library(TSP)

rm(list = ls())

## Functions
# This returns the centre of the ploygon, but if it's on land it will create a new centroid

st_centroid_within_poly <- function (poly) { #returns true centroid if inside polygon otherwise makes a centroid inside the polygon
  
  # check if centroid is in polygon
  centroid <- poly %>% st_centroid() 
  in_poly <- st_within(centroid, poly, sparse = F)[[1]] 
  
  # if it is, return that centroid
  if (in_poly) return(centroid) 
  
  # if not, calculate a point on the surface and return that
  centroid_in_poly <- st_point_on_surface(poly) 
  return(centroid_in_poly)
}


#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
m_dir <- paste(working.dir, "Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")

model.name <- "ningaloo"


#### LOAD FILES ####

## Data
#* Spatial Data ####
setwd(sp_dir)

# Map of WA Coastline
wa_map <- st_read("WACoastline.shp")%>%
  st_transform(4283)%>%
  st_make_valid() %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5)
plot(wa_map$geometry, col="#eeeeeeff")

# Map of original State Marine Park
old_MP <- st_read("1987_NMP_boundary.shp") %>% 
  st_transform(4283)%>%
  st_make_valid %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) %>% 
  dplyr::filter(DESCRIPTOR == "sanctuary") %>% # Only pick the areas that are sanctuaries 
  mutate(Year.Sanct = 1987)
plot(old_MP$geometry)

old_MP <- old_MP %>%  # Make sure that the columns are the same as the other shape files to rbind them
  dplyr::select(Year.Sanct, NAME, geometry) %>% 
  rename(Name = "NAME") %>% 
  mutate(Name = paste0("Old", sep=" ", Name))

# Map of Current State Marine Parks
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
plot(NTZ)

# Map of Commonwealth Marine Park
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
NTZ <- rbind(NTZ, AMP_NTZ, old_MP) # Put all the different NTZs together into one object
NTZ <- st_make_valid(NTZ) 

plot(NTZ$geometry) # Check

#### SETTING UP SMALLER MODEL ####
small_model <- st_read("mangrove.gpkg") %>% 
  st_transform(4283)%>%
  st_make_valid(.)
plot(small_model$geometry)

NTZ <- NTZ %>% 
  filter(Name %in% c("Mangrove Sanctuary Zone")| Name %in% c("Old Mangrove"))
plot(NTZ$geometry, add=T)


#* Habitat Files ####
setwd(sp_dir)
habitat.types <- c("Reef", "Lagoon", "RockReef", "Pelagic") # These need to be the same as the habitat types in the names of the shapefiles you want to load, if it's all in one then ignore the loop and load the file

habitat <- data.frame()

for (i in 1:4){
  
  hab.file <- st_read(paste0(habitat.types[i], "Habitat.gpkg", sep="")) %>% 
    st_transform(4283)%>%
    st_make_valid%>%
    st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5)%>% # Has to match the area you are modelling
    mutate(type = habitat.types[i]) %>% 
    dplyr::select(type, geom)
  
  habitat <- rbind(habitat, hab.file)
  
}
plot(habitat$geom)

habitat.union <- habitat %>%  
  st_union(.) %>% 
  st_make_valid()

#### SPATIAL DATA ####

## Create extent of area you want to cover 
# site_map <- st_crop(wa_map, xmin=113.8826, xmax=113.9967, ymin=-22.0725, ymax=-21.84793) # Change the crop dimensions to whatever you want
# plot(site_map$geometry, col="#8ec3ca")

# Make grid cells for fish to live in
grd <- st_make_grid(habitat.union, cellsize=0.05, square=FALSE) %>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5)#make sure extent of grid is the same as the polygon
plot(grd)

water <- st_intersection(grd, habitat.union) %>% 
  st_make_valid()
plot(water, border="#8ec3ca")

# Make a smaller grid for the cells that are closer to the shore
GrdSmall <- st_make_grid(habitat.union, cellsize=0.025, square=FALSE) %>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5)%>%  #make sure extent of grid is the same as the polygon
  st_make_valid()
plot(GrdSmall)

HabitatSmall <- habitat %>% 
  filter(type!="Pelagic") %>% # Remove any habitats you don't want the small grids over 
  st_union(.) %>% 
  st_make_valid()
plot(HabitatSmall) #check

SmallGrd <- st_intersection(GrdSmall, HabitatSmall)
plot(SmallGrd, add = T) #check that they line up

# Merge the small grid with the grid with larger cells
BigGrd <- st_difference(water, HabitatSmall)
plot(BigGrd)
plot(SmallGrd, add=T)

water <- append(BigGrd, SmallGrd) %>% 
  st_make_valid(water)
plot(water) # check that it looks right

water_tbl <- water %>% 
  st_as_sf() %>% 
  mutate(cell_area=st_area(.)) %>% 
  filter(as.numeric(cell_area)>1) # get rid of any tiny cells that are less than 1m^2

water <- water_tbl %>% 
  dplyr::select(x)

#### SAVE FILES FOR NEXT STEP ####
setwd(data_dir)
saveRDS(NTZ, file=paste0(model.name, sep="_", "full_NTZ"))

setwd(sp_dir)
saveRDS(water, file=paste0(model.name, sep="_", "water"))


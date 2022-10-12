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

data_dir <- paste(working.dir, "Data", sep="/") # For raw data
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/") # For spatial dataframes that we're still working on
sg_dir <- paste(working.dir, "Staging", sep="/") # For finished dataframes that we aren't going to change any more

model.name <- "ningaloo"


#### LOAD FILES ####

## Data
setwd(data_dir)
NTZ <- readRDS(paste0(model.name, sep="_", "full_NTZ")) # Read in your cleaned NTZ file, whether this was done in the cleaning script or hasn't really changed over time 

setwd(sp_dir)
water <- readRDS(paste0(model.name, sep="_", "water"))
plot(water)

NTZ.timing <- c("1960", "1987", "2005", "2017") # Put in when NTZs have changed over time, if they have

# Buffer to restrict the amount of water
buffer_85m <- st_read("85m_Buffer_Layer.shp")%>% # To about 150m
  st_transform(4283)%>%
  st_make_valid

# Create a grid layer that is just for NTZs so we can adjust fishing mortality later
NTZarea <- st_intersection(NTZ, water) # This messes up the attributes and the spatial stuff so we have to fix this later on
NTZarea <- st_make_valid(NTZarea)
plot(NTZarea$geometry)

# Change the water layer so that it excludes the NTZs because we want to be able to differentiate them
NTZ_union <- st_union(NTZ)
plot(NTZ)

Fished_area <- st_difference(water, NTZ_union) %>% 
  st_make_valid()
plot(Fished_area)

# turn areas into data frames for easier use
Fished_area <- st_sf(Fished_area) 

for (i in 1:length(NTZ.timing)){
  
  colname <- paste0("Fished", sep="_", NTZ.timing[i])
  
  Fished_area <- Fished_area %>% 
    mutate(!!colname := "Y")
}

Fished_area <- Fished_area %>% 
  mutate(Year.Sanct = NA) %>% 
  mutate(Name = NA)

names(NTZarea)[names(NTZarea) == 'geometry'] <- 'Spatial'
st_geometry(NTZarea) <- "Spatial"

names(Fished_area)[names(Fished_area) == 'x'] <- 'Spatial'
st_geometry(Fished_area) <- "Spatial"

plot(Fished_area$Spatial, border="blue")

#### SORTING OUT NO TAKE CELLS IF THERE ARE ISSUES ####
NoTake <- st_join(NTZarea, NTZ) # Somewhere along the lines all of the comments got separated from their IDs/Comments so now we need to put things back together to get the right labels with the right cells
NoTake <- NoTake %>% 
  dplyr::select(Name.y, Year.Sanct.y, Spatial) %>% 
  rename(Name = "Name.y",
         Year.Sanct = "Year.Sanct.y") %>% 
  mutate(Old.Sanct = ifelse(Year.Sanct==1987, "Old", "New")) %>% 
  distinct(Old.Sanct, Spatial, .keep_all=T) %>% 
  dplyr::select(-Old.Sanct) %>% 
  mutate(ID=row_number())

## Fixing Cloates sanctuary zone
setwd(sp_dir)
AMP_NTZ <- st_read("NTZ and Fished areas for status.shp") %>% 
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) %>% 
  filter(ResName == "Ningaloo"|COMMENTS == "Cloates Sanctuary Zone") %>% # Select just Cloates
  mutate(Year.Sanct = ifelse(ResName %in% c("Ningaloo"), 2017, 2005)) %>% 
  rename(Name = "COMMENTS") %>% 
  mutate(Name = ifelse(is.na(Name), "Comm Cloates", Name)) %>% 
  dplyr::select(Name, Year.Sanct, geometry) 

NTZ_state <- NoTake %>% 
  filter(Year.Sanct==2005)
plot(NTZ_state$Spatial)

NTZ_comm <- NoTake %>% 
  st_intersects(., AMP_NTZ[1, ]) %>% 
  as.data.frame()
NTZ_comm <- NoTake[c(as.numeric(NTZ_comm$row.id)), ]
plot(NTZ_comm$Spatial)

Cloates_overlap <- st_touches(NTZ_comm, NTZ_state)
Cloates_overlap <- as.data.frame(Cloates_overlap) # The row ID from NTZ_comm are the cells that need to be reassigned back to the state NTZ, not commonwealth
Cloates_overlap <- NTZ_comm[as.numeric(c(Cloates_overlap$row.id)),]
plot(Cloates_overlap$Spatial)

NoTake <- NoTake %>% 
  mutate(Name = ifelse(ID %in% Cloates_overlap$ID, "Cloates Sanctuary Zone", Name)) %>% 
  mutate(Year.Sanct = ifelse(ID %in% Cloates_overlap$ID, 2005, Year.Sanct))

check <- NoTake%>% 
  filter(Name %in% ("Cloates Sanctuary Zone"))
plot(check$Spatial)

## Fixing 1987 sanctuary zones
setwd(sp_dir)

old_MP <- st_read("1987_NMP_boundary.shp") %>% 
  st_transform(4283)%>%
  st_make_valid %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) %>% 
  dplyr::filter(DESCRIPTOR == "sanctuary") %>% # Only pick the areas that are sanctuaries 
  mutate(Year.Sanct = 1987) %>% # Make sure that the columns are the same as the other shape files to rbind them
  dplyr::select(Year.Sanct, NAME, geometry) %>% 
  rename(Name = "NAME") %>% 
  mutate(Name = paste0("Old", sep=" ", Name))

NTZ_state <- NoTake %>% 
  filter(Year.Sanct==1987)
NTZ_state_points <- st_centroid(NTZ_state)
plot(NTZ_state_points$Spatial)

NTZ_old_poly <- st_cast(st_geometry(old_MP), "LINESTRING")
NTZ_old_poly <- st_polygonize(NTZ_old_poly)
plot(NTZ_old_poly)

inside_NTZ <- st_covered_by(NTZ_state_points, NTZ_old_poly) %>% 
  as.data.frame()

old_NTZ_ID <- NTZ_state[-as.numeric(c(inside_NTZ$row.id)),]
plot(NTZ_state$Spatial)
plot(old_NTZ_ID$Spatial, add=T, col="Blue")

NoTake <- NoTake %>% 
  mutate(Year.Sanct = ifelse(ID %in% old_NTZ_ID$ID, 2005, Year.Sanct)) 
plot(NoTake$Spatial)

NTZ_old_check <- NoTake %>% 
  filter(Year.Sanct==1987) %>% 
  ggplot(.)+
  geom_sf()+
  theme_void()

NTZ_2005_check <- NoTake %>% 
  filter(Year.Sanct==2005) %>% 
  ggplot(.)+
  geom_sf()+
  theme_void()

NTZ_2017_check <- NoTake %>% 
  filter(Year.Sanct==2017) %>% 
  ggplot(.)+
  geom_sf()+
  theme_void()

## Put this back together with the fished area to create "water" again

for (i in 1:length(NTZ.timing)){
  
  colname <- paste0("Fished", sep="_", NTZ.timing[i])
  
  if(NTZ.timing[i] %in% c("1960")){
    NoTake <- NoTake %>% 
      mutate(!!colname := "Y")
  } else if(NTZ.timing[i] %in% c("1987")){
    NoTake <- NoTake %>% 
      mutate(!!colname := ifelse(Year.Sanct==1987, "N", "Y"))
  } else if(NTZ.timing[i] %in% c("2005")){
    NoTake <- NoTake %>% 
      mutate(!!colname := ifelse(Year.Sanct==1987|Year.Sanct==2005, "N", "Y"))
  } else if(NTZ.timing[i] %in% c("2017")){
    NoTake <- NoTake %>% 
      mutate(!!colname := "N")
  }
}  

NoTake <- NoTake %>% 
  dplyr::select(-ID)

water <- rbind(NoTake, Fished_area)
plot(water$Spatial)

## Check that the NTZs are where you expect them to be
ggplot(water)+
  geom_sf(aes(fill=Fished_2005))+
  theme_void()+
  scale_fill_manual(values=c("#f3c1af", "#c8dfe3"))

water <- st_make_valid(water) %>% 
  st_as_sf()

#water <- water[-299, ] # SPECIFIC FOR WHOLE NINGALOO MODEL

## Remove cells that are greater than a certain depth
buffer_85m2 <- st_crop(buffer_85m, xmin=112.5, xmax=114.65, ymin=-24, ymax=-21.1)

buffer_85m2 <- st_union(buffer_85m2) %>% 
  st_make_valid()

plot(water$Spatial)
plot(buffer_85m2, add=T)

buffered <- st_intersection(water, buffer_85m2)
plot(buffered)

water <- st_make_valid(buffered)

water <- water %>% 
  mutate(Old.Sanct = ifelse(Year.Sanct==1987, "A", "B")) %>%
  arrange(Old.Sanct, Name) %>%
  distinct(Spatial, .keep_all=T) %>%  # Remove any duplicated polygons if they happen to have cropped up
  st_make_valid() %>%
  mutate(ID=row_number()) %>%
  dplyr::select(-Old.Sanct)

check <- water %>% 
  filter(Fished_2005 %in% c("N")) %>% # Change this to check the different NTZs are where you think they should be
  ggplot(.)+
  geom_sf()+
  theme_void() 

## Remove polygons with a 0 area ##
water <- water %>% 
  mutate(cell_area = st_area(Spatial)) %>% 
  filter(as.numeric(cell_area)>1)

#### Create a list of the NTZ IDs ####
NoTakeList <- list()

NoTake87 <- water %>% 
  st_drop_geometry() %>% 
  filter(Fished_1987 %in% c("N")) %>% 
  dplyr::select(ID) 
NoTake87 <- as.numeric(NoTake87$ID)

NoTake05 <- water %>% 
  st_drop_geometry() %>% 
  filter(Fished_2005 %in% c("N")) %>% 
  dplyr::select(ID) 
NoTake05 <- as.numeric(NoTake05$ID)

NoTake17 <- water %>% 
  st_drop_geometry() %>% 
  filter(Fished_2017 %in% c("N")) %>% 
  dplyr::select(ID) 
NoTake17 <- as.numeric(NoTake17$ID)

NoTakeList[[1]] <- NoTake87
NoTakeList[[2]] <- NoTake05
NoTakeList[[3]] <- NoTake17

#### SAVE FILES FOR NEXT STEP ####
setwd(sg_dir)
saveRDS(NoTakeList, file=paste0(model.name, sep="_", "NoTakeList"))

setwd(sp_dir)
saveRDS(water, file=paste0(model.name, sep="_","water"))








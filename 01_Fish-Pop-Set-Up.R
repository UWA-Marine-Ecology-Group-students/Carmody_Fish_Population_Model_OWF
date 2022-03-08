library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(stringr)
library(forcats)
library(RColorBrewer)
library(geosphere)



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

## Create colours for the plot
cols <- brewer.pal(8, "RdBu")
levels_water <- data.frame(c("<1000", "1000-1100", "1100-1200", "1200-1300",
                             "1300-1400", "1400-1500", "1500-1600", ">1600"))
names(levels_water)[1] <- "Levels"
levels_water$Levels <- as.factor(levels_water$Levels)
names(cols) <- levels(levels_water$Levels)


#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
m_dir <- paste(working.dir, "Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")

#### LOAD FILES ####

## Data
setwd(data_dir)
boat_days <- read.csv("Boat_Days_Gascoyne.csv")

boat_days <- boat_days%>%
  mutate(NumMonth = as.numeric(NumMonth)) %>% 
  mutate(Month = as.factor(Month)) %>% 
  mutate(Gascoyne_Boat_Days = as.numeric(Gascoyne_Boat_Days)) %>% 
  mutate(Month = fct_relevel(Month, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))

## Spatial Data
setwd(sp_dir)

# Map of WA Coastline
wa_map <- st_read("WACoastline.shp")%>%
  st_transform(4283)%>%
  st_make_valid
plot(wa_map$geometry, col="#eeeeeeff")

# Map of original State Marine Park
old_MP <- st_read("1987_NMP_boundary.shp") %>% 
  st_transform(4283)%>%
  st_make_valid %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) %>% 
  dplyr::filter(DESCRIPTOR == "sanctuary")
plot(old_MP$geometry)

old_MP <- old_MP %>%  # Make sure that the columns are the same as the other shape files to rbind them
  dplyr::select(RECNUM, NAME, COMMENTS, geometry) %>% 
  rename(PA_ID = "RECNUM") %>% 
  mutate(COMMENTS = "Old") %>% 
  mutate(COMMENTS = paste(NAME, COMMENTS, sep=" ")) %>% 
  dplyr::select(-NAME)

# Map of Current State Marine Parks
MP <- st_read("WA_MPA_2018.shp")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5)

NTZ <- MP%>%
  filter(IUCN == "IA") %>% 
  filter(!COMMENTS == "Cloates Sanctuary Zone") %>% 
  dplyr::select(PA_ID,COMMENTS, geometry) %>% 
  filter(!COMMENTS %in% c("North Muiron Conservation Area","South Muiron Conservation Area",
                          "Previously PA_ID WA_42756 in terrestrial CAPAD", "Conservation Area"))
plot(NTZ)

# Map of Commonwealth Marine Park
AMP <- st_read("NTZ and Fished areas for status.shp") %>% 
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5)
plot(AMP)

AMP_NTZ <- AMP %>% 
  filter(ResName == "Ningaloo"|COMMENTS == "Cloates Sanctuary Zone") %>% 
  dplyr::select(PA_ID, COMMENTS, geometry) %>% 
  mutate(COMMENTS = ifelse(is.na(COMMENTS), "Comm Cloates", COMMENTS))
plot(AMP_NTZ$geometry)

NTZ <- rbind(NTZ, AMP_NTZ, old_MP) 
NTZ <- st_make_valid(NTZ)

plot(NTZ$geometry)

# Habitat Layers
reef <- st_read("ReefHabitat.gpkg")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.3, ymin=-24, ymax=-21)%>% 
  dplyr::select(geom)

reef$type <- "reef"

lagoon <- st_read("LagoonHabitat.gpkg")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.3, ymin=-24, ymax=-21) %>% 
  dplyr::select(geom)

lagoon$type <- "lagoon"

rocky <- st_read("RockReefHabitat.gpkg")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.3, ymin=-24, ymax=-21)%>% 
  dplyr::select(geom)

rocky$type <- "rocky"

pelagic <- st_read("PelagicHabitat.gpkg")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.3, ymin=-24, ymax=-21)

pelagic$type <- "pelagic"

#### SPATIAL DATA ####

## Create extent of area you want to cover 

ningaloo_map <- st_crop(wa_map, xmin=100, xmax=116, ymin=-24, ymax=-21.1) #NEED TO CHANGE THIS TO BE WHOLE GULF DOWN TO CORAL BAY
plot(ningaloo_map$geometry, col="#8ec3ca")

# Make grid cells for fish to live in
grd <- st_make_grid(ningaloo_map, cellsize=0.05, square=FALSE) %>%
  st_crop(xmin=100, xmax=114.65, ymin=-24, ymax=-21.1) #make sure extent of grid is the same as the polygon
plot(grd, add=TRUE)

water <- st_difference(grd, ningaloo_map)
plot(water, border="#8ec3ca")

# Adjust the sizes of the grid cells so that the ones closer to the shore are smaller 
GrdSmall <- st_make_grid(ningaloo_map, cellsize=0.025, square=FALSE) %>%
  st_crop(xmin=100, xmax=114.65, ymin=-24, ymax=-20) #make sure extent of grid is the same as the polygon
plot(GrdSmall)

HabitatSmall <- rbind(reef, rocky, lagoon)
HabitatSmall <- st_union(HabitatSmall) %>% 
  st_make_valid()
plot(HabitatSmall)

SmallGrd <- st_intersection(GrdSmall, HabitatSmall)
plot(SmallGrd, add = T)

# Merge the small grid with the grid with larger cells
BigGrd <- st_difference(water, HabitatSmall)
plot(BigGrd)
plot(SmallGrd, add=T)

BigGrd <- st_make_valid(BigGrd)

water <- append(BigGrd, SmallGrd)
plot(water)

# Create a grid layer that is just for NTZs so we can adjust fishing mortality later
NTZarea <- st_intersection(NTZ, water) # Gives you just the hexagons in each NTZ with the comments
NTZarea <- st_make_valid(NTZarea)
plot(NTZarea$geometry)

# Change the water layer so that it excludes the NTZs because we want to be able to differentiate them
NTZ <- st_union(NTZ)
plot(NTZ)

Fished_area <- st_difference(water, NTZ) 
plot(Fished_area)

# turn areas into data frames for easier use
Fished_area <- st_sf(Fished_area) %>% 
  mutate(Fished = "Y") %>% 
  rename(Spatial = "Fished_area") %>% 
  mutate(PA_ID = NA) %>% 
  mutate(COMMENTS = NA)

NTZ_area <- st_sf(NTZarea) %>% 
  mutate(Fished = "N")

names(NTZ_area)[names(NTZ_area) == 'geometry'] <- 'Spatial'
st_geometry(NTZ_area) <- "Spatial"

water <- rbind(Fished_area, NTZ_area)

water <- water%>%
  mutate(ID=row_number())%>%
  mutate(ID=factor(ID)) %>% 
  dplyr::select(-PA_ID)

## Check that the NTZs are where you expect them to be
ggplot(water)+
  geom_sf(aes(fill=Fished))+
  theme_void()+
  scale_fill_manual(values=c("#f3c1af", "#c8dfe3"))

water <- st_make_valid(water)

# Get centroids for the grid cells 
centroids <- st_centroid_within_poly(water)
plot(centroids, cex=0.3) 

points <- as.data.frame(st_coordinates(centroids))%>% #The points start at the bottom left and then work their way their way right
  mutate(ID=row_number())

NCELL <- nrow(points) #Set the number of cells in the model

## Add habitat data and assign it to the grid cells e.g. % cover 
## Can then put a coefficient on the different types of habitat to change how the fish move (i.e. they move towards reef)

# Join all the habitat together
reef <- reef%>%
  dplyr::select(geom, type)
lagoon <- lagoon%>%
  dplyr::select(geom, type)
rocky <- rocky%>%
  dplyr::select(geom, type)
pelagic <- pelagic%>%
  dplyr::select(geom, type)

habitat <- rbind(reef, lagoon, rocky, pelagic)
plot(habitat$geom)
plot(water, add=T)

# Calculate the percentage of each habitat type in each of the different cells 
intersection <- st_intersection(habitat, water) #This gives you all the individual areas where there is overlap, the ID of the original grid cell and the habitat in each of the new little cells
plot(intersection$geom)

intersection <- st_make_valid(intersection)

intersection$area <- st_area(intersection) #gives you the area each of the small cells covers 
water$area <- st_area(water) #gives you the area of each of the cells in the grid

hab_perc <- merge(st_drop_geometry(intersection), st_drop_geometry(water), by.x="ID", by.y="ID")%>%
  mutate(area.x = as.numeric(area.x))%>%
  mutate(area.y = as.numeric(area.y))%>%
  mutate(perc_habitat = ((area.x/area.y)*100)) # This tells you how much of each water grid cell is made up of the different habitat types 

# Create separate data frames for each habitat type and fill in for where cells don't have a certain habitat type

pelagic_perc <- hab_perc%>%
  filter(type=="pelagic")%>%
  dplyr::select(ID, type, perc_habitat)%>%
  mutate(ID = as.numeric(ID))%>%
  complete(ID = 1:NCELL, fill = list(perc_habitat=0))%>%
  mutate(type = replace_na("pelagic"))

reef_perc <- hab_perc%>%
  filter(type=="reef")%>%
  dplyr::select(ID, type, perc_habitat)%>%
  mutate(ID = as.numeric(ID))%>%
  complete(ID = 1:NCELL, fill = list(perc_habitat=0))%>%
  mutate(type = replace_na("reef"))

rocky_perc <- hab_perc%>%
  filter(type=="rocky")%>%
  dplyr::select(ID, type, perc_habitat)%>%
  mutate(ID = as.numeric(ID))%>%
  complete(ID = 1:NCELL, fill = list(perc_habitat=0))%>%
  mutate(type = replace_na("rocky"))

lagoon_perc <- hab_perc%>%
  filter(type=="lagoon")%>%
  dplyr::select(ID, type, perc_habitat)%>%
  mutate(ID = as.numeric(ID))%>%
  complete(ID = 1:NCELL, fill = list(perc_habitat=0))%>%
  mutate(type = replace_na("lagoon"))

#### CREATE MATRIX OF CONNVECTIVITY FOR FISH MOVEMENT ####
## Calculate a distance matrix between the cells 
# To make them go round the land, create visibility graph (grass function), then distance between nodes, then pick the
# shortest path
dist_matrix <- as.matrix(dist(points)) #This is just the linear distance  between the centroids of the grid cells

## Calculate the probability a fish moves to this site in a given time step using the swimming speed we set earlier 
# This creates a dispersal kernel based on the negative exponential distribution

pDist <- matrix(NA, ncol=NCELL, nrow=NCELL)

for(r in 1:NCELL){
  for(c in 1:NCELL){
    p <- exp(-dist_matrix[r,c]/SwimSpeed)
    pDist[r,c] <- p
  }
}

## Calculate the difference in habitat types between each of the cells i.e. will there be an increase in reef % if you go from site 1 to site 2

pPelagic <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    p <- as.numeric((pelagic_perc[c,3]) - (pelagic_perc[r,3]))
    pPelagic[r,c] <- p
  }
}

pReef <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    p <- as.numeric((reef_perc[c,3]) - (reef_perc[r,3]))
    pReef[r,c] <- p
  }
}

pLagoon <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    p <- as.numeric((lagoon_perc[c,3]) - (lagoon_perc[r,3]))
    pLagoon[r,c] <- p
  }
}

pRocky<- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    p <- as.numeric((rocky_perc[c,3]) - (rocky_perc[r,3]))
    pRocky[r,c] <- p
  }
}

#### CREATE PROBABILITY OF MOVEMENT USING UTILITY FUNCTION ####
# First determine the utility of each of the sites 
# This is very sensitive to changes in the values particularly for reef 
# PROBABLY ALSO NEED TO PUT DEPTH IN HERE

a = 1
b = 1
c = 1
d = 1
e = 1

Vj <- (a*pDist) + (b*pReef) + (c*pLagoon) + (d*pRocky) + (e*pPelagic)

# Calculate the summed utility across the rows 
rowU <- matrix(NA, ncol=1, nrow=NCELL)
cellU <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    U <- exp(Vj[r,c])
    cellU[r,c] <- U
  }
}

rowU <- as.data.frame(rowSums(cellU))

# Calculate the probability that the fish will move to this site

Pj <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    Pj[r,c] <- (exp(Vj[r,c]))/rowU[r,1]
  }
}
rowSums(Pj)

#### RECRUITMENT MATRIX ####
## Want the recruits to be in the lagoons and then move out from there 
dispersal <- lagoon_perc %>% 
  dplyr::select(perc_habitat)
dispersal$area <- as.vector(water$area)

temp <- array(0, dim=c(NCELL, 1))
recruitment <- array(0, dim=c(NCELL, 1))

for (CELL in 1:NCELL){
  
  for (cell in 1:NCELL){
  temp[cell, 1] <- as.numeric((dispersal[cell,2]/dispersal[cell,1]))
  }
  
  temp[which(!is.finite(temp))] <- 0
  summed <- sum(temp)
  
  recruitment[CELL,1] <- as.numeric(temp[CELL,1]/summed)
}

#### RECRUIT MOVEMENT ####
## Want the recruits to stay in the lagoon until they mature and move to the reef

a = 1
b = 0.8
c = 3
d = 0.6
e = 0.01

Recj <- (a*pDist) + (b*pReef) + (c*pLagoon) + (d*pRocky) + (e*pPelagic)

# Calculate the summed utility across the rows 
rowU <- matrix(NA, ncol=1, nrow=NCELL)
cellU <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    U <- exp(Recj[r,c])
    cellU[r,c] <- U
  }
}

rowU <- as.data.frame(rowSums(cellU))

# Calculate the probability that the fish will move to this site

ProbRec <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    ProbRec[r,c] <- (exp(Recj[r,c]))/rowU[r,1]
  }
}
rowSums(ProbRec)

#### SAVE FILES ####
setwd(sg_dir)
saveRDS(Pj, file="movement")
saveRDS(water, file="water")
saveRDS(ProbRec, file="juvmove")

#### SST FOR RECRUITMENT ####













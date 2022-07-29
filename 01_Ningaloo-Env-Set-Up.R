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

SwimSpeed <- 1

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

# Buffer to restrict the amount of water
buffer_85m <- st_read("85m_Buffer_Layer.shp")%>%
  st_transform(4283)%>%
  st_make_valid

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

water <- water[-299, ]

## Remove cells that are greater than a certain depth
buffer_85m2 <- st_crop(buffer_85m, xmin=112.5, xmax=114.65, ymin=-24, ymax=-21.1)

buffer_85m2 <- st_union(buffer_85m2) %>% 
  st_make_valid()

plot(water$Spatial)
plot(buffer_85m2, add=T)

buffered <- st_intersection(water, buffer_85m2)
plot(buffered)

buffered <- st_make_valid(buffered)

water <- buffered %>% 
  mutate(ID=row_number())%>%
  mutate(ID=factor(ID))

# Get centroids for the grid cells 
centroids <- st_centroid_within_poly(water)
plot(centroids, cex=0.3) 

points <- as.data.frame(st_coordinates(centroids))%>% #The points start at the bottom left and then work their way their way right
  mutate(ID=row_number())

NCELL <- nrow(points) #Set the number of cells in the model


#### CREATE THE DISTANCE MATRIX ####
## This method allows us to create a distance matrix that takes into account the land and forces the fish to swim around

## Convert the points in the centroids of the ploygon to a spatial points file
points$ID <- as.integer(points$ID)
points2 <- st_as_sf(points, coords = c("X", "Y"))
points2 <- st_cast(st_geometry(points2), "POINT")

## Figure out which six points are the closest - won't all have six close ones but can start from there
dist.mat <- st_distance(points2) # Great Circle distance since in lat/lon
# Number within 1.5km: Subtract 1 to exclude the point itself
num.5km <- apply(dist.mat, 1, function(x) {
  sum(x < 0.05) - 1
})

nn.dist <- apply(dist.mat, 1, function(x) {
  return(sort(x, partial = 2)[2])
})

## Get the IDS for the cells that are neighbours
closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[2] })
second.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[3] })
third.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[4] })
fourth.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[5] })
fifth.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[6] })
sixth.closest <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[7] })

neighbours <- water %>% 
  dplyr::select(ID) %>% 
  mutate(`closest` = closest) %>% 
  mutate(second = second.closest) %>% 
  mutate(third = third.closest) %>% 
  mutate(fourth = fourth.closest) %>% 
  mutate(fifth = fifth.closest) %>% 
  mutate(sixth = sixth.closest)

## Give the neighbouring points geometry based on the original set of points
closest <- as.data.frame(closest) %>%
  rename(ID = "closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y")) 
closest <- st_cast(st_geometry(closest), "POINT") 

second.closest <- as.data.frame(second.closest) %>%
  rename(ID = "second.closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y"))
second.closest <- st_cast(st_geometry(second.closest), "POINT")

third.closest <- as.data.frame(third.closest) %>%
  rename(ID = "third.closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y"))
third.closest <- st_cast(st_geometry(third.closest), "POINT")

fourth.closest <- as.data.frame(fourth.closest) %>%
  rename(ID = "fourth.closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y"))
fourth.closest <- st_cast(st_geometry(fourth.closest), "POINT")

fifth.closest <- as.data.frame(fifth.closest) %>%
  rename(ID = "fifth.closest") %>% 
  inner_join(., points, by="ID") %>%  
  st_as_sf(., coords = c("X", "Y"))
fifth.closest <- st_cast(st_geometry(fifth.closest), "POINT")

sixth.closest <- as.data.frame(sixth.closest) %>%
  rename(ID = "sixth.closest") %>% 
  inner_join(., points, by="ID") %>% 
  st_as_sf(., coords = c("X", "Y"))
sixth.closest <- st_cast(st_geometry(sixth.closest), "POINT")

### CONNECT THE POINTS TO THEIR NEIGHBOURS TO FORM A NETWORK ###
n <- nrow(points)

## Form linestrings and then multilinestrings
# Closest
linestrings.closest <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points2[x], closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.closest <- st_multilinestring(do.call("rbind", linestrings.closest))

# Second closest
linestrings.second <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points2[x], second.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.second <- st_multilinestring(do.call("rbind", linestrings.second))

# Third closest
linestrings.third <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points2[x], third.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.third<- st_multilinestring(do.call("rbind", linestrings.third))

# Fourth closest
linestrings.fourth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points2[x], fourth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.fourth <- st_multilinestring(do.call("rbind", linestrings.fourth))

# Fifth closest
linestrings.fifth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points2[x], fifth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.fifth <- st_multilinestring(do.call("rbind", linestrings.fifth))

# Sixth closest
linestrings.sixth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points2[x], sixth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.sixth <- st_multilinestring(do.call("rbind", linestrings.sixth))

connected <- st_combine(c(multilinestring.closest, multilinestring.second, multilinestring.third, multilinestring.fourth, multilinestring.fifth, multilinestring.sixth))
connected <- st_cast(connected, "LINESTRING") # Needs to be a line string rather than multiline for the next step

### SET UP THE SF NETWORK AND CREATE A DISTANCE MATRIX ###
network <- as_sfnetwork(connected, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

## Calculate the distances from each point to every other point on the network
net <- activate(network, "nodes")
cost_matrix <- st_network_cost(net)
dim(cost_matrix) # Check that the dimensions match up to how many points you think you should have in the network

dist_matrix <- dist(points)

#### ADD HABITAT TO THE ENVIRONMENT ####

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
intersection <- st_intersection(habitat, water) #This gives you all the individual areas where there is overlap, the ID of the original grid cell and the habitat in each of the new little cells=plot(intersection$geom)

intersection <- st_make_valid(intersection)

intersection$hab_area <- st_area(intersection) #gives you the area each of the small cells covers 
water$cell_area <- st_area(water) #gives you the area of each of the cells in the grid

hab_perc <- merge(st_drop_geometry(intersection), st_drop_geometry(water), by.x="ID", by.y="ID")%>%
  mutate(hab_area = as.numeric(hab_area))%>%
  mutate(cell_area = as.numeric(cell_area))%>%
  mutate(perc_habitat = ((hab_area/cell_area)*100)) %>%  # This tells you how much of each water grid cell is made up of the different habitat types 
  mutate(perc_habitat = ifelse(perc_habitat>100, 100, perc_habitat)) # The intersection has given a couple of places where the % is just over 100 so just round these back to 100

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


#### SAVE FILES FOR USE IN NEXT STEP ####
setwd(sp_dir)

saveRDS(dist_matrix, file="dist_matrix")
saveRDS(pelagic_perc, file="pelagic_perc")
saveRDS(reef_perc, file="reef_perc")
saveRDS(lagoon_perc, file="lagoon_perc")
saveRDS(rocky_perc, file="rocky_perc")













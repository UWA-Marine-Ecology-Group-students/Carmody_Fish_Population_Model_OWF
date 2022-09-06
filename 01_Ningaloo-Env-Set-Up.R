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
#* Spatial Data ####
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
  dplyr::filter(DESCRIPTOR == "sanctuary") # Only pick the areas that are sanctuaries
plot(old_MP$geometry)

old_MP <- old_MP %>%  # Make sure that the columns are the same as the other shape files to rbind them
  dplyr::select(RECNUM, NAME, COMMENTS, geometry) %>% 
  mutate(COMMENTS_2005 = NA,
         COMMENTS_2017 = NA) %>% 
  rename(PA_ID = "RECNUM") %>% 
  mutate(COMMENTS_OLD = "Old") %>% 
  mutate(COMMENTS_OLD = paste(NAME, COMMENTS_OLD, sep=" ")) %>% 
  dplyr::select(-NAME, -COMMENTS)

# Map of Current State Marine Parks
MP <- st_read("WA_MPA_2018.shp")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5)

NTZ <- MP%>%
  filter(IUCN == "IA") %>% 
  filter(!COMMENTS == "Cloates Sanctuary Zone") %>% # Cloates is included in the Australian marine parks shapefile because I couldn't get them to line up nicely using the State MP shapefile
  mutate(COMMENTS_OLD = NA, # Add the correct columns so they al bind
         COMMENTS_2005 = NA,
         COMMENTS_2017 = NA) %>% 
  mutate(COMMENTS_2005 = paste(COMMENTS, sep=" ")) %>% 
  dplyr::select(PA_ID, COMMENTS_OLD, COMMENTS_2005, COMMENTS_2017, geometry) %>% 
  filter(!COMMENTS_2005 %in% c("Previously PA_ID WA_42756 in terrestrial CAPAD", "Conservation Area")) # Remove everything that isn't a sanctuary zone
plot(NTZ)

# Map of Commonwealth Marine Park
AMP <- st_read("NTZ and Fished areas for status.shp") %>% 
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5)
plot(AMP)

#
AMP_NTZ <- AMP %>% 
  filter(ResName == "Ningaloo"|COMMENTS == "Cloates Sanctuary Zone") %>% # Select just Cloates
  mutate(COMMENTS_OLD = NA,
         COMMENTS_2005 = NA,
         COMMENTS_2017 = NA) %>% 
  mutate(COMMENTS_2017 = paste(COMMENTS, sep=" ")) %>% 
  dplyr::select(PA_ID, COMMENTS_OLD, COMMENTS_2005, COMMENTS_2017, geometry) %>% 
  mutate(COMMENTS_2017 = ifelse(COMMENTS_2017 %in% c("NA"), "Comm Cloates", COMMENTS_2017))
plot(AMP_NTZ$geometry)

NTZ <- rbind(NTZ, AMP_NTZ, old_MP) # Put all the different NTZs together into one object
NTZ <- st_make_valid(NTZ) %>% 
  mutate(COMMENTS_2005 = ifelse(COMMENTS_2017 %in% c("Cloates Sanctuary Zone"), "Cloates Sanctuary Zone", COMMENTS_2005)) %>% 
  mutate(COMMENTS_2017 = ifelse(is.na(COMMENTS_2017), paste(COMMENTS_2005), COMMENTS_2017)) %>% # Make sure the comments carry over so that the "Old" NTZ polygons are labelled as old in 2005/2017 because they overlap with the expanded NTZ in 2005
  mutate(COMMENTS_2005 = ifelse(is.na(COMMENTS_2005), paste(COMMENTS_OLD), COMMENTS_2005)) %>% 
  mutate(COMMENTS_2017 = ifelse(COMMENTS_2017 %in% c("NA"), paste(COMMENTS_2005), COMMENTS_2017))

plot(NTZ$geometry) # Check


# Buffer to restrict the amount of water
buffer_85m <- st_read("85m_Buffer_Layer.shp")%>% # To about 150m
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

ningaloo_map <- st_crop(wa_map, xmin=100, xmax=116, ymin=-24, ymax=-21.1) 
plot(ningaloo_map$geometry, col="#8ec3ca")

# Make grid cells for fish to live in
grd <- st_make_grid(ningaloo_map, cellsize=0.05, square=FALSE) %>%
  st_crop(xmin=100, xmax=114.65, ymin=-24, ymax=-21.1) #make sure extent of grid is the same as the polygon
plot(grd, add=TRUE)

water <- st_difference(grd, ningaloo_map)
plot(water, border="#8ec3ca")

# Make a smaller grid for the cells that are closer to the shore
GrdSmall <- st_make_grid(ningaloo_map, cellsize=0.025, square=FALSE) %>%
  st_crop(xmin=100, xmax=114.65, ymin=-24, ymax=-20) #make sure extent of grid is the same as the polygon
plot(GrdSmall)

HabitatSmall <- rbind(reef, rocky, lagoon)
HabitatSmall <- st_union(HabitatSmall) %>% 
  st_make_valid()
plot(HabitatSmall) #check

SmallGrd <- st_intersection(GrdSmall, HabitatSmall)
plot(SmallGrd, add = T) #check that they line up

# Merge the small grid with the grid with larger cells
BigGrd <- st_difference(water, HabitatSmall)
plot(BigGrd)
plot(SmallGrd, add=T)

BigGrd <- st_make_valid(BigGrd)

water <- append(BigGrd, SmallGrd)
plot(water) # check that it looks right

# Create a grid layer that is just for NTZs so we can adjust fishing mortality later
NTZarea <- st_intersection(NTZ, water) # This messes up the attributes and the spatial stuff so we have to fix this later on
NTZarea <- st_make_valid(NTZarea)
plot(NTZarea$geometry)

# Change the water layer so that it excludes the NTZs because we want to be able to differentiate them
NTZ_union <- st_union(NTZ)
plot(NTZ)

Fished_area <- st_difference(water, NTZ_union) %>% 
  st_make_valid()

# turn areas into data frames for easier use
Fished_area <- st_sf(Fished_area) %>% 
  mutate(Fished_60 = "Y") %>% 
  mutate(Fished_87 = "Y") %>% 
  mutate(Fished_05 = "Y") %>% 
  mutate(Fished_17 = "Y") %>% 
  mutate(COMMENTS_OLD = NA,
         COMMENTS_2005 = NA,
         COMMENTS_2017 = NA) 

NTZ_area <- st_sf(NTZarea) %>% 
  mutate(Fished_60 = "Y") %>% 
  mutate(Fished_87 = ifelse(is.na(COMMENTS_OLD), "Y", "N")) %>% 
  mutate(Fished_05 = ifelse(COMMENTS_2005 %in% c("NA"), "Y", "N")) %>% 
  mutate(Fished_17 = ifelse(COMMENTS_2017 %in% c("NA"), "Y", "N")) %>% 
  dplyr::select(-PA_ID)


names(NTZ_area)[names(NTZ_area) == 'geometry'] <- 'Spatial'
st_geometry(NTZ_area) <- "Spatial"

names(Fished_area)[names(Fished_area) == 'Fished_area'] <- 'Spatial'
st_geometry(Fished_area) <- "Spatial"


#### SORTING OUT NO TAKE CELLS ####
NoTake <- st_join(NTZ_area, NTZ) # Somewhere along the lines all of the comments got separated from their IDs/Comments so now we need to put things back together to get the right labels with the right cells
NoTake <- NoTake %>% 
  dplyr::select(COMMENTS_OLD.y, COMMENTS_2005.y, COMMENTS_2017.y, Spatial) %>% 
  rename(COMMENTS_OLD = "COMMENTS_OLD.y",
         COMMENTS_2005 = "COMMENTS_2005.y",
         COMMENTS_2017 = "COMMENTS_2017.y") %>% 
  distinct(COMMENTS_OLD, Spatial, .keep_all=T) %>% 
  mutate(ID=row_number())
plot(NoTake$Spatial)

## Fixing Cloates sanctuary zone
NTZ_state <- NoTake %>% 
  filter(!COMMENTS_2005 %in% c("NA"))

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
  mutate(COMMENTS_2017 = ifelse(ID %in% Cloates_overlap$ID, "Cloates Sanctuary Zone", COMMENTS_2017)) %>% 
  mutate(COMMENTS_2005 = ifelse(ID %in% Cloates_overlap$ID, "Cloates Sanctuary Zone", COMMENTS_2005))

NTZ_state <- NoTake %>% 
  filter(COMMENTS_2005 %in% ("Cloates Sanctuary Zone"))
plot(NTZ_state$Spatial)

## Fixing 1987 sanctuary zones
NTZ_state <- NoTake %>% 
  filter(!is.na(COMMENTS_OLD))
NTZ_state_points <- st_centroid(NTZ_state)

NTZ_old_poly <- st_cast(st_geometry(old_MP), "LINESTRING")
NTZ_old_poly <- st_polygonize(NTZ_old_poly)

inside_NTZ <- st_covered_by(NTZ_state_points, NTZ_old_poly) %>% 
  as.data.frame()

old_NTZ_ID <- NTZ_state[as.numeric(c(inside_NTZ$row.id)),]

NoTake <- NoTake %>% 
  mutate(COMMENTS_OLD = ifelse(ID %in% old_NTZ_ID$ID, COMMENTS_OLD, NA)) 

NTZ_old_check <- NoTake %>% 
  filter(!is.na(COMMENTS_OLD)) %>% 
  ggplot(.)+
  geom_sf()+
  theme_void()

NTZ_2005_check <- NoTake %>% 
  filter(!COMMENTS_2005 %in% c("NA")) %>% 
  ggplot(.)+
  geom_sf()+
  theme_void()

NTZ_2017_check <- NoTake %>% 
  filter(!COMMENTS_2017 %in% c("NA")) %>% 
  ggplot(.)+
  geom_sf()+
  theme_void()

## Put this back together with the fished area to create "water" again

NoTake <- NoTake %>% 
  mutate(Fished_60 = "Y") %>% 
  mutate(Fished_87 = ifelse(is.na(COMMENTS_OLD), "Y", "N")) %>% 
  mutate(Fished_05 = ifelse(COMMENTS_2005 %in% c("NA"), "Y", "N")) %>% 
  mutate(Fished_17 = ifelse(COMMENTS_2017 %in% c("NA"), "Y", "N")) %>% 
  dplyr::select(-ID)

water <- rbind(Fished_area, NoTake)

## Check that the NTZs are where you expect them to be
ggplot(water)+
  geom_sf(aes(fill=Fished_17))+
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
  arrange(COMMENTS_2005) %>% 
  distinct(Spatial, .keep_all=T) %>%  # Remove any duplicated polygons if they happen to have cropped up
  st_make_valid() %>% 
  mutate(ID=row_number()) 

check <- water %>% 
  filter(str_detect(COMMENTS_2005, "Old")) %>% # Change this to check the different NTZs are where you think they should be
  ggplot(.)+
  geom_sf()+
  theme_void() 

## Create a list of the NTZ IDs
NoTakeList <- list()

NoTake87 <- water %>% 
  st_drop_geometry() %>% 
  filter(str_detect(COMMENTS_2005, "Old")) %>% 
  dplyr::select(ID) 
NoTake87 <- as.numeric(NoTake87$ID)

NoTake05 <- water %>% 
  st_drop_geometry() %>% 
  filter(Fished_05 %in% c("N")) %>% 
  dplyr::select(ID) 
NoTake05 <- as.numeric(NoTake05$ID)

NoTake17 <- water %>% 
  st_drop_geometry() %>% 
  filter(Fished_17 %in% c("N")) %>% 
  dplyr::select(ID) 
NoTake17 <- as.numeric(NoTake17$ID)

NoTakeList[[1]] <- NoTake87
NoTakeList[[2]] <- NoTake05
NoTakeList[[3]] <- NoTake17


#### CREATE THE DISTANCE MATRIX ####
## This method allows us to create a distance matrix that takes into account the land and forces the fish to swim around

# Get centroids for the grid cells 
centroids <- st_centroid_within_poly(water)
plot(centroids, cex=0.3) 

points <- as.data.frame(st_coordinates(centroids))%>% #The points start at the bottom left and then work their way their way right
  mutate(ID=row_number()) 

NCELL <- nrow(points) #Set the number of cells in the model

## Convert the points in the centroids of the ploygon to a spatial points file
points$ID <- as.integer(points$ID)
points_sf <- st_as_sf(points, coords = c("X", "Y"))
points_sp <- st_cast(st_geometry(points_sf), "POINT")

## Figure out which six points are the closest - won't all have six close ones but can start from there
dist.mat <- st_distance(points_sp) # Great Circle distance since in lat/lon

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

## Give the neighbouring points geometry based on the original set of points - don't judge me, I had intended to put this into a loop and just never got around to it
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

## Form linestrings and then multilinestrings - same thing here...
# Closest
linestrings.closest <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.closest <- st_multilinestring(do.call("rbind", linestrings.closest))

# Second closest
linestrings.second <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], second.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.second <- st_multilinestring(do.call("rbind", linestrings.second))

# Third closest
linestrings.third <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], third.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.third<- st_multilinestring(do.call("rbind", linestrings.third))

# Fourth closest
linestrings.fourth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], fourth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.fourth <- st_multilinestring(do.call("rbind", linestrings.fourth))

# Fifth closest
linestrings.fifth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], fifth.closest[x]))
  line <- st_cast(pair, "LINESTRING")
  return(line)
})
multilinestring.fifth <- st_multilinestring(do.call("rbind", linestrings.fifth))

# Sixth closest
linestrings.sixth <- lapply(X = 1:n, FUN = function(x) {
  pair <- st_combine(c(points_sp[x], sixth.closest[x]))
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
network_matrix <- st_network_cost(net, from=points_sf, to=points_sf)
network_matrix <- network_matrix * 111 # Multiple by 111 to get from degrees to kms
dim(network_matrix) # Check that the dimensions match up to how many points you think you should have in the network

## Checks to make sure it's done what you want
# Checking that short distances between points are the same e.g. 1 -> 2
test <- st_distance(points_sf[1,2], points_sf[2,2])
test * 111 # 9km which is sensible given that each hexagon is about 5km
network_matrix[1,2] # 10.92 km which is also sensible

# Checking that points that cross the land are not the same - should be longer in our network matrix
# Find two points that are either side of the land and should have a long distances between them
plot(points_sf, col = ifelse(points_sf$ID==930 | points_sf$ID==1000, "red", "black")) # These are on opposite sides of the cape

test <- st_distance(points_sf[930,2], points_sf[1000,2])
test*111 # 28km

test_matrix <- st_network_cost(net, from=points_sf[930, 2], to=points_sf[1000,2])
test_matrix*111 # 88km


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

#* Calculate the percentage of each habitat type in each of the different cells ####

## Add area of rocky habitat in each cell to water data frame - these take a long time to run!
rock_inter <- st_intersection(water, rocky)
reef_inter <- st_intersection(water, reef)
lagoon_inter <- st_intersection(water, lagoon)
pelagic_inter <- st_intersection(water, pelagic)

rock_inter <- rock_inter %>% 
  st_make_valid() %>% 
  mutate(ID = as.numeric(ID)) %>% 
  mutate(hab_area = st_area(.))
  
rocky_cells <- rock_inter %>% 
  st_make_valid(.) %>% 
  st_centroid()
plot(rocky_cells$Spatial)

rocky_cells_id <- st_intersects(water, rocky_cells) %>% # rocky_cells is the column, water is the row - tells you which water cells contain which rocky points 
  as.data.frame() %>% 
  rename(inter.row.ID = "col.id",
         water.row.ID = "row.id") %>% 
  mutate(inter.row.ID = as.numeric(inter.row.ID)) # join it back to intersection to get the habitat percentages

rocky_cells_hab <- rock_inter[as.numeric(c(rocky_cells_id$inter.row.ID)), ] %>% 
  st_drop_geometry(.) %>% 
  mutate(water.row.ID = rocky_cells_id$water.row.ID) %>% 
  group_by(water.row.ID) %>% 
  summarise(rocky.area=sum(hab_area, na.rm=T))  # sum of habitat by water ID
  
water <- water %>% # join this back to the water dataframe
  mutate(ID = as.numeric(ID)) %>% 
  mutate(rocky.area = ifelse(ID %in% as.numeric(rocky_cells_hab$water.row.ID), rocky_cells_hab$rocky.area, 0))

## Add area of reef habitat in each cell to water data frame
reef_inter <- reef_inter %>% 
  st_make_valid() %>% 
  mutate(ID = as.numeric(ID)) %>% 
  mutate(hab_area = st_area(.))

reef_cells <- reef_inter %>% 
  st_make_valid(.) %>% 
  st_centroid()
plot(reef_cells$Spatial)

reef_cells_id <- st_intersects(water, reef_cells) %>% # reef_cells is the column, water is the row - tells you which water cells contain which reef points 
  as.data.frame() %>% 
  rename(inter.row.ID = "col.id",
         water.row.ID = "row.id") %>% 
  mutate(inter.row.ID = as.numeric(inter.row.ID)) # join it back to intersection to get the habitat percentages

reef_cells_hab <- reef_inter[as.numeric(c(reef_cells_id$inter.row.ID)), ] %>% 
  st_drop_geometry(.) %>% 
  mutate(water.row.ID = reef_cells_id$water.row.ID) %>% 
  group_by(water.row.ID) %>% 
  summarise(reef.area=sum(hab_area, na.rm=T))  # sum of habitat by water ID

water <- water %>% # join this back to the water dataframe
  mutate(ID = as.numeric(ID)) %>% 
  mutate(reef.area = ifelse(ID %in% as.numeric(reef_cells_hab$water.row.ID), reef_cells_hab$reef.area, 0))

## Add area of lagoon habitat in each cell to water data frame
lagoon_inter <- lagoon_inter %>% 
  st_make_valid() %>% 
  mutate(ID = as.numeric(ID)) %>% 
  mutate(hab_area = st_area(.))

lagoon_cells <- lagoon_inter %>% 
  st_make_valid(.) %>% 
  st_centroid()
plot(lagoon_cells$Spatial)

lagoon_cells_id <- st_intersects(water, lagoon_cells) %>% # lagoon_cells is the column, water is the row - tells you which water cells contain which lagoon points 
  as.data.frame() %>% 
  rename(inter.row.ID = "col.id",
         water.row.ID = "row.id") %>% 
  mutate(inter.row.ID = as.numeric(inter.row.ID)) # join it back to intersection to get the habitat percentages

lagoon_cells_hab <- lagoon_inter[as.numeric(c(lagoon_cells_id$inter.row.ID)), ] %>% 
  st_drop_geometry(.) %>% 
  mutate(water.row.ID = lagoon_cells_id$water.row.ID) %>% 
  group_by(water.row.ID) %>% 
  summarise(lagoon.area=sum(hab_area, na.rm=T))  # sum of habitat by water ID

water <- water %>% # join this back to the water dataframe
  mutate(ID = as.numeric(ID)) %>% 
  mutate(lagoon.area = ifelse(ID %in% as.numeric(lagoon_cells_hab$water.row.ID), lagoon_cells_hab$lagoon.area, 0))

## Add area of pelagic habitat in each cell to water data frame
pelagic_inter <- pelagic_inter %>% 
  st_make_valid() %>% 
  mutate(ID = as.numeric(ID)) %>% 
  mutate(hab_area = st_area(.))

pelagic_cells <- pelagic_inter %>% 
  st_make_valid(.) %>% 
  st_centroid()
plot(pelagic_cells$Spatial)

pelagic_cells_id <- st_intersects(water, pelagic_cells) %>% # pelagic_cells is the column, water is the row - tells you which water cells contain which pelagic points 
  as.data.frame() %>% 
  rename(inter.row.ID = "col.id",
         water.row.ID = "row.id") %>% 
  mutate(inter.row.ID = as.numeric(inter.row.ID)) # join it back to intersection to get the habitat percentages

pelagic_cells_hab <- pelagic_inter[as.numeric(c(pelagic_cells_id$inter.row.ID)), ] %>% 
  st_drop_geometry(.) %>% 
  mutate(water.row.ID = pelagic_cells_id$water.row.ID) %>% 
  group_by(water.row.ID) %>% 
  summarise(pelagic.area=sum(hab_area, na.rm=T))  # sum of habitat by water ID

water <- water %>% # join this back to the water dataframe
  mutate(ID = as.numeric(ID)) %>% 
  mutate(pelagic.area = ifelse(ID %in% as.numeric(pelagic_cells_hab$water.row.ID), pelagic_cells_hab$pelagic.area, 0))

check <- water %>% 
  mutate(lagoon.area = as.numeric(lagoon.area)) %>% 
  filter(lagoon.area > 0) %>% 
  ggplot(.)+
  geom_sf()+
  theme_void()

hab_perc <- water %>% 
  mutate(cell_area = st_area(.)) %>% 
  #mutate(hab_area = as.numeric(hab_area))%>%
  mutate(cell_area = as.numeric(cell_area)) %>%
  st_drop_geometry() %>% 
  pivot_longer(cols = c(rocky.area, reef.area, lagoon.area, pelagic.area), names_to="type", values_to="hab_area") %>% 
  mutate(perc_habitat = ((hab_area/cell_area)*100)) %>%  # This tells you how much of each water grid cell is made up of the different habitat types 
  mutate(perc_habitat = ifelse(perc_habitat>100, 100, perc_habitat)) # The intersection has given a couple of places where the % is just over 100 so just round these back to 100

# Create separate data frames for each habitat type and fill in for where cells don't have a certain habitat type
pelagic_perc <- hab_perc%>%
  filter(type=="pelagic.area")%>%
  dplyr::select(ID, type, perc_habitat)%>%
  mutate(ID = as.numeric(ID))%>%
  complete(ID = 1:NCELL, fill = list(perc_habitat=0))%>%
  mutate(type = replace_na("pelagic"))

reef_perc <- hab_perc%>%
  filter(type=="reef.area")%>%
  dplyr::select(ID, type, perc_habitat)%>%
  mutate(ID = as.numeric(ID))%>%
  complete(ID = 1:NCELL, fill = list(perc_habitat=0))%>%
  mutate(type = replace_na("reef"))

rocky_perc <- hab_perc%>%
  filter(type=="rocky.area")%>%
  dplyr::select(ID, type, perc_habitat)%>%
  mutate(ID = as.numeric(ID))%>%
  complete(ID = 1:NCELL, fill = list(perc_habitat=0))%>%
  mutate(type = replace_na("rocky"))

lagoon_perc <- hab_perc%>%
  filter(type=="lagoon.area")%>%
  dplyr::select(ID, type, perc_habitat)%>%
  mutate(ID = as.numeric(ID))%>%
  complete(ID = 1:NCELL, fill = list(perc_habitat=0))%>%
  mutate(type = replace_na("lagoon"))


#### SAVE FILES FOR USE IN NEXT STEP ####
setwd(sp_dir)

saveRDS(network_matrix, file="network_matrix")
saveRDS(pelagic_perc, file="pelagic_perc")
saveRDS(reef_perc, file="reef_perc")
saveRDS(lagoon_perc, file="lagoon_perc")
saveRDS(rocky_perc, file="rocky_perc")
saveRDS(water, file="water")

setwd(sg_dir)
saveRDS(NoTakeList, file="NoTakeList")









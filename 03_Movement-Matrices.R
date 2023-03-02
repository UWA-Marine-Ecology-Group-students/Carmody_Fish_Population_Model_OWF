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

model.name <- "small"

#### READ FILES ####
setwd(sp_dir)
Water <- readRDS(paste0(model.name, sep="_", "water"))

# Habitat Layers
setwd(sp_dir)
habitat.types <- c("Reef", "Lagoon", "RockReef" ,"Pelagic") # These need to be the same as the habitat types in the names of the shapefiles you want to load, if it's all in one then ignore the loop and load the file

habitat <- data.frame()

for (i in 1:4){
  
  hab.file <- st_read(paste0(habitat.types[i], "Habitat.gpkg", sep="")) %>% 
    st_transform(4283)%>%
    st_make_valid%>%
    st_crop(xmin=113.8826, xmax=113.9967, ymin=-22.0725, ymax=-21.84793)%>% # This needs to be the same as in the first script
    mutate(type = habitat.types[i]) %>% 
    dplyr::select(type, geom)
  
  habitat <- rbind(habitat, hab.file)
  
}


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
dist.mat <- st_distance(points_sp) # Great Circle distance since in lat/long

## Get the IDS for the cells that are neighbours
n.closest <- 6 # Decide how many neighbours you want to use

neighbours <- as.data.frame(array(0, dim=c(NCELL, n.closest)))

for (i in 1:n.closest){
  neighbours[,i] <- apply(dist.mat, 1, function(x) { order(x, decreasing=F)[i+1] })
}


## Give the neighbouring points geometry based on the original set of points
point.list <- list()

for (i in 1:n.closest){
  
  temp1 <- as.data.frame(neighbours[,i])
  
  temp2 <- temp1 %>%
    rename(ID = "neighbours[, i]") %>% 
    inner_join(., points, by="ID") %>% 
    st_as_sf(., coords = c("X", "Y")) 
  
  temp3 <- st_cast(st_geometry(temp2), "POINT")
  
  point.list[[i]] <- temp3
  
}

### CONNECT THE POINTS TO THEIR NEIGHBOURS TO FORM A NETWORK ###
n <- nrow(points)

## Form linestrings and then multilinestrings
multilinestrings <- list()

for(i in 1:n.closest){
  
  linestring <- point.list[[i]]
  
  temp <- lapply(X = 1:n, FUN = function(x) {
    pair <- st_combine(c(points_sp[x], linestring[x]))
    line <- st_cast(pair, "LINESTRING")
    return(line)
    })
    
  temp2 <- st_multilinestring(do.call("rbind", temp))
  
  multilinestrings[[i]] <- temp2
  
}

connected <- st_combine(c(multilinestrings[[1]], multilinestrings[[2]], multilinestrings[[3]], multilinestrings[[4]], multilinestrings[[5]], multilinestrings[[6]]))
connected <- st_cast(connected, "LINESTRING") # Needs to be a line string rather than multiline for the next step

setwd(working.dir)
st_write(connected, paste0(model.name, sep="_", "network.shapefile.shp"))

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
#* Calculate the percentage of each habitat type in each of the different cells ####

## Add area of habitat type in each cell to water data frame 
habitat.types <- c("Reef", "Lagoon", "RockReef" ,"Pelagic")

habitat_inter <- list()

for(i in 1:4){
  
  temp <- habitat %>% 
    filter(type %in% habitat.types[i])
  
  inter <- st_intersection(water, temp) %>% 
    st_make_valid() %>% 
    mutate(ID = as.numeric(ID)) %>%
    mutate(hab_area = as.numeric(st_area(.)*0.000001)) 
  
  habitat_inter[[i]] <- inter
  
}

check_hab <- rbind(habitat_inter[[1]], habitat_inter[[2]], habitat_inter[[3]], habitat_inter[[4]]) %>% 
  mutate(cell_area = as.numeric(cell_area*0.000001))
plot(check_hab$Spatial)
water1 <- water %>% 
  mutate(cell_area = cell_area*0.000001)

cell1_water <- water1 %>% 
  filter(ID==9)
plot(water1$Spatial)
plot(cell1_water$Spatial, add=T, col="blue")
cell1_hab <- check_hab %>% 
  filter(ID==9)
plot(cell1_hab$Spatial, add=T, col="red")

water_hab <- check_hab %>% 
  pivot_wider(id_col=ID, names_from=type, values_from=hab_area) %>% 
  left_join(., water, by="ID") %>% 
  mutate(cell_area = (as.numeric(cell_area)*0.000001)) %>% 
  rowwise() %>% 
  mutate(hab_sum = (sum(Reef, Lagoon, RockReef, Pelagic, na.rm = T))) %>% 
  mutate(correct_hab = ifelse(hab_sum>as.numeric(cell_area), "N", "Y")) %>% 
  mutate(difference = hab_sum-cell_area) %>% 
  ungroup()
  

check <- water_hab %>% 
  filter(Lagoon>0)

plot(check$Spatial)

hab_perc <- water_hab %>% 
  dplyr::select(-Spatial) %>% 
  pivot_longer(cols = c(RockReef, Reef, Lagoon, Pelagic), names_to="type", values_to="hab_area") %>% 
  mutate(perc_habitat = ((hab_area/cell_area)*100)) %>%  # This tells you how much of each water grid cell is made up of the different habitat types 
  mutate(perc_habitat = ifelse(perc_habitat>100, 100, perc_habitat)) # The intersection has given a couple of places where the % is just over 100 so just round these back to 100

# Create separate data frames for each habitat type and fill in for where cells don't have a certain habitat type
perc_by_hab <- list()

for(i in 1:length(habitat.types)){
  
  temp <- hab_perc%>%
    filter(type==paste0(habitat.types[i]))%>%
    dplyr::select(ID, type, perc_habitat)%>%
    mutate(ID = as.numeric(ID))%>%
    complete(ID = 1:NCELL, fill = list(perc_habitat=0))%>%
    mutate(type = replace_na(habitat.types[i]))
  
  perc_by_hab[[i]] <- temp
}
check <- cbind(perc_by_hab[[1]], perc_by_hab[[2]], perc_by_hab[[3]], perc_by_hab[[4]])

#### SAVE FILES FOR USE IN NEXT STEP ####
setwd(sp_dir)

#### NEED TO REDO THE NINGALOO ONES BECAUSE YOU'VE OVERWRITTEN THEM

saveRDS(network_matrix, file=paste0(model.name, sep="_", "network_matrix"))
saveRDS(perc_by_hab[[1]], file=paste0(model.name, sep="_", "reef_perc"))
saveRDS(perc_by_hab[[2]], file=paste0(model.name, sep="_", "lagoon_perc"))
saveRDS(perc_by_hab[[3]], file=paste0(model.name, sep="_", "rocky_perc"))
saveRDS(perc_by_hab[[4]], file=paste0(model.name, sep="_", "pelagic_perc"))
saveRDS(water, file=paste0(model.name, sep="_", "water"))

#### CREATE MATRIX OF CONNVECTIVITY FOR FISH MOVEMENT ####
## Assume the fish will try and swim the shortest path between locations
# Calculate the probability a fish moves to this site in a given time step using the swimming speed we set earlier 
# This creates a dispersal kernel based on the negative exponential distribution
# Some of these loops take quite a while to run
setwd(sp_dir)

network_matrix <- readRDS(paste0(model.name, sep="_", "network_matrix"))
reef_perc <- readRDS(paste0(model.name, sep="_", "reef_perc"))
lagoon_perc <- readRDS(paste0(model.name, sep="_", "lagoon_perc"))
rocky_perc <- readRDS(paste0(model.name, sep="_", "rocky_perc"))
pelagic_perc <- readRDS(paste0(model.name, sep="_", "pelagic_perc"))



pDist <- matrix(NA, ncol=NCELL, nrow=NCELL)

for(r in 1:NCELL){
  for(c in 1:NCELL){
    p <- network_matrix[r,c]*1
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

### Save all the files you need to remake these matrices ###
setwd(sp_dir)

saveRDS(pDist, file=paste0(model.name, sep="_", "pDist"))
saveRDS(pPelagic, file=paste0(model.name, sep="_","pPelagic"))
saveRDS(pReef, file=paste0(model.name, sep="_","pReef"))
saveRDS(pLagoon, file=paste0(model.name, sep="_","pLagoon"))
saveRDS(pRocky, file=paste0(model.name, sep="_","pRocky"))

#### CREATE PROBABILITY OF MOVEMENT USING UTILITY FUNCTION ####
## Load files if need be ##
setwd(sp_dir)

pDist <- readRDS(paste0(model.name, sep="_", "pDist"))
pPelagic <- readRDS(paste0(model.name, sep="_","pPelagic"))
pReef <- readRDS(paste0(model.name, sep="_","pReef"))
pLagoon <- readRDS(paste0(model.name, sep="_","pLagoon"))
pRocky <- readRDS(paste0(model.name, sep="_","pRocky"))

# First determine the utility of each of the sites 
# This is very sensitive to changes in the values particularly for reef 
# PROBABLY ALSO NEED TO PUT DEPTH IN HERE

SwimSpeed <- 1.6 #2.5

a = -SwimSpeed
b =  0.15 #0.1 # Reef
c =  0.09 #0.09 #Lagoon
d =  0.003 #0.003 #Rocky Reef
e =  0.001 #0.001 #Pelagic

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
  dplyr::select(perc_habitat, ID)
dispersal$area <- as.vector(Water$cell_area*0.000001)

dispersal <- dispersal %>% 
  filter(perc_habitat!=0) %>% 
  mutate(area = area*0.1,
         perc_habitat = perc_habitat*0.5)

recruitment <- array(0, dim=c(nrow(dispersal), 2))
recruitment[ ,2] <- as.numeric(dispersal$ID) 

cellU <- matrix(0, ncol=2, nrow=(nrow(dispersal)))
cellU[,2] <- as.numeric(dispersal$ID) 


for(cell in 1:nrow(recruitment)){
  U <- exp(dispersal[cell,1])
  cellU[cell,1] <- as.numeric(U)
}

rowU <- as.data.frame(sum(cellU[,1]))

for (cell in 1:nrow(recruitment)){
  recruitment[cell,1] <- cellU[cell, 1]/rowU[1,1]
}
colSums(recruitment)

recruitment <- as.data.frame(recruitment) %>% 
  rename(ID = "V2")

recruitment <- merge(recruitment, lagoon_perc, by="ID", all=T) %>% #check that cells with no lagoon habitat have 0 probability of recruitment
  mutate_all(~replace(., is.na(.), 0)) #For cells where there was no lagoon habitat put probability of recruitment as 0

recruitment <- as.vector(recruitment[,2])

#### RECRUIT MOVEMENT ####
## Want the recruits to stay in the lagoon until they mature and move to the reef

a = -SwimSpeed #-SwimSpeed
b = 0.01 #0.01
c = 0.09 #0.09
d = 0.05 #0.05
e = 0.001 #0.001

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
saveRDS(Pj, file=paste0(model.name, sep="_", "movement"))
saveRDS(ProbRec, file=paste0(model.name, sep="_","juvmove"))
saveRDS(recruitment, file=paste0(model.name, sep="_","recruitment"))
saveRDS(water, file=paste0(model.name, sep="_", "water"))


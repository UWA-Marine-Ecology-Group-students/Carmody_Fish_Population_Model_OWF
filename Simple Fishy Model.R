library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(stringr)
library(lwgeom)


#DATA + DIRECTORIES
working.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

## Save these directory names to use later----
data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
cm_dir <- paste(working.dir, "Connectivity_Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial Data", sep="/")

#### PARAMETER VALUES ####
StartN <- 2000 # Number of fish in each cell at the beginning of the model
M <- 0.146 # Natural mortality rate, Marriot et al. 2011
totaltime <- 20 # Number of time steps, each time step is half a day

#Ricker recruitment model parameters (these are currently just made up values)
a <- 7
b <- 0.0017
M50 <- 2 # From Grandcourt et al. 2010
M95 <- 5 # From Grandcourt et al. 2010 (technically M100)

#Fishing mortality parameters
q <- 0.3 # Made this value up
E <- 10 # Can get real value in days per year from Marriott et al. 2012
A50 <- 4 # For L. miniatus from Williams et al. 2010
A95 <- 6 # For L. miniatus from Williams et al. 2012

Fishing <- E*q

#### SPATIAL DATA ####

# THis returns the centre of the ploygon, but if it's on land it will create a new centroid

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

## Load spatial files for area your model covers
setwd(sp_dir)
wa_map <- st_read("WACoastline.shp")%>%
  st_transform(4283)%>%
  st_make_valid

ningaloo_map <- st_crop(wa_map, xmin=112.5, xmax=114.3, ymin=-24, ymax=-21)
plot(ningaloo_map$geometry)

# Make grid cells for fish to live in
grd <- st_make_grid(ningaloo_map, cellsize=0.1, square=FALSE)%>%
  st_crop(xmin=112.5, xmax=114.3, ymin=-24, ymax=-21) #make sure extent of grid is the same as the polygon
plot(grd, add=TRUE)

water <- st_difference(grd, ningaloo_map)
plot(water)
# turn water into a data frame for easier use
water <- st_sf(water)

water <- water%>%
  mutate(ID=row_number())%>%
  mutate(ID=factor(ID))

# Get centroids for the grid cells 
centroids <- st_centroid_within_poly(water)
plot(centroids, add=TRUE) #Something is wrong with where the centroids are now

points <- as.data.frame(st_coordinates(centroids))%>% #The points start at the bottom left and then work their way their way right
  mutate(ID=row_number())
  
# Calculate a distance matrix between the cells 
  # To make them go round the land, create visibility graph (grass function), then distance between nodes, then pick the
  # shortest path
dist_matrix <- as.matrix(dist(points)) #This is just the linear distance  between the centroids of the grid cells

# Bring in the habitat data 
setwd(sp_dir)
habitat <- st_read("SeamapAus_WA_DPAW_marine_habitatsPolygon.shp")%>%
  st_transform(4283)%>%
  st_make_valid()

ningaloo_habitat <- st_crop(habitat, xmin=112.5, xmax=114.3, ymin=-24, ymax=-21)
plot(ningaloo_habitat$geometry, add=TRUE)

# Set the different habitat types to be a value between 0 and 4
  # 0 is places fish are unlikely to be found e.g.  mudflats and saltmarsh
  # 2 is unnatractive places for fish including mobile sand, deep water habitat, pelagic
  # 5 is macroalgae and mangroves where you're likely to get juvenile fish not targeted by fishers
  # 10 is coral reef and subtidal bare reef (subtidal and intertidal)

ningaloo_habitat_cls <- ningaloo_habitat%>%
  mutate(classification = ifelse(str_detect(SM_HAB_CLS, "Mudflat|Saltmarsh"), 0,
                                ifelse(str_detect(SM_HAB_CLS, "Mobile sand|Pelagic"), 2,
                                       ifelse(str_detect(SM_HAB_CLS, "Macroalgae|Mangroves"), 5, 10))))

# Create an intersection of the habitat and the grid cells to see what overlaps where, then
# calculate the area of each polygon
intersection <- st_intersection(ningaloo_habitat_cls, water)
intersection$area <- st_area(intersection)

#calculate the area of the water grid cells
water$area <- st_area(water)

#calculate the percentage of each habitat type found in the different water grid cells
hab_perc <- intersection%>%
  st_drop_geometry()%>%
  as.data.frame()%>%
  mutate(classification=factor(classification))%>%
  mutate(ID=factor(ID))%>%
  dplyr::select(classification, ID, area)%>%
  group_by(ID, classification)%>%
  summarise(total_area = sum(area))

# Join the habitat percentage data to the grid cell data and calculate percent reef
per_reef <- water%>%
  st_drop_geometry()%>%
  full_join(hab_perc)%>%
  filter(classification==10)%>%
  mutate(percent_reef=(total_area/area)*100)%>%
  full_join(water, by="ID")%>%
  dplyr::select(ID, percent_reef)


# For cells where we are missing data about habitat take raster data and then calculate
# the average rugosity, for now we'll just say anything over XXX is reef
setwd(sp_dir)
bathy <- raster("Carnarvon_Shelf_Bathymetry_3_2008_3m_cog.tiff")
bathy <- flip(bathy, direction="y")
proj4string(bathy) #check projection



water_ras <- st_transform(water, "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

grd_transform <- st_transform(grd, "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

plot(bathy)
plot(water_ras$water, add=TRUE)

st_crs(grd_transform)

#### CREATE CONNECTIVITY COEFFICIENT ####

# This needs to be habitat matrix*coefficient + distance*coefficient in some form or another. Need to figure out how to decide
# what the coefficients should be so we can weight the importance of habitat relative to distance

## Settlement matrix
setwd(cm_dir)
Settlement <- read.csv("Settlement_Matrix.csv", header=F)

Settlement <- as.matrix(Settlement)
rownames(Settlement) <- c("Site 1", "Site 2", " Site 3")
colnames(Settlement) <- "Abundance"

## Movement matrix for 1 year old fish
Conn_Year_1 <- read.csv("Conn_Year_1.csv", header=F)

Conn_Year_1 <- as.matrix(Conn_Year_1)
colnames(Conn_Year_1) <- c("Site 1", "Site 2", " Site 3")
rownames(Conn_Year_1) <- c("Site 1", "Site 2", " Site 3")

## Movement matrix for fish 2-5 years
Conn_Year_2 <- read.csv("Conn_Year_2.csv", header=F)

Conn_Year_2 <- as.matrix(Conn_Year_2)
colnames(Conn_Year_2) <- c("Site 1", "Site 2", " Site 3")
rownames(Conn_Year_2) <- c("Site 1", "Site 2", " Site 3")

## Movement matrix for fish 5-7 years
Conn_Year_3 <- read.csv("Conn_Year_3.csv", header=F)

Conn_Year_3 <- as.matrix(Conn_Year_3)
colnames(Conn_Year_3) <- c("Site 1", "Site 2", " Site 3")
rownames(Conn_Year_3) <- c("Site 1", "Site 2", " Site 3")


## Movement matrix
Movement <- array(c(Conn_Year_1, Conn_Year_2, Conn_Year_3), dim= c(3,3,3))

## Set up the model
Days <- array(0, dim = c(3,7,7)) # This would have to be t+1 as you need the matrix to run for as many years as the model will run
rownames(Days) <- rownames(Fish)      
colnames(Days) <- seq(1,7)
Days[,1,1] <- 200
Days[,1,2] <- c(10,25,15)
Days[,1,3] <- c(35,20,15)
Days[,1,4] <- c(100,75,145)
Days[,1,5] <- c(160,55,20)
Days[,1,6] <- c(10,5,25)
Days[,1,7] <- c(5,5,40)
Recruits <- 200

## Run the model
for(t in 2:7){
  for(d in 1:dim(Days)[3]){
    if(d==1){Days[,t-1,d] <- Recruits}
    else if (d>=6) {Days[,t,d] <-  Movement[,,3] %*% Days[,t-1,d-1]}
    else {Days[,t,d] <- Movement[,,2] %*% Days[,t-1,d-1]}
  }
  Recs <- matrix(0, nrow=1, ncol=1)
  # Mortality
  for(d in 1:dim(Days)[3]){
    sa <- 1/(1+(exp(-log(19)*((d-A95/A95-A50)))))
    Days[,t,d] <- Days[,t,d]*exp(-M)*exp(-sa*Fishing)
  }
  # Ricker recruitment model with maturation
  Recruitment <- as.matrix(colSums(Days[,t,], dims=1)) 
  for(Age in 1:dim(Recruitment)[1]){
    Mature <- 1/(1+(exp(-log(19)*((Age-M95)/(M95-M50)))))
    S <- colSums(Recruitment)
    Rec <- a*S*exp(-b*S)*exp(-Mature)
    Recs <- rbind(Recs, Rec)
  }
  R <- colSums(Recs)
  Recruits <- as.matrix(Settlement*R)
}
Days


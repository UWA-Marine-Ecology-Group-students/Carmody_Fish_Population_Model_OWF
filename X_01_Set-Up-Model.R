###################################################

# Script for setting up the test fishing surface 
# Area has been made smaller so things can be run faster
# Sets up all the files from scripts 01 and 02 

###################################################
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


#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
m_dir <- paste(working.dir, "Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")

#### LOAD FILES ####
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
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-23.2) %>% 
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
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-23.2)

NTZ <- MP%>%
  filter(IUCN == "IA") %>% 
  filter(!COMMENTS == "Cloates Sanctuary Zone") %>% 
  dplyr::select(PA_ID,COMMENTS, geometry) %>% 
  filter(!COMMENTS %in% c("North Muiron Conservation Area","South Muiron Conservation Area",
                          "Previously PA_ID WA_42756 in terrestrial CAPAD", "Conservation Area"))
plot(NTZ)

# # Map of Commonwealth Marine Park - NOT NEEDED AS I'VE CROPPED THIS OUT
# AMP <- st_read("NTZ and Fished areas for status.shp") %>% 
#   st_transform(4283)%>%
#   st_make_valid%>%
#   st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-23.2)
# plot(AMP)

# AMP_NTZ <- AMP %>% 
#   filter(ResName == "Ningaloo"|COMMENTS == "Cloates Sanctuary Zone") %>% 
#   dplyr::select(PA_ID, COMMENTS, geometry) %>% 
#   mutate(COMMENTS = ifelse(is.na(COMMENTS), "Comm Cloates", COMMENTS))
# plot(AMP_NTZ$geometry)

NTZ <- rbind(NTZ, old_MP) 
NTZ <- st_make_valid(NTZ)

plot(NTZ$geometry)

# Habitat Layers
reef <- st_read("ReefHabitat.gpkg")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-23.2)%>% 
  dplyr::select(geom)

reef$type <- "reef"

lagoon <- st_read("LagoonHabitat.gpkg")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-23.2) %>% 
  dplyr::select(geom)

lagoon$type <- "lagoon"

rocky <- st_read("RockReefHabitat.gpkg")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.3, ymin=-24, ymax=-23.2)%>% 
  dplyr::select(geom)

rocky$type <- "rocky"

pelagic <- st_read("PelagicHabitat.gpkg")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.3, ymin=-24, ymax=-23.2)

pelagic$type <- "pelagic"

#### SPATIAL DATA ####

## Create extent of area you want to cover 

ningaloo_map <- st_crop(wa_map, xmin=100, xmax=114.65, ymin=-24, ymax=-23.2) #NEED TO CHANGE THIS TO BE TALLER TO STOP DROPPING CELLS
plot(ningaloo_map$geometry, col="#8ec3ca")

# Make grid cells for fish to live in
grd <- st_make_grid(ningaloo_map, cellsize=0.05, square=FALSE) %>%
  st_crop(xmin=100, xmax=114.64, ymin=-24, ymax=-23.21) #make sure extent of grid is the same as the polygon
plot(grd, add=TRUE)

water <- st_difference(grd, ningaloo_map)
plot(water, border="#8ec3ca")

# Adjust the sizes of the grid cells so that the ones closer to the shore are smaller 
GrdSmall <- st_make_grid(ningaloo_map, cellsize=0.025, square=FALSE) %>%
  st_crop(xmin=100, xmax=114.64, ymin=-24, ymax=-23.21) #make sure extent of grid is the same as the polygon
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
water <- st_make_valid(water)
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

a = 0.08
b = 0.15
c = 0.1
d = 0.04
e = 0.01

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
dispersal$area <- as.vector(water$area)

dispersal <- dispersal %>% 
  filter(perc_habitat!=0)

temp <- array(0, dim=c(NCELL, 1))
recruitment <- array(0, dim=c(nrow(dispersal), 2))
recruitment[ ,2] <- as.numeric(dispersal$ID) 

for (CELL in 1:nrow(dispersal)){
  
  for (cell in 1:nrow(dispersal)){
    temp[cell, 1] <- as.numeric((dispersal[cell,1]/dispersal[cell,3]))    
  }

  summed <- sum(temp)
  
  recruitment[CELL,1] <- as.numeric(temp[CELL,1]/summed)
}

recruitment <- as.data.frame(recruitment)
colnames(recruitment)[colnames(recruitment)=="V2"] <- "ID"

recruitment <- merge(recruitment, lagoon_perc, by="ID", all=T) %>% #check that cells with no lagoon habitat have 0 probability of recruitment
  mutate_all(~replace(., is.na(.), 0)) #For cells where there was no lagoon habitat put probability of recruitment as 0
  


#### RECRUIT MOVEMENT ####
## Want the recruits to stay in the lagoon until they mature and move to the reef

a = 0.1
b = 0.01
c = 0.09
d = 0.05
e = 0.001

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
saveRDS(Pj, file="test_movement")
saveRDS(water, file="test_water")
saveRDS(ProbRec, file="test_juvmove")
saveRDS(recruitment, file="recruitment")

#### SET UP THE FISHING SURFACE - SAME AS SCRIPT 02 ####
## Data
setwd(data_dir)
boat_days <- read.csv("Boat_Days_Gascoyne.csv")

boat_days <- boat_days%>%
  mutate(NumMonth = as.numeric(NumMonth)) %>% 
  mutate(Month = as.factor(Month)) %>% 
  mutate(Gascoyne_Boat_Days = as.numeric(Gascoyne_Boat_Days)) %>% 
  mutate(Month = fct_relevel(Month, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))

## Spatial Data
setwd(sg_dir)
water <- readRDS("test_water")

# Locations of the boat ramps
setwd(sp_dir)

BR <- st_read("Boat_Ramps.shp") %>% 
  st_transform(4283)%>%
  st_make_valid

#### SET UP FISHING SURFACE ####

# Work out the number of boat days for the Exmouth region
# For this we're using data from the whole of the Gascoyne and then splitting up the number of boat days in Exmouth compared
# to Carnarvon/Shark Bay and then using visitor numbers to work out what percentage of the baot days were in Exmouth
# LGA statistics for visitors to the biggest shires where people would go fishing in the Gascoyne indicate that trips to
# Exmouth are about 33% of the visitation to the region so we'll allocate 33% of the Boat Days to Exmouth
# We then need to create a model for hindcasting fishing effort back to something like the 1960s to get an estimate of what 
# effort might have been like for the years where we don't have any data 

## Working out the fishing effort in the Gascoyne
# Plotting the data to see what it looks like
TotalYear <- boat_days %>% 
  group_by(Year) %>% 
  summarise(., Total_Boat_Days=sum(Gascoyne_Boat_Days))

YearPlot <- ggplot(TotalYear) + 
  geom_point(aes(x=Year, y=Total_Boat_Days))

YearModel <- lm(Total_Boat_Days~Year, data=TotalYear)
summary(YearModel)

#### HINDCASTING ####
Year2011_1990 <- as.data.frame(array(0, dim=c(21,1))) %>% 
  mutate(V1=seq(1990, 2010, by=1)) %>% 
  rename("Year"=V1)

predictions <- predict(YearModel, newdata=Year2011_1990)

Year2011_1990 <- Year2011_1990 %>% 
  mutate(Total_Boat_Days = predictions)

boat_days_hind <- rbind(TotalYear, Year2011_1990)

effort <- seq(0, 118466, length=30)
years <- seq(1960, 1989, by=1)

Years_1960_1989 <- as.data.frame(cbind(years, effort)) %>% 
  rename("Year" = years) %>% 
  rename("Total_Boat_Days" = effort)

#### JOIN FISHING EFFORT ####
boat_days_hind <- rbind(boat_days_hind, Years_1960_1989)

YearPlot <- ggplot(boat_days_hind) + 
  geom_point(aes(x=Year, y=Total_Boat_Days))

#### ALLOCATING MONLTHY EFFORT ####

# Work out proportion of fishing effort in each month for each year
boat_days <- boat_days %>% 
  group_by(Year) %>% 
  mutate(Year_Sum = sum(Gascoyne_Boat_Days)) %>%
  mutate(Month_Prop = Gascoyne_Boat_Days/Year_Sum) %>% 
  dplyr::select(-Year_Sum)

# Work out the average fishing effort for each month
boat_days <- boat_days %>% 
  group_by(Month) %>% 
  mutate(Ave_Month_Prop = mean(Month_Prop))

Month_Prop_Ave <- boat_days[1:12,c(2,6)]

check <- boat_days %>% 
  group_by(Year) %>% 
  summarise(Ave_Month_Prop = sum(Ave_Month_Prop))

check <- boat_days %>% 
  group_by(Year) %>% 
  summarise(Ave_Month = mean(Gascoyne_Boat_Days))

# Split up the hindcast data into monthly means 
boat_days_hind <- boat_days_hind %>% 
  bind_rows(replicate(11, boat_days_hind, simplify = FALSE)) %>% # Make a row for each month of each year in the hindcast
  mutate(Year = as.integer(Year)) %>% 
  arrange(-Year) %>% 
  mutate(NumMonth = (rep(1:12, 59))) %>% # Add numerical month
  mutate(Month = rep(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), 59)) %>% #Add word month
  filter(Year<2011) 


boat_days_hind <- boat_days_hind %>% 
  mutate(Gascoyne_Boat_Days = 0) %>% 
  left_join(., Month_Prop_Ave) %>% 
  group_by(Year) %>% 
  mutate(Gascoyne_Boat_Days = Total_Boat_Days*Ave_Month_Prop) %>% 
  ungroup() %>% 
  dplyr::select(-c(Total_Boat_Days))

check <- boat_days_hind %>% 
  group_by(Year) %>% 
  summarise(., sum(Gascoyne_Boat_Days))

# Put it all together into one big dataframe

Full_Boat_Days <- boat_days[,1:4] %>% 
  rbind(boat_days_hind) %>% 
  arrange(-Year)


# Plot and check that it looks right
MonthPlot <- Full_Boat_Days %>%
  mutate(Unique = paste(Year, NumMonth, sep="_")) %>%
  filter(Month %in% c("Jun")) %>%
  ggplot() +
  geom_point(aes(x=Unique, y=Gascoyne_Boat_Days))

#### ALLOCATION EFFORT TO EXMOUTH AND BOAT RAMPS ####
## We're using tourism numbers to allocate to Exmouth - saying about 65% of the trips to the region are to Ningaloo

# Create data frame that is just Exmouth Boat Days 
Exmouth_Boat_Days <- Full_Boat_Days %>% 
  mutate(Boat_Days = Gascoyne_Boat_Days*0.65) %>% 
  dplyr::select(-Gascoyne_Boat_Days)

# Split up the effort to Boat Ramps according to the information we have collected from Exmouth about how often people
# Use the boat ramps - the Exmouth Marina wasn't constructed until 1997 but there was a ramp at town beah from the 60s 
# onwards, this was built at the same time that the Tantabiddi ramp was built 
# Bundegi boat ramp wasn't built until 2008 and before that was just a beach with some concrete covered in sand
# So might be a good idea to reduce the number of people we assume launched from there as it was effectively a beach launch
# You also need to make sure that you standardise by number of hours spent at each boat ramp to account for sampling effort

# Tantabidi had 149.48 hours of sampling effort with 224 trips
# Bundegi had 133.9 hours of sampling effort with 157
# Exmouth Marina had 166.15 hours of sampling effort with 198 trips
# Coral Bay had 146.07 hours of sampling effort with 151 trips

BR_Trips <- data.frame(BoatRamp = c("Tantabiddi", "Bundegi", "ExmouthMar", "CoralBay"),
                       Effort = c(194.48, 133.9, 166.15, 146.07),
                       Trips = c(224, 157, 198, 151))

BR_Trips <- BR_Trips %>% 
  mutate(Trip_per_Hr = as.numeric(unlist((Trips/Effort)))) %>% # Standardise the no. trips based on how much time you spent sampling
  mutate(BR_Prop = Trip_per_Hr/sum(Trip_per_Hr)) %>% #Then work out the proportion of trips each house that leave from each boat ramp
  mutate(BR_Prop_08 = c(0.3031544, 0.1577101, 0.3119251, 0.2272104)) # Create separate proportions for Bundegi before 2008, have 
# allocated 10% of its boats to the other Exmouth Boat Ramps

# Create a loop that allocates the correct proportion of boat effort to each of the BRs accounting for reduced effort at Bundegi before ther ramp was built
for(Y in 1:708){
  if(Exmouth_Boat_Days[Y,1]<2008){ # This is for when Bundegi wasn't a proper ramp so probably would have had less effort
    Exmouth_Boat_Days$Tb_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[1,6]
    Exmouth_Boat_Days$Bd_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[2,6]
    Exmouth_Boat_Days$ExM_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[3,6] 
    Exmouth_Boat_Days$CrB_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[4,6]
  } else {
    Exmouth_Boat_Days$Tb_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[1,5] 
    Exmouth_Boat_Days$Bd_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[2,5] 
    Exmouth_Boat_Days$ExM_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[3,5] 
    Exmouth_Boat_Days$CrB_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[4,5]
  }
}

check <- Exmouth_Boat_Days %>% 
  mutate(Total = Tb_BR+Bd_BR+ExM_BR+CrB_BR)

#### ALLOCATING EFFORT TO CELLS ####

## Work out the probability of visiting a cell from each boat ramp based on distance and size
BR <- st_sf(BR)

water <- water %>% 
  mutate(DistBR = 0)

# Working out distance from each BR to each cell
DistBR <- as.data.frame(array(0, dim = c(NCELL,3))) # This will contain the distance from each boat ramp to every cell

for(CELL in 1:NCELL){
  
  for(RAMP in 1:4){
    x <- st_distance(centroids[CELL,1], BR[RAMP,3])
    DistBR[CELL,RAMP] <- (x/1000) #to get the distance in km
    
  }
}

# Create a data frame with both the distances and the areas of the cells
Cell_Vars <- DistBR %>% 
  mutate(Area = as.vector((water$area)/1000)) %>% #Cells are now in km^2 but with no units
  rename("Bd_BR"=V1) %>% 
  rename("ExM_BR" = V2) %>% 
  rename("Tb_BR" = V3) %>% 
  rename("CrB_BR"=V4)

## Now need to create a separate fishing surface for each month of each year based on distance to boat ramp, size of each
## cell and multiply that by the effort in the cell to spatially allocate the effort across the area
BR_U <- as.data.frame(matrix(0, nrow=NCELL, ncol=4)) #Set up data frame to hold utilities of cells

BR_U <- BR_U %>% #make sure you give the columns good names so that you know what they are
  rename("Bd_U"=V1) %>% 
  rename("ExM_U" = V2) %>% 
  rename("Tb_U" = V3) %>% 
  rename("CrB_U"=V4)


for(RAMP in 1:4){
  
  Tot <- sum(Cell_Vars$Area/Cell_Vars[,RAMP])
  
  for(CELL in 1:NCELL){
    BR_U[CELL,RAMP] <- (Cell_Vars[CELL,RAMP]/Tot)
  } 
} 

# Now we have BR_U which has the "utilities" for each cells based on it's size and distance from BR  

BR_Trips <- Exmouth_Boat_Days %>% # This is just the trips from each boat ramp
  ungroup() %>% 
  mutate(NumYear = rep(59:1, each=12)) %>% #This is to turn the years into a count for the loop
  dplyr::select(NumYear, NumMonth, Bd_BR, ExM_BR, Tb_BR, CrB_BR)  


Fishing <- array(0, dim=c(NCELL, 12, 59)) #This array has a row for every cell, a column for every month, and a layer for every year
Months <- array(0, dim=c(NCELL, 12))
Ramps <- array(0, dim=c(NCELL, 4))

for(YEAR in 1:59){
  
  for(MONTH in 1:12){
    
    for(RAMP in 1:4){
      
      temp <- BR_Trips %>% 
        filter(NumYear==YEAR) %>% 
        dplyr::select(-c(NumYear, NumMonth))
      
      temp <- as.matrix(temp) 
      
      Ramps[ ,RAMP] <- BR_U[ ,RAMP] * temp[MONTH,RAMP]
    }
    
    Months[,MONTH] <-  rowSums(Ramps)
  }
  Fishing[ , ,YEAR] <- Months 
}

Fishing <- Fishing*q # multiply by catchability

## Add NTZs
NoTake <- st_sf(water) %>% 
  st_drop_geometry() %>% 
  dplyr::select(Fished, ID, COMMENTS)

# Creating new columns for each of our year groups when new NTZs were put in place (easier than coding it in the model)
NoTake <- NoTake %>% 
  rename(Fished60_87 = "Fished") %>% 
  mutate(Fished60_87 = "Y") %>% 
  mutate(Fished87_05 = ifelse(str_detect(COMMENTS, c("Old")), "N", "Y")) %>% 
  mutate(Fished05_18 = ifelse(str_detect(COMMENTS, "[:alpha:]") & !str_detect(COMMENTS, ("Comm Cloates")), "N", "Y")) %>% 
  mutate(Fished18_21 = ifelse(str_detect(COMMENTS, "[:alpha:]"), "N", "Y")) %>% 
  dplyr::select(ID, Fished60_87, Fished87_05, Fished05_18, Fished18_21) %>% 
  mutate_at(vars(contains("Fished")), ~replace(., is.na(.), "Y"))

## Set the effort in the corresponding cells to be 0s in the years where the sanctuary zones are in place
for(YEAR in 1:59){
  for(CELL in 1:NCELL){
   
    if(Year>=27&Year<=44){
      if(NTZ[CELL, 3]=="N"){Fishing[CELL, ,YEAR] <- 0}
    }
    else if(Year>=45&Year<=57){
      if(NTZ[CELL, 4]=="N"){Fishing[CELL, ,YEAR] <- 0}
    }
    else if (Year>57) {
      if(NTZ[CELL, 5]=="N"){Fishing[CELL, ,YEAR] <- 0}
    } 
  }
}

#### SAVE DATA ####
setwd(sg_dir)

saveRDS(Fishing, file="test_fishing")
#saveRDS(NoTake, file="test_NoTake")


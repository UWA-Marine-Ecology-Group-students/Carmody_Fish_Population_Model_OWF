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
levels_water <- data.frame(c("<1000",  "1000-1100", "1100-1200", "1200-1300",
                             "1300-1400", "1400-1500", "1500-1600", ">1600"))
names(levels_water)[1] <- "Levels"
levels_water$Levels <- as.factor(levels_water$Levels)
names(cols) <- levels(levels_water$Levels)


#### SET DIRECTORIES ####
working.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
cm_dir <- paste(working.dir, "Connectivity_Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial Data", sep="/")

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

# Map of State Marine Parks
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

NTZ <- rbind(NTZ, AMP_NTZ) 

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

# Locations of the boat ramps - NEED TO FIX
BR <- st_read("Boat_Ramps.shp") %>% 
  st_transform(4283)%>%
  st_make_valid

#### PARAMETER VALUES ####

## Natural Mortality
# We have instantaneous mortality from Marriot et al 2011 and we need to convert that into monthly mortality
yearly_surv=exp(-0.146)
monthly_mort=1-(yearly_surv^(1/12))

M <- monthly_mort # Natural mortality rate per month

#Ricker recruitment model parameters (these are currently just made up values)
a <- 7
b <- 0.0017
M50 <- 2 # From Grandcourt et al. 2010
M95 <- 5 # From Grandcourt et al. 2010 (technically M100)

#Fish movement parameters
SwimSpeed <- 1.0 # Swim 1km in a day - this is completely made up 

## Fishing mortality parameters
A50 <- 4 # For L. miniatus from Williams et al. 2010 # L. miniatus becomes vulnerable to fishing at about age two
A95 <- 6 # For L. miniatus from Williams et al. 2012
q <- 0.5 # Apparently this is what lots of stock assessments set q to be...
#Biomass <- c(15, 0.55) # This is kg per 500m^2 (both inside and outside NTZ) from Comunity-Assisted Scientific Assessment and Management of WA MPAs, Ningaloo Reef Life Survey

#### SPATIAL DATA ####

## Create extent of area you want to cover 

ningaloo_map <- st_crop(wa_map, xmin=112.5, xmax=114.7, ymin=-24, ymax=-21.1) #NEED TO CHANGE THIS TO BE WHOLE GULF DOWN TO CORAL BAY
plot(ningaloo_map$geometry, col="#8ec3ca")
plot(BR, add=T, col="black", bg="#edcfcc", pch=21, cex=3)

# Make grid cells for fish to live in
grd <- st_make_grid(ningaloo_map, cellsize=0.05, square=FALSE)%>%
  st_crop(xmin=112.5, xmax=114.65, ymin=-24, ymax=-21.1) #make sure extent of grid is the same as the polygon
plot(grd, add=TRUE)

water <- st_difference(grd, ningaloo_map)
plot(water, border="#8ec3ca")
plot(BR$geometry)

# Adjust the sizes of the grid cells so that the ones closer to the shore are smaller 
GrdSmall <- st_make_grid(ningaloo_map, cellsize=0.025, square=FALSE)%>%
  st_crop(xmin=112.5, xmax=114.65, ymin=-24, ymax=-21.1) #make sure extent of grid is the same as the polygon
plot(GrdSmall)

HabitatSmall <- rbind(reef, rocky, lagoon)
HabitatSmall <- st_union(HabitatSmall) %>% 
  st_make_valid()
plot(HabitatSmall)

SmallGrd <- st_intersection(GrdSmall, HabitatSmall)
plot(SmallGrd)

# Merge the small grid with the grid with larger cells
BigGrd <- st_difference(water, HabitatSmall1)
plot(BigGrd)
plot(SmallGrd, add=T)

BigGrd <- st_make_valid(BigGrd)

water <- append(BigGrd, SmallGrd)
plot(water)

# Create a grid layer that is just for NTZs so we can adjust fishing mortality later
# Need some way to deal with the expansion of sanctuary zones over time once you get Cresswell's map
NTZgrd <- st_make_grid(NTZ, cellsize=0.05, square=FALSE) %>% 
  st_crop(xmin=112.5, xmax=114.2, ymin=-24, ymax=-21.1)
NTZarea <- st_intersection(NTZ, water) # Gives you just the hexagons in each NTZ with the comments
plot(NTZarea$geometry)

# Change the water layer so that it excludes the NTZs because we want to be able to differentiate them
NTZ <- st_union(NTZ)
plot(NTZ)

Fished_area <- st_difference(water, NTZ) 
plot(Fished_area$Spatial)


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

## Snap all the nodes together to eliminate 

# Get centroids for the grid cells 
centroids <- st_centroid_within_poly(water)
plot(centroids, cex=0.3) #Something is wrong with where the centroids are now

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

### CREATE MATRIX OF CONNVECTIVITY FOR FISH MOVEMENT ####
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
AverageYear <- boat_days %>% 
  group_by(Year) %>% 
  summarise(., Mean_Boat_Days=mean(Gascoyne_Boat_Days))
  
YearPlot <- ggplot(AverageYear) + 
  geom_point(aes(x=Year, y=Mean_Boat_Days))

YearModel <- lm(Gascoyne_Boat_Days~Year, data=boat_days)
summary(YearModel)

# We are currently predicting back to 1990 and then we'll use a different line from 1990 to 1960
Year2011_1990 <- as.data.frame(array(0, dim=c(21,1))) %>% 
  mutate(V1=seq(1990, 2010, by=1)) %>% 
  rename("Year"=V1)

predictions <- predict(YearModel, newdata=Year2011_1990)

Year2011_1990 <- Year2011_1990 %>% 
  mutate(Mean_Boat_Days = predictions)

boat_days_hind <- rbind(AverageYear, Year2011_1990)

# Creating a straight line from 1990 back to 1960
effort <- seq(0, 9872, length=30)
years <- seq(1960, 1989, by=1)

Years_1960_1989 <- as.data.frame(cbind(years, effort)) %>% 
  rename("Year" = years) %>% 
  rename("Mean_Boat_Days" = effort)

# Create the full fishing effort from 1960 to 2018
boat_days_hind <- rbind(boat_days_hind, Years_1960_1989)

YearPlot <- ggplot(boat_days_hind) + 
  geom_point(aes(x=Year, y=Mean_Boat_Days))

## Allocating the fishing effort in each year to months based on the proportions that we already have 
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
  group_by(Year) %>% 
  mutate(Gascoyne_Boat_Days = as.integer(unlist((Mean_Boat_Days*Month_Prop_Ave[,2])))) %>% 
  ungroup() %>% 
  dplyr::select(-c(Mean_Boat_Days))

check <- boat_days_hind %>% 
  group_by(Year) %>% 
  summarise(., sum(Gascoyne_Boat_Days))

# Put it all together into one big dataframe

Full_Boat_Days <- boat_days[,1:4] %>% 
  rbind(boat_days_hind) %>% 
  arrange(-Year)

# # Show this plot to Matt because things look hella janky and I'm not sure what we can do about it 
# MonthPlot <- Full_Boat_Days %>%
#   mutate(Unique = paste(Year, NumMonth, sep="_")) %>%
#   filter(Month %in% c("Jan","Jun")) %>%
#   ggplot() +
#   geom_point(aes(x=Unique, y=Gascoyne_Boat_Days))

## Next step is to allocate the effort to Exmouth and then to the specific boat ramps
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

Fishing <- Fishing*q

## Add NTZs
NoTake <- st_sf(water) %>% 
  st_drop_geometry() %>% 
  dplyr::select(Fished, ID, COMMENTS)

# Creating new columns for each of our year groups when new NTZs were put in place (easier than coding it in the model)
NoTake2 <- NoTake %>% 
  rename(Fished60_87="Fished") %>% 
  mutate(Fished87_05 = ifelse(str_detect(COMMENTS, c("Comm Cloates|Lighthouse Bay Sanctuary Zone| 
                                                     Gnaraloo Bay Sanctuary Zone|Cape Farquhar Sanctuary Zone")), "Y", "N")) %>% 
  mutate(Fished05_17 = ifelse(str_detect(COMMENTS, "Comm Cloates"), "Y", "N")) %>% 
  mutate(Fished17_21 = ifelse(str_detect(COMMENTS, "[:alpha:]"), "N", "Y")) %>% 
  dplyr::select(ID, Fished60_87, Fished87_05, Fished05_17, Fished17_21) %>% 
  mutate_at(vars(contains("Fished")), ~replace(., is.na(.), "Y"))


#### RUN MODEL ####

# Fill the array with random numbers of fish for now
# Need to change the loop so that rather than filling in a column it adds a new column each time 

YearlyTotal <- array(0, dim = c(NCELL,12,100)) #This is our yearly population split by age category (every layer is an age group)
for(d in 1:100){
  YearlyTotal[,1,] <- matrix(floor(runif(NCELL, 1, 2000)))
}

PopTotal <- array(0, dim=c(NCELL, 12, 59)) # This is our total population, all ages are summed and each column is a month (each layer is a year)

Movement <- NULL

# The time steps now represent 1 month (need to make sure you change the distance your fish can swim to match the timestep)
# NEED TO PUT IN NATURAL MORTALITY EACH MONTH AND ALSO THE FACT THAT THERE HASN'T ALWAYS BEEN NO TAKE ZONES

for(YEAR in 2:59){
  
  for(MONTH in 2:12){
    
    ### SUMMER 
    if(MONTH >=2 & MONTH<4|MONTH==12){
      ## Movement - this is where you'd change things to suit the months
      for(A in 2:dim(YearlyTotal)[3]){
        if(A==2 & MONTH==12){YearlyTotal[,12,A-1] <- matrix(floor(runif(NCELL, 1, 2000)))} # This is putting in random recruitment 
        else{ }
        
        Pop <- matrix(YearlyTotal[ ,MONTH-1,A-1]) #This gives you the fish in all the sites at timestep t-1 of age A-1
        Movers <- sum(Pop)
        Movement <- NULL
        Movement2 <- NULL
        
        for(s in 1:NCELL){
          Pop2 <- as.matrix(Pj[s,] * Pop[s,1]) #This should give you the number of fish that move from site s to all the other sites
          Pop2 <- t(Pop2)
          Movement <- rbind(Movement, Pop2) #This should give you an array with 143 rows representing the fish that move from each site to all the other sites 
        } #End bracket for movement in each cell
        
        Movement2 <- as.matrix(colSums(Movement))
        Moved <- sum(Movement2)
        print(isTRUE(all.equal(Movers,Moved)))
        if((isTRUE(all.equal(Movers,Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
        {print(Movers)
          print(Moved)}
        else{ }
        
        
        YearlyTotal[ ,MONTH,A-1] <- Movement2
        
      } #End bracket for movement in each month
      ## Fishing Mortality
      for(A in 1:dim(YearlyTotal)[3]){
        sa <- 1/(1+(exp(-log(19)*((A-A95/A95-A50))))) # This puts in size selectivity
        
        for(CELL in 1:NCELL){
          
          # Up until 1987 when there were no sanctuaries
          if(YEAR>0&YEAR<=26){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
          
          # From 1987 to 2005 when the first sanctuaries were in place
          else if(YEAR>=27&YEAR<=44){
            # If it's a NTZ then we don't have any fishing mortality in that cell
            if(NoTake[CELL, 3]=="Y"){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
            else { }  
          }
          # From 2005 to 2017 when we had all the state sanctuaries but no commonwealth one
          else if(YEAR>=45&YEAR<=57){
            # If it's a NTZ then we don't have any fishing mortality in that cell
            if(NoTake[CELL, 4]=="Y"){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
            else { }  
          }
          
          # From 2017 onwards when we had all the sanctuaries 
          else if(YEAR>57){
            # If it's a NTZ then we don't have any fishing mortality in that cell
            if(NoTake[CELL, 5]=="Y"){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
            else { }  
          }
        } # End fishing mortality for each cell
      } # End calculating selectivity for each age group
      
      # Natural mortality in each month
      YearlyTotal[,MONTH,] <- YearlyTotal[,MONTH,]*(1-M)
      
    } # SUMMER COMPLETED #
    
    ### AUTUMN
    else if(MONTH>=4 & MONTH<=5){
        ## Movement - this is where you'd change things to suit the months
        for(A in 2:dim(YearlyTotal)[3]){
          if(A==2 & MONTH==12){YearlyTotal[,12,A-1] <- matrix(floor(runif(NCELL, 1, 2000)))} # This is putting in random recruitment 
          else{ }
          
          Pop <- matrix(YearlyTotal[ ,MONTH-1,A-1]) #This gives you the fish in all the sites at timestep t-1 of age A-1
          Movers <- sum(Pop)
          Movement <- NULL
          Movement2 <- NULL
          
          for(s in 1:NCELL){
            Pop2 <- as.matrix(Pj[s,] * Pop[s,1]) #This should give you the number of fish that move from site s to all the other sites
            Pop2 <- t(Pop2)
            Movement <- rbind(Movement, Pop2) #This should give you an array with 143 rows representing the fish that move from each site to all the other sites 
          } #End bracket for movement in each cell
          
          Movement2 <- as.matrix(colSums(Movement))
          Moved <- sum(Movement2)
          print(isTRUE(all.equal(Movers,Moved)))
          if((isTRUE(all.equal(Movers,Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
          {print(Movers)
            print(Moved)}
          else{ }
          
          
          YearlyTotal[ ,MONTH,A-1] <- Movement2
          
        } #End bracket for movement in each month
      ## Fishing Mortality
      for(A in 1:dim(YearlyTotal)[3]){
        sa <- 1/(1+(exp(-log(19)*((A-A95/A95-A50))))) # This puts in size selectivity
        
        for(CELL in 1:NCELL){
          
          # Up until 1987 when there were no sanctuaries
          if(YEAR>0&YEAR<=26){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
          
          # From 1987 to 2005 when the first sanctuaries were in place
          else if(YEAR>=27&YEAR<=44){
            # If it's a NTZ then we don't have any fishing mortality in that cell
            if(NoTake[CELL, 3]=="Y"){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
            else { }  
          }
          # From 2005 to 2017 when we had all the state sanctuaries but no commonwealth one
          else if(YEAR>=45&YEAR<=57){
            # If it's a NTZ then we don't have any fishing mortality in that cell
            if(NoTake[CELL, 4]=="Y"){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
            else { }  
          }
          
          # From 2017 onwards when we had all the sanctuaries 
          else if(YEAR>57){
            # If it's a NTZ then we don't have any fishing mortality in that cell
            if(NoTake[CELL, 5]=="Y"){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
            else { }  
          }
        }
      } 
      
      # Natural mortality in each month
      YearlyTotal[,MONTH,] <- YearlyTotal[,MONTH,]*(1-M)
      
    } # AUTUMN COMPLETED #
    
    ### WINTER
    else if (MONTH >=6 & MONTH<=8){
        ## Movement - this is where you'd change things to suit the months
        for(A in 2:dim(YearlyTotal)[3]){
          if(A==2 & MONTH==12){YearlyTotal[,12,A-1] <- matrix(floor(runif(NCELL, 1, 2000)))} # This is putting in random recruitment 
          else{ }
          
          Pop <- matrix(YearlyTotal[ ,MONTH-1,A-1]) #This gives you the fish in all the sites at timestep t-1 of age A-1
          Movers <- sum(Pop)
          Movement <- NULL
          Movement2 <- NULL
          
          for(s in 1:NCELL){
            Pop2 <- as.matrix(Pj[s,] * Pop[s,1]) #This should give you the number of fish that move from site s to all the other sites
            Pop2 <- t(Pop2)
            Movement <- rbind(Movement, Pop2) #This should give you an array with 143 rows representing the fish that move from each site to all the other sites 
          } #End bracket for movement in each cell
          
          Movement2 <- as.matrix(colSums(Movement))
          Moved <- sum(Movement2)
          print(isTRUE(all.equal(Movers,Moved)))
          if((isTRUE(all.equal(Movers,Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
          {print(Movers)
            print(Moved)}
          else{ }
          
          
          YearlyTotal[ ,MONTH,A-1] <- Movement2
          
        }
      ## Fishing Mortality
      for(A in 1:dim(YearlyTotal)[3]){
        sa <- 1/(1+(exp(-log(19)*((A-A95/A95-A50))))) # This puts in size selectivity
        
        for(CELL in 1:NCELL){
          
          # Up until 1987 when there were no sanctuaries
          if(YEAR>0&YEAR<=26){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
          
          # From 1987 to 2005 when the first sanctuaries were in place
          else if(YEAR>=27&YEAR<=44){
            # If it's a NTZ then we don't have any fishing mortality in that cell
            if(NoTake[CELL, 3]=="Y"){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
            else { }  
          }
          # From 2005 to 2017 when we had all the state sanctuaries but no commonwealth one
          else if(YEAR>=45&YEAR<=57){
            # If it's a NTZ then we don't have any fishing mortality in that cell
            if(NoTake[CELL, 4]=="Y"){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
            else { }  
          }
          
          # From 2017 onwards when we had all the sanctuaries 
          else if(YEAR>57){
            # If it's a NTZ then we don't have any fishing mortality in that cell
            if(NoTake[CELL, 5]=="Y"){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
            else { }  
          }
        }
      } 
      
      # Natural mortality in each month
      YearlyTotal[,MONTH,] <- YearlyTotal[,MONTH,]*(1-M)
      
    } # WINTER COMPLETED #
    
    ### SPRING
    else if(MONTH>=9 & MONTH<=11){
        ## Movement - this is where you'd change things to suit the months
        for(A in 2:dim(YearlyTotal)[3]){
          if(A==2 & MONTH==12){YearlyTotal[,12,A-1] <- matrix(floor(runif(NCELL, 1, 2000)))} # This is putting in random recruitment 
          else{ }
          
          Pop <- matrix(YearlyTotal[ ,MONTH-1,A-1]) #This gives you the fish in all the sites at timestep t-1 of age A-1
          Movers <- sum(Pop)
          Movement <- NULL
          Movement2 <- NULL
          
          for(s in 1:NCELL){
            Pop2 <- as.matrix(Pj[s,] * Pop[s,1]) #This should give you the number of fish that move from site s to all the other sites
            Pop2 <- t(Pop2)
            Movement <- rbind(Movement, Pop2) #This should give you an array with 143 rows representing the fish that move from each site to all the other sites 
          } #End bracket for movement in each cell
          
          Movement2 <- as.matrix(colSums(Movement))
          Moved <- sum(Movement2)
          print(isTRUE(all.equal(Movers,Moved)))
          if((isTRUE(all.equal(Movers,Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
          {print(Movers)
            print(Moved)}
          else{ }
          
          
          YearlyTotal[ ,MONTH,A-1] <- Movement2
          
        } #End bracket for movement in each month
      ## Fishing Mortality
      for(A in 1:dim(YearlyTotal)[3]){
        sa <- 1/(1+(exp(-log(19)*((A-A95/A95-A50))))) # This puts in size selectivity
        
        for(CELL in 1:NCELL){
          
          # Up until 1987 when there were no sanctuaries
          if(YEAR>0&YEAR<=26){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
          
          # From 1987 to 2005 when the first sanctuaries were in place
          else if(YEAR>=27&YEAR<=44){
            # If it's a NTZ then we don't have any fishing mortality in that cell
            if(NoTake[CELL, 3]=="Y"){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
            else { }  
          }
          # From 2005 to 2017 when we had all the state sanctuaries but no commonwealth one
          else if(YEAR>=45&YEAR<=57){
            # If it's a NTZ then we don't have any fishing mortality in that cell
            if(NoTake[CELL, 4]=="Y"){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
            else { }  
          }
          
          # From 2017 onwards when we had all the sanctuaries 
          else if(YEAR>57){
            # If it's a NTZ then we don't have any fishing mortality in that cell
            if(NoTake[CELL, 5]=="Y"){YearlyTotal[CELL,MONTH,A] <- YearlyTotal[CELL,MONTH,A]*exp(-sa*Fishing[CELL,MONTH,YEAR])}
            else { }  
          }
        }
      } 
      
      # Natural mortality in each month
      YearlyTotal[,MONTH,] <- YearlyTotal[,MONTH,]*(1-M)
      
    } # SPRING COMPLETED #
    
    # ## Recruitment
    # Recs <- matrix(0, nrow=1, ncol=1) #Blank matrix to add the recruits to
    # Recruitment <- as.matrix(colSums(Pop[,t,], dims=1)) 
    # for(A in 1:dim(Recruitment)[1]){
    #   Mature <- 1/(1+(exp(-log(19)*((A-M95)/(M95-M50))))) #Number of mature individuals in each age class
    #   S <- colSums(Recruitment)*Mature #Spawning stock
    #   Rec <- a*S*exp(-b*S) #Number of recruits from that age class
    #   Recs <- rbind(Recs, Rec) #Combine recruits from each age class into one dataframe
    # }
    
  PopTotal[ , , YEAR] <- rowSums(YearlyTotal)
    
  }
  print(YEAR)
  water$pop <- rowSums(PopTotal[ , , YEAR])
  
  water <- water%>%
    mutate(Population = ifelse(pop < 1000, "<1000",
                               ifelse (pop>1000 & pop<110, "1000-1100",
                                       ifelse (pop>1100 & pop<1200, "1100-1200",
                                               ifelse (pop>1200 & pop<1300, "1200-1300",
                                                       ifelse(pop>1300 & pop<1400, "1300-1400",
                                                              ifelse(pop>1400 & pop<1500, "1400-1500",
                                                                     ifelse(pop>1500 & pop<1600, "1500-1600", ">1600"
                                                                            
                                                                     ))))))))%>%
    mutate(Population = factor(Population))
  water$Population <- fct_relevel(water$Population, "<1000",  "1000-1100", "1100-1200", "1200-1300",
                                  "1300-1400", "1400-1500", "1500-1600", ">1600")
  
  print(map <- ggplot(water)+
          geom_sf(aes(fill=Population, color=Fished))+
          scale_fill_manual(name="Population", values= cols, drop=FALSE)+
          scale_color_manual(name="Fished", values=c("white", "black"))+
          theme_void())
  Sys.sleep(3)
}


### Things I need to try and do
## Make recruitment equilibrate to natural mortality - need to figure out what numbers mean that the number of recruits 
# produced each year equals the number of fish that die every year
# Can also work out a and b using the equations in Eva's sea cucumber paper
## Also need to determine the size of the cells and the swimming speed so we can figure out how far the fish can move in
# one time step
# 50% kernal density has a mean of 2.3km and 95% is 8.6 
# In the lagoon they moved on average 2.92km in 30 days in the lagoon and 4.21km on the reef slope (Pillans 2014)
# Also show seasonal migration based on spawning
# Also move from the sargassum to the reefs as part of an ontogenic shift - need to have a different matrix for the one year old fish to make them move to the reef
## Will need to make sure when you put mortaility back in that you're subsetting the population between the ones that are in/outside a sanctuary zone and apply the correct mortality
# this will have to be based on the intersection of our grid with a shape file of the sanctuary zones

#### OLD THINGS THAT I'VE TRIED AND DON'T USE ANYMORE ####

##INITIAL POPULATION MODEL##
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

# 
# # For cells where we are missing data about habitat take raster data and then calculate
# # the average rugosity, for now we'll just say anything over XXX is reef
# setwd(sp_dir)
# bathy <- raster("Carnarvon_Shelf_Bathymetry_3_2008_3m_cog.tiff")
# bathy <- flip(bathy, direction="y")
# proj4string(bathy) #check projection
# 
# 
# 
# water_ras <- st_transform(water, "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# 
# grd_transform <- st_transform(grd, "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# 
# plot(bathy)
# plot(water_ras$water, add=TRUE)
# 
# st_crs(grd_transform)

# This needs to be habitat matrix*coefficient + distance*coefficient in some form or another. Need to figure out how to decide
# what the coefficients should be so we can weight the importance of habitat relative to distance

## Settlement
Settlement <- matrix(runif(143, 1, 10),nrow=143) # Random for now

# Need to normalise so that it sums to 1 so we're not adding or subtracting any fish anywhere
Settlement_norm <- (Settlement[,1] = (Settlement[,1]/sum(Settlement[,1])))
Settlement_norm <- as.matrix(Settlement_norm)

# ## Movement 
# Movement <- array(1, dim= c(143,143,1))
# for(c in 1:143){
#   Movement_norm[,c,] <- (Movement[,c,] = (Movement[,c,]/sum(Movement[,c,])))
# }
# Movement_norm <- as.matrix(Movement_norm)

## Set up the model
Pop <- array(0, dim = c(143,100,10)) #This is our population

# Fill the array with random numbers of fish for now
for(d in 1:10){
  Pop[,1,d] <- matrix(runif(143, 1, 10000))
}
Pop <- round(Pop, digits=0)
Recruits <- 2000


##Population Model###
## Run the model
for(t in 2:5){
  for(A in 1:dim(Pop)[3]){
    if(A==1){Pop[,t-1,A] <- Recruits}
    else {Pop[,t,A] <-  Pop[,t-1,A-1]}
    #else if (d>=6) {Days[,t,d] <-  Movement[,,1] %*% Days[,t-1,d-1]}
    #else {Days[,t,d] <- Movement[,,1] %*% Days[,t-1,d-1]}
  }
  
  
  # Mortality
  for(A in 1:dim(Pop)[3]){
    sa <- 1/(1+(exp(-log(19)*((A-A95/A95-A50)))))
    Pop[,t,A] <- Pop[,t,A]*exp(-M)*exp(-sa*Fishing)
  }
  
  
  # Ricker recruitment model with maturation
  Recs <- matrix(0, nrow=1, ncol=1) #Blank matrix to add the recruits to
  Recruitment <- as.matrix(colSums(Pop[,t,], dims=1)) 
  for(A in 1:dim(Recruitment)[1]){
    Mature <- 1/(1+(exp(-log(19)*((A-M95)/(M95-M50))))) #Number of mature individuals in each age class
    S <- colSums(Recruitment)*Mature #Spawning stock
    Rec <- a*S*exp(-b*S) #Number of recruits from that age class
    Recs <- rbind(Recs, Rec) #Combine recruits from each age class into one dataframe
  }
  R <- colSums(Recs) #Add up the number of recruits produced from all age classes
  Recruits <- as.matrix(Settlement_norm*R) #Distribute the recruits among the 143 sites
  
  # Plotting
  water$pop <- rowSums(Pop[,t,])
  
  water <- water%>%
    mutate(Population = ifelse(pop < 1000, "<1000",
                               ifelse (pop>1000 & pop<5000, "1000-5000",
                                       ifelse (pop>5000 & pop<10000, "5000-10000",
                                               ifelse (pop>10000 & pop<15000, "10000-15000",
                                                       ifelse(pop>15000 & pop<20000, "15000-20000",
                                                              ifelse(pop>20000 & pop<25000, "20000-25000",
                                                                     ifelse(pop>25000 & pop<30000, "25000-30000", ">30000"
                                                                            
                                                                     ))))))))%>%
    mutate(Population = factor(Population))
  water$Population <- fct_relevel(water$Population, "<1000",  "1000-5000", "5000-10000", "10000-15000",
                                  "15000-20000", "20000-25000", "25000-30000", ">30000")
  
  print(map <- ggplot(water)+
          geom_sf(aes(fill=Population))+
          scale_fill_manual(name="Population", values= cols, drop=FALSE)+
          theme_void())
  Sys.sleep(3)
  
}

## Making the fish move south
# Matrix that makes more southerly grid cells more preferable
latlong <- points[,1:2]
bearings <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (i in 1:NCELL){
  b <- bearing(latlong[i,], latlong)
  bearings[i, 1:NCELL] <- b
}

# We want all the bearings between 90-180 and -90 - -180 to be positive and everything else to be negative
# These are the bearings that indicate a southerly movement
bearings <- as.data.frame(bearings)%>%
  mutate_all(~ ifelse(.<(-90), .*-1,
                    ifelse(.>-90 & .<0, .*1,
                           ifelse(.>0 & .<90, .*-1, .))))

bearings <- as.matrix(bearings)

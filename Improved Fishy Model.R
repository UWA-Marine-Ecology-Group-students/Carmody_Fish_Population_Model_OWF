library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(stringr)
library(forcats)
library(RColorBrewer)
library(geosphere)
library(mgcv)


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
  mutate(Year = as.factor(Year)) %>% 
  mutate(NumMonth = as.numeric(NumMonth)) %>% 
  mutate(Month = as.factor(Month)) %>% 
  mutate(Gascoyne_Boat_Days = as.numeric(Gascoyne_Boat_Days)) %>% 
  mutate(Exmouth_Boat_Days = Gascoyne_Boat_Days*0.33) %>% 
  mutate(Month = fct_relevel(Month, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))

## Spatial Data
setwd(sp_dir)

# Map of WA Coastline
wa_map <- st_read("WACoastline.shp")%>%
  st_transform(4283)%>%
  st_make_valid

# Map of State Marine Parks
MP <- st_read("WA_MPA_2018.shp")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.3, ymin=-24, ymax=-21)

NTZ <- MP%>%
  filter(IUCN == "IA") %>% 
  dplyr::select(PA_ID, geometry)

# Habitat Layers
reef <- st_read("ReefHabitat.gpkg")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.3, ymin=-24, ymax=-21)

reef$type <- "reef"

lagoon <- st_read("LagoonHabitat.gpkg")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.3, ymin=-24, ymax=-21)

lagoon$type <- "lagoon"

rocky <- st_read("RockReefHabitat.gpkg")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.3, ymin=-24, ymax=-21)

rocky$type <- "rocky"

pelagic <- st_read("PelagicHabitat.gpkg")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.3, ymin=-24, ymax=-21)

pelagic$type <- "pelagic"

# Locations of the boat ramps
BR <- st_read("Boat_Ramps.shp") %>% 
  st_transform(4283)%>%
  st_make_valid

#### PARAMETER VALUES ####
M <- 0.146 # Natural mortality rate, Marriot et al. 2011

#Ricker recruitment model parameters (these are currently just made up values)
a <- 7
b <- 0.0017
M50 <- 2 # From Grandcourt et al. 2010
M95 <- 5 # From Grandcourt et al. 2010 (technically M100)

#Fish movement parameters
SwimSpeed <- 1.0 # Swim 1km in a day - this is completely made up 

## Fishing mortality parameters
E <- 25000/365 # This is in number of fishers out on the water per day, on average across a year for Ningaloo from 
               # A 12-month survey of recreational fishing in the Gascoyne bioregion of Western Australia during 1998-99 Summner 2002
A50 <- 4 # For L. miniatus from Williams et al. 2010
A95 <- 6 # For L. miniatus from Williams et al. 2012
q <- 0.5 # Apparently this is what lots of stock assessments set q to be...
#Biomass <- c(15, 0.55) # This is kg per 500m^2 (both inside and outside NTZ) from Comunity-Assisted Scientific Assessment and Management of WA MPAs, Ningaloo Reef Life Survey

#### SPATIAL DATA ####

## Create extent of area you want to cover 

ningaloo_map <- st_crop(wa_map, xmin=112.5, xmax=114.3, ymin=-24, ymax=-21)
plot(ningaloo_map$geometry)

# Make grid cells for fish to live in
grd <- st_make_grid(ningaloo_map, cellsize=0.05, square=FALSE)%>%
  st_crop(xmin=112.5, xmax=114.2, ymin=-24, ymax=-21.1) #make sure extent of grid is the same as the polygon
plot(grd, add=TRUE)

water <- st_difference(grd, ningaloo_map)
plot(water)

# Create a grid layer that is just for NTZs so we can adjust fishing mortality later
NTZgrd <- st_make_grid(NTZ, cellsize=0.05, square=FALSE) %>% 
  st_crop(xmin=112.5, xmax=114.2, ymin=-24, ymax=-21.1)
NTZarea <- st_intersection(NTZgrd, NTZ) # Gives you where the grid and the NTZs overlap
plot(NTZarea)

# Change the water layer so that it excludes the NTZs because we want to be able to differentiate them
NTZ <- st_union(NTZ)
plot(NTZ)

Fished_area <- st_difference(water, NTZ)
plot(Fished_area)

# turn areas into data frames for easier use
Fished_area <- st_sf(Fished_area) %>% 
  mutate(Fished = "Y") %>% 
  rename("Spatial" = Fished_area)

NTZarea <- st_sf(NTZarea) %>% 
  mutate(Fished = "N") %>% 
  rename("Spatial" = NTZarea)

water <- rbind(Fished_area, NTZarea)
plot(water, col=water$Fished)

water <- water%>%
  mutate(ID=row_number())%>%
  mutate(ID=factor(ID)) 

## Check that the NTZs are where you expect them to be
ggplot(water)+
  geom_sf(aes(fill=Fished))+
  theme_void()

# Get centroids for the grid cells 
centroids <- st_centroid_within_poly(water)
plot(centroids, cex=0.3 ,add=TRUE) #Something is wrong with where the centroids are now

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
# Need to put density dependence into the movement!

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


# Plotting the data to see what it looks like
Effort_Plot <- ggplot(boat_days)+
  geom_point(aes(Month, Exmouth_Boat_Days))+
  facet_wrap( ~ Year, strip.position = "bottom", scales = "free_x") +
  facet_grid(~Year, switch = "x", scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside")  




# Work out the probability of visiting a cell from each boat ramp
BR <- st_sf(BR)

water <- water %>% 
  mutate(DistBR = 0)

DistBR <- as.data.frame(array(0, dim = c(NCELL,3))) # This will contain the distance from each boat ramp to every cell

for(CELL in 1:NCELL){
  
  for(RAMP in 1:3){
    x <- st_distance(centroids[CELL,1], BR[RAMP,3])
    DistBR[CELL,RAMP] <- (x/1000) #to get the distance in km
    
  }
}



#### RUN MODEL ####

# Fill the array with random numbers of fish for now
# Need to change the loop so that rather than filling in a column it adds a new column each time 

PopTotal <- array(0, dim = c(NCELL,100,10)) #This is our population
for(d in 1:10){
  PopTotal[,1,] <- matrix(floor(runif(NCELL, 1, 2000)))
}

Movement <- NULL

# The time steps now represent 1 month (need to make sure you change the distance your fish can swim to match the timestep)
for(t in 2:12){

  ### SUMMER 
  if(t >=2 & t<4|t==12){
    ## Movement
    for(A in 2:dim(PopTotal)[3]){
      if(A==2){PopTotal[,t,A-1] <- matrix(floor(runif(NCELL, 1, 2000)))}
      else{ }
      
      Pop <- matrix(PopTotal[ ,t-1,A-1]) #This gives you the fish in all the sites at timestep t-1 of age A-1
      Movers <- sum(Pop)
      Movement <- NULL
      Movement2 <- NULL
      
      for(s in 1:NCELL){
        Pop2 <- as.matrix(Pj[s,] * Pop[s,1]) #This should give you the number of fish that move from site s to all the other sites
        Pop2 <- t(Pop2)
        Movement <- rbind(Movement, Pop2) #This should give you an array with 143 rows representing the fish that move from each site to all the other sites 
      }
      
      Movement2 <- as.matrix(colSums(Movement))
      Moved <- sum(Movement2)
      print(isTRUE(all.equal(Movers,Moved)))
      if((isTRUE(all.equal(Movers,Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
      {print(Movers)
        print(Moved)}
      else{ }
      
      
      PopTotal[ ,t,A] <- Movement2
      
    }
    
    ## Natural Mortality
    for(A in 1:dim(PopTotal)[3]){
      PopTotal[,t,A] <- PopTotal[,t,A]*exp(-M)
    }
    
    ## Fishing Mortality
    for(A in 1:dim(PopTotal)[3]){
      sa <- 1/(1+(exp(-log(19)*((A-A95/A95-A50))))) # This puts in size selectivity
      
      for(CELL in 1:NCELL){
        # If it's a NTZ then we don't have any fishing mortality in that cell
        if(Fishing[CELL, 2]=="Y"){PopTotal[CELL,t,A] <- PopTotal[CELL,t,A]*exp(-sa*Fishing[CELL,1])}
        else { }
        
      }
  } 
   

  }
  
  ### AUTUMN
  else if(t>=4 & t<=5){
    ## Movement
    for(A in 2:dim(PopTotal)[3]){
      if(A==2){PopTotal[,t,A-1] <- matrix(floor(runif(NCELL, 1, 2000)))}
      else{ }
      
      Pop <- matrix(PopTotal[ ,t-1,A-1]) #This gives you the fish in all the sites at timestep t-1 of age A-1
      Movers <- sum(Pop)
      Movement <- NULL
      Movement2 <- NULL
      
      for(s in 1:NCELL){
        Pop2 <- as.matrix(Pj[s,] * Pop[s,1]) #This should give you the number of fish that move from site s to all the other sites
        Pop2 <- t(Pop2)
        Movement <- rbind(Movement, Pop2) #This should give you an array with 143 rows representing the fish that move from each site to all the other sites 
      }
      
      Movement2 <- as.matrix(colSums(Movement))
      Moved <- sum(Movement2)
      print(isTRUE(all.equal(Movers,Moved)))
      if((isTRUE(all.equal(Movers,Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
      {print(Movers)
        print(Moved)}
      else{ }
      
      
      PopTotal[ ,t,A] <- Movement2
      
    }
    
    ## Natural Mortality
    for(A in 1:dim(PopTotal)[3]){
      PopTotal[,t,A] <- PopTotal[,t,A]*exp(-M)
    }
    
    ## Fishing Mortality
    for(A in 1:dim(PopTotal)[3]){
      sa <- 1/(1+(exp(-log(19)*((A-A95/A95-A50))))) # This puts in size selectivity
      
      for(CELL in 1:NCELL){
        # If it's a NTZ then we don't have any fishing mortality in that cell
        if(Fishing[CELL, 2]=="Y"){PopTotal[CELL,t,A] <- PopTotal[CELL,t,A]*exp(-sa*Fishing[CELL,1])}
        else { }
        
      }
    } 
  }
  
  ### WINTER
  else if (t >=6 & t<=8){
    ## Movement
    for(A in 2:dim(PopTotal)[3]){
      if(A==2){PopTotal[,t,A-1] <- matrix(floor(runif(NCELL, 1, 2000)))}
      else{ }
      
      Pop <- matrix(PopTotal[ ,t-1,A-1]) #This gives you the fish in all the sites at timestep t-1 of age A-1
      Movers <- sum(Pop)
      Movement <- NULL
      Movement2 <- NULL
      
      for(s in 1:NCELL){
        Pop2 <- as.matrix(Pj[s,] * Pop[s,1]) #This should give you the number of fish that move from site s to all the other sites
        Pop2 <- t(Pop2)
        Movement <- rbind(Movement, Pop2) #This should give you an array with 143 rows representing the fish that move from each site to all the other sites 
      }
      
      Movement2 <- as.matrix(colSums(Movement))
      Moved <- sum(Movement2)
      print(isTRUE(all.equal(Movers,Moved)))
      if((isTRUE(all.equal(Movers,Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
      {print(Movers)
        print(Moved)}
      else{ }
      
      
      PopTotal[ ,t,A] <- Movement2
      
    }
    
    ## Natural Mortality
    for(A in 1:dim(PopTotal)[3]){
      PopTotal[,t,A] <- PopTotal[,t,A]*exp(-M)
    }
    
    ## Fishing Mortality
    for(A in 1:dim(PopTotal)[3]){
      sa <- 1/(1+(exp(-log(19)*((A-A95/A95-A50))))) # This puts in size selectivity
      
      for(CELL in 1:NCELL){
        # If it's a NTZ then we don't have any fishing mortality in that cell
        if(Fishing[CELL, 2]=="Y"){PopTotal[CELL,t,A] <- PopTotal[CELL,t,A]*exp(-sa*Fishing[CELL,1])}
        else { }
        
      }
    } 
  }
  
  ### SPRING
  else if(t>=9 & t<=11){
    ## Movement
    for(A in 2:dim(PopTotal)[3]){
      if(A==2){PopTotal[,t,A-1] <- matrix(floor(runif(NCELL, 1, 2000)))}
      else{ }
      
      Pop <- matrix(PopTotal[ ,t-1,A-1]) #This gives you the fish in all the sites at timestep t-1 of age A-1
      Movers <- sum(Pop)
      Movement <- NULL
      Movement2 <- NULL
      
      for(s in 1:NCELL){
        Pop2 <- as.matrix(Pj[s,] * Pop[s,1]) #This should give you the number of fish that move from site s to all the other sites
        Pop2 <- t(Pop2)
        Movement <- rbind(Movement, Pop2) #This should give you an array with 143 rows representing the fish that move from each site to all the other sites 
      }
      
      Movement2 <- as.matrix(colSums(Movement))
      Moved <- sum(Movement2)
      print(isTRUE(all.equal(Movers,Moved)))
      if((isTRUE(all.equal(Movers,Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
      {print(Movers)
        print(Moved)}
      else{ }
      
      
      PopTotal[ ,t,A] <- Movement2
      
    }
    
    ## Natural Mortality
    for(A in 1:dim(PopTotal)[3]){
      PopTotal[,t,A] <- PopTotal[,t,A]*exp(-M)
    }
    
    ## Fishing Mortality
    for(A in 1:dim(PopTotal)[3]){
      sa <- 1/(1+(exp(-log(19)*((A-A95/A95-A50))))) # This puts in size selectivity
      
      for(CELL in 1:NCELL){
        # If it's a NTZ then we don't have any fishing mortality in that cell
        if(Fishing[CELL, 2]=="Y"){PopTotal[CELL,t,A] <- PopTotal[CELL,t,A]*exp(-sa*Fishing[CELL,1])}
        else { }
        
      }
    } 
  }
  
  # ## Recruitment
  # Recs <- matrix(0, nrow=1, ncol=1) #Blank matrix to add the recruits to
  # Recruitment <- as.matrix(colSums(Pop[,t,], dims=1)) 
  # for(A in 1:dim(Recruitment)[1]){
  #   Mature <- 1/(1+(exp(-log(19)*((A-M95)/(M95-M50))))) #Number of mature individuals in each age class
  #   S <- colSums(Recruitment)*Mature #Spawning stock
  #   Rec <- a*S*exp(-b*S) #Number of recruits from that age class
  #   Recs <- rbind(Recs, Rec) #Combine recruits from each age class into one dataframe
  # }
  
  
  print(t)
  water$pop <- rowSums(PopTotal[ ,t-1, ])
  
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

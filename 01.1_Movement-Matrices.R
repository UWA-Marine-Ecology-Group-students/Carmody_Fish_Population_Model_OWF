library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(stringr)
library(forcats)
library(RColorBrewer)
library(geosphere)

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

SwimSpeed <- 0.01

#### READ FILES ####
setwd(sp_dir)

dist_matrix <- readRDS("network_matrix")
pelagic_perc <- readRDS("pelagic_perc")
reef_perc <- readRDS("reef_perc")
lagoon_perc <- readRDS("lagoon_perc")
rocky_perc <- readRDS("rocky_perc")
water <- readRDS("water")

NCELL <- nrow(water)

water$cell_area <- st_area(water)

#### CREATE MATRIX OF CONNVECTIVITY FOR FISH MOVEMENT ####
## Assume the fish will try and swim the shortest path between locations
# Calculate the probability a fish moves to this site in a given time step using the swimming speed we set earlier 
# This creates a dispersal kernel based on the negative exponential distribution
# Some of these loops take quite a while to run

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

### Save all the files you need to remake these matrices ###

saveRDS(pDist, file="pDist")
saveRDS(pPelagic, file="pPelagic")
saveRDS(pReef, file="pReef")
saveRDS(pLagoon, file="pLagoon")
saveRDS(pRocky, file="pRocky")

#### CREATE PROBABILITY OF MOVEMENT USING UTILITY FUNCTION ####
## Load files if need be ##
setwd(sp_dir)
# 
# pDist <- readRDS("pDist")
pPelagic <- readRDS("pPelagic")
pReef <- readRDS("pReef")
pLagoon <- readRDS("pLagoon")
pRocky <- readRDS("pRocky")

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
dispersal$area <- as.vector(water$cell_area)

dispersal <- dispersal %>% 
  filter(perc_habitat!=0)

temp <- array(0, dim=c(NCELL, 1))
recruitment <- array(0, dim=c(nrow(dispersal), 2))
recruitment[ ,2] <- as.numeric(dispersal$ID) 

for (CELL in 1:NCELL){
  
  for (cell in 1:NCELL){
    temp[cell, 1] <- as.numeric((dispersal[cell,2]/dispersal[cell,1]))
  }
  
  temp[which(!is.finite(temp))] <- 0
  summed <- sum(temp)
  
  recruitment[CELL,1] <- as.numeric(temp[CELL,1]/summed) # This loop throws an error saying that the subscript is out of bounds but as far as I can tell it's not a problem, I just made the loop run for too long and haven't fixed it up yet because
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
saveRDS(Pj, file="movement")
saveRDS(water, file="water")
saveRDS(ProbRec, file="juvmove")
saveRDS(recruitment, file="recruitment")


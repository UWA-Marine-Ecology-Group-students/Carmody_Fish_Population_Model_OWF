###################################################

# Script for trying to figure out the MSY of the 
# population and then creating Kobe plots 

###################################################
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(forcats)
library(RColorBrewer)
library(MQMF)
library(Rcpp)
library(RcppArmadillo)
library(raster)
library(sfnetworks)
library(abind)

rm(list = ls())
#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")
pop_dir <-  paste(working.dir, "Output_Population", sep="/")
sim_dir <- paste(working.dir, "Simulations", sep="/")

model.name <- "ningaloo"

## Read in functions
setwd(working.dir)
sourceCpp("X_Model_RccpArm.cpp")
source("X_Functions.R")

#### READ IN DATA ####
setwd(sg_dir)
AdultMove <- readRDS(paste0(model.name, sep="_", "movement"))
Settlement <- readRDS(paste0(model.name, sep="_","recruitment")) 
Effort <- readRDS(paste0(model.name, sep="_", "fishing"))
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
water <- readRDS(paste0(model.name, sep="_","water"))
YearlyTotal <- readRDS(paste0(model.name, sep="_", "BurnInPop"))
Selectivity <- readRDS("selret")
Mature <- readRDS("maturity")
Weight <- readRDS("weight")


# Fishing effort surface
month.ave <- readRDS("Average_Monthly_Effort")

setwd(sp_dir)
BR <- st_read("Boat_Ramps.shp") %>% 
  st_transform(4283)%>%
  st_make_valid()

network <- st_read(paste0(model.name, sep="_","network.shapefile.shp"))

NCELL <- nrow(water)

#### RUN MSY MODEL ####
#* PARAMETER VALUES ####
## Natural Mortality
# We have instantaneous mortality from Marriott et al 2011 and we need to convert that into monthly mortality
NatMort <- 0.146
step <- 1/12 # We're doing a monthly time step here

# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
BHa = 0.4344209 #0.4344209
BHb = 0.0009398152 #0.01889882
PF = 0.5

# PlotTotal <- T #This is whether you want a line plot of the total or the map

Pop.Groups <- seq(1,12)

#### SET UP FISHING EFFORT FOR EACH LEVEL OF F ####

q <- 0.01 # Value often used when we don't know what our value of q is meant to be

## Effort values 
F_finite_values <- seq(0, 0.9, 0.05) # These are the values of F that we want to cycle through to see where our MSY is
E_values <- as.data.frame(array(0, dim = c(19, 13))) %>% 
  mutate(V1 = -log(1-F_finite_values)) %>%  # These are my big F instantaneous values
  mutate(V1 = V1/q)

names(E_values)[1:13] <- c("Yearly_Total", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

## Allocation to months
Monthly_effort <- E_values[,2:13]

for(E in 1:19){
  for(M in 1:12){
    Monthly_effort[E,M] <- month.ave[M,2] * E_values[E,1] 
  }
}

Monthly_effort <- Monthly_effort %>% 
  mutate(FM = F_finite_values)

## Allocation to boat ramps
BR_Trips <- data.frame(BoatRamp = c("Tantabiddi", "Bundegi", "ExmouthMar", "CoralBay"),
                       Effort = c(194.48, 133.9, 166.15, 146.07),
                       Trips = c(224, 157, 198, 151))

BR_Trips <- BR_Trips %>% 
  mutate(Trip_per_Hr = as.numeric(unlist((Trips/Effort)))) %>% # Standardise the no. trips based on how much time you spent sampling
  mutate(BR_Prop = Trip_per_Hr/sum(Trip_per_Hr)) #Then work out the proportion of trips each hour that leave from each boat ramp

## Format the data so that you can work out per month and per boat ramp
Full_Effort <- Monthly_effort %>% 
  pivot_longer(!FM, names_to = "Month", values_to = "Effort")

## Combine with boat ramp data
for(fm in 1:132){
  Full_Effort$Tb_BR = Full_Effort$Effort*BR_Trips[1,5]
  Full_Effort$Bd_BR = Full_Effort$Effort*BR_Trips[2,5]
  Full_Effort$ExM_BR = Full_Effort$Effort*BR_Trips[3,5] 
  Full_Effort$CrB_BR = Full_Effort$Effort*BR_Trips[4,5]
}

check <- Full_Effort %>% 
  mutate(Total = Tb_BR+Bd_BR+ExM_BR+CrB_BR) #Looks good!

## Allocating effort to the cells

## Work out the probability of visiting a cell from each boat ramp based on distance and size
## Work out the probability of visiting a cell from each boat ramp based on distance and size
BR <- st_as_sf(BR)
st_crs(BR) <- NA 

BR <- BR[1:4, ] %>% 
  mutate(name = c("Bundegi","Exmouth","Tantabiddi","CoralBay"))

centroids <- st_centroid_within_poly(water)
points <- as.data.frame(st_coordinates(centroids))%>% #The points start at the bottom left and then work their way their way right
  mutate(ID=row_number()) 
points_sf <- st_as_sf(points, coords = c("X", "Y")) 
st_crs(points_sf) <- NA


network <- as_sfnetwork(network, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

net <- activate(network, "nodes")
network_matrix <- st_network_cost(net, from=BR, to=points_sf)
network_matrix <- network_matrix*111
dim(network_matrix)

DistBR <- as.data.frame(t(network_matrix)) %>% 
  rename("Bd_BR"=V1) %>% 
  rename("ExM_BR" = V2) %>% 
  rename("Tb_BR" = V3) %>% 
  rename("CrB_BR"=V4)

# Create a data frame with both the distances and the areas of the cells
Cell_Vars <- DistBR %>% 
  mutate(Area = as.vector((water$cell_area)/100000)) #Cells are now in km^2 but with no units

## Work out the utilities for each cell
BR_U <- as.data.frame(matrix(0, nrow=NCELL, ncol=4)) #Set up data frame to hold utilities of cells

BR_U <- BR_U %>% #make sure you give the columns good names so that you know what they are
  rename("Bd_U"=V1) %>% 
  rename("ExM_U" = V2) %>% 
  rename("Tb_U" = V3) %>% 
  rename("CrB_U"=V4)

# Add coefficients to each variable - all the BRs are the same but make them negative because the further away they are the less likely people are to visit
a = 0.1
b = -0.1 
c = -0.1 
d = -0.1
e = -0.1

Vj <- Cell_Vars %>% 
  mutate(Bd_BR = Bd_BR*b,
         ExM_BR = ExM_BR*c,
         Tb_BR = Tb_BR*d,
         CrB_BR = CrB_BR*e,
         Area= Area*a)
  
rowU <- matrix(NA, ncol=1, nrow=NCELL)
cellU <- matrix(NA, ncol=4, nrow=NCELL)

for(RAMP in 1:4){
  for(cell in 1:NCELL){
    U <- exp(Vj[cell,5]/Vj[cell,RAMP])
    cellU[cell, RAMP] <- U
  }
}

rowU <- as.data.frame(colSums(cellU))

Pj <- matrix(NA, ncol=4, nrow=NCELL)

for (RAMP in 1:4){
  for (cell in 1:NCELL){
    Pj[cell,RAMP] <- (exp(Vj[cell,5]/Vj[cell,RAMP]))/rowU[RAMP,1]
  }
}
colSums(Pj)

# Allocate effort to cells
Fishing_MSY <- array(0, dim=c(NCELL, 12,19)) #This array has a row for every cell, a column for every month, and a layer for every value of F
Months <- array(0, dim=c(NCELL, 12))
Ramps <- array(0, dim=c(NCELL, 4))
layer <- 1

for(fm in 1:19){
  
  for(MONTH in 1:12){
    
    for(RAMP in 1:4){
      
      temp <- Full_Effort %>% 
        filter(FM == as.numeric(F_finite_values[fm])) %>% 
        dplyr::select(-c(Month, FM, Effort))
      
      temp <- as.matrix(temp) 
      
      for(CELL in 1:NCELL){
        Ramps[CELL,RAMP] <- Pj[CELL,RAMP] * temp[MONTH,RAMP]
      }
    }
    
    Months[,MONTH] <-  rowSums(Ramps)
  }
  Fishing_MSY[ , ,layer] <- Months 
  layer <- layer+1
}

Fishing_MSY <- Fishing_MSY*q
Effort <- Fishing_MSY

sum(Effort[,,8])

#### SET UP INITIAL POPULATION ####
NCELL <- nrow(water)
Ages <- seq(1,30) #These are the ages you want to plot 
MaxAge <- 30
MaxYear <- 59
MaxCell <- NCELL
PF <- 0.5

survived.age <- array(0, dim=c(MaxAge,19))
YPR <- array(0, dim=c(MaxAge, 1))
Biomass <- array(0, dim=c(MaxAge, 1))
SB <- array(0, dim=c(MaxAge, 1))
YPR.F <- array(0, dim=c(19, 2))
YPR.F[,1] <- seq(0,0.9,0.05)
Biomass.F <- array(0, dim=c(19, 2))
Biomass.F[,1] <- seq(0,0.9,0.05)
SB.F <- array(0, dim=c(19, 2))
SB.F[,1] <- seq(0,0.9,0.05)

bio.catch <- array(0, dim=c(MaxAge, 19))
yearly.catch <- array(0, dim=(c(MaxYear, 1)))
catch.by.age <- array(0, dim=(c(MaxAge, 19)))
Flevel.catch <- array(0, dim=c(19, 1))

Selectivity <- Selectivity[,,45:59] # Have to modify this to just be years where selectivity is the same as it is now, otherwise we get the selectivity where everything gets fished

##### RUN MODEL #####
# Need to get it to produce the plots we want as well 
# Want it to run for a year and then get the values for the population at the end of the year 
setwd(sg_dir)

for(FM in 1:19){
  
  print(FM)
  
  start.pop <- readRDS(paste0(model.name, sep="_","Starting_Pop"))
  Total <- array(0, dim=c(MaxYear,1))
  
  YearlyTotal <- array(0, dim = c(MaxCell,12,30)) #This is our yearly population split by age category (every layer is an age group)
  
  start.pop.year <- start.pop %>% 
    slice(which(row_number() %% 12 == 1)) # Gives you the total in each age group at the end of the year
  
  for(d in 1:dim(YearlyTotal)[3]){ # This allocates fish to cells randomly with the fish in age group summing to the total we calculated above - beware the numbers change slightly due to rounding errors
    for(N in 1:start.pop.year[d,1]){
      cellID <- ceiling((runif(n=1, min = 0, max = 1))*MaxCell)
      YearlyTotal[cellID,1,d] <- YearlyTotal[cellID,1,d]+1
    }
  }
  
  
  Effort <- abind(Fishing_MSY[,,FM],Fishing_MSY[,,FM], along=3)
  Effort <- abind(Effort, Effort, along=3)
  Effort <- abind(Effort, Effort, along=3)
  Effort <- abind(Effort, Effort, along=3)
  
  
  for (YEAR in 1:10){ # Start of model year loop - set it as 45 to make sure the selectivity is the same as present day
    
    ## Loop over all the Rcpp functions in the model
    
    ModelOutput <- RunModelfunc_cpp(YEAR, MaxAge, MaxYear, MaxCell, NatMort, BHa, BHb, PF, 
                                    AdultMove, Mature, Weight, Settlement, 
                                    YearlyTotal, Selectivity, Effort)
    
    survived.age[,FM] <- colSums(ModelOutput$YearlyTotal[,12,]) #No. fish that survived the year in each age class
    
    
    total.catch <- ModelOutput$Month_catch # This gives you catch by *weight* in each cell for each month, which each layer representing an age class
    bio.catch[,FM] <- colSums(total.catch[,,1:MaxAge], dim=2) # Biomass of fish caught in each age group 
    monthly.catch <- ModelOutput$age_catch #This is the *number* of fish in each age class caught in each cell
    catch.by.age[,FM] <- colSums(monthly.catch[,,1:MaxAge], dim=2) 
    
    
  } # End of model year loop
  
  for(A in 1:30){
    YPR[A,1] <- ypr.func(survived.age[,FM], F_finite_values, M, Weight, Age=A, Selectivity)
    YPR.F[FM,2] <- sum(YPR)
  }      
      
      Biomass <- bio.func(survived.age[,FM], Weight, YearlyTotal)
      
      ## Currently not right because you're using month 12 but you need to change how you get survived.age to get month 10 as well
      SB <- SB.func(survived.age[,FM], Weight, YearlyTotal, Mature)
      
      
      
      Biomass.F[FM,2] <- sum(Biomass)
      SB.F[FM,2] <- sum(SB)

    

  Flevel.catch[FM,1] <- sum(monthly.catch) # Total catch across the year for each level of F
  
  # Create plots at the very end of the loop with all of the different values of F
    MSY.plots <- msy.plot.func(YPR.F, Biomass.F, SB.F) # List with my three different plots in it
}

MSY.plots[[1]]
MSY.plots[[2]]
MSY.plots[[3]]

dead <- ModelOutput$age_died
dead30 <- sum(dead[,,30])
survived <- ModelOutput$age_survived
survived30 <- sum(survived[,,30])
catch <- ModelOutput$age_catch
catch30 <- sum(catch[,,30])


#### KOBE PLOT ####
F_SB <- cbind(TotMatBio, YearlyFishing) %>% 
  mutate(Year = c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018")) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) 

F_SB_Plot <- F_SB %>% 
  mutate(Group = 1) %>% 
  ggplot(., aes(x=TotalBio, y=YearlyFishing, label=Year, color=NTGroup, group=Group))+
  geom_point()+
  geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_path()+
  scale_x_continuous(breaks=seq(0,12,1), limits=c(0,12))+
  theme_classic()+
  xlab("Total Spawning Biomass")+
  ylab("Yearly Fishing Effort")



##### SORTING THIS OUT #####
tot_survived_monthly <- array(0, dim=c(30,12))
tot_survived <- array(0, dim=c(30,19))

MONTH <- 1
YEAR <- 1
FM <- 19

start.pop <- readRDS(paste0(model.name, sep="_","Starting_Pop"))
Total <- array(0, dim=c(MaxYear,1))

YearlyTotal <- array(0, dim = c(MaxCell,12,30)) #This is our yearly population split by age category (every layer is an age group)

start.pop.year <- start.pop %>% 
  slice(which(row_number() %% 12 == 1)) # Gives you the total in each age group at the end of the year

for(d in 1:dim(YearlyTotal)[3]){ # This allocates fish to cells randomly with the fish in age group summing to the total we calculated above - beware the numbers change slightly due to rounding errors
  for(N in 1:start.pop.year[d,1]){
    cellID <- ceiling((runif(n=1, min = 0, max = 1))*MaxCell)
    YearlyTotal[cellID,1,d] <- YearlyTotal[cellID,1,d]+1
  }
}

Month_effort_cehck <- Monthly_effort*q

for(FM in 1:19){
  for(AGE in 1:30){
    tot_survived[AGE, 2] = sum(exp(-(Selectivity[15, ,45]*Effort[, ,45])))
  }
}


temp <- Selectivity[15, ,45]*Effort[, ,45]

sum(temp)

MONTH <- 6
YEAR <- 58

for(AGE in 0:29){
  ModelOutput <- mortalityfunc_cpp(AGE, MaxCell, MONTH, YEAR, NatMort,
                                   Weight, Selectivity, YearlyTotal, Effort)
}

temp <- array(0, dim=c(MaxCell,1))
for(CELL in 1:MaxCell){
  temp[CELL,1] <- 1 - (exp(-Selectivity[AGE,6,40]*Effort[CELL,6,40]))
}

sum(PopTotal[,6,59])
sum(ModelOutput$tot_survived)
sum(ModelOutput$tot_died)
sum(ModelOutput$Fish_catch)
sum(ModelOutput$finite_f)

survived <- sum(ModelOutput$tot_survived) #14046.86
died <- sum(ModelOutput$tot_died) #173.1827
start <- survived + died #14220.04
caught <- sum(ModelOutput$Fish_catch) #0.0177405 but should be = to died.from.fish which is 1.235338
natural <- start *(1-exp(-(NatMort/12.0))) #171.9623
died.from.fish <- sum(ModelOutput$Pop_check)




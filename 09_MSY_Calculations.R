library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(forcats)
library(RColorBrewer)
library(MQMF)

rm(list = ls())
#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")
pop_dir <-  paste(working.dir, "Output_Population", sep="/")
sim_dir <- paste(working.dir, "Simulations", sep="/")

model.name <- "small"

#### READ IN DATA ####
setwd(sg_dir)

NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
water <- readRDS(paste0(model.name, sep="_","water"))
StartPop <- readRDS(paste0(model.name, sep="_", "BurnInPop"))
selectivity <- readRDS("selret")
maturity <- readRDS("maturity")
weight <- readRDS("weight")

# Fishing effort surface
month.ave <- readRDS("Average_Monthly_Effort")

setwd(sp_dir)
BR <- st_read("Boat_Ramps.shp") %>% 
  st_transform(4283)%>%
  st_make_valid()

## Read in functions
setwd(working.dir)
source("X_Functions.R")

NCELL <- nrow(water)

#### RUN MSY MODEL ####
#* PARAMETER VALUES ####
## Natural Mortality
# We have instantaneous mortality from Marriott et al 2011 and we need to convert that into monthly mortality
M <- 0.146
step <- 1/12 # We're doing a monthly time step here

# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
alpha <- 0.4344209 #0.4344209
beta <-	0.01889882 #0.01889882

NCELL <- nrow(water)
Ages <- seq(1,30) #These are the ages you want to plot 
Time <- seq(1,59) #This is how long you want the model to run for
# PlotTotal <- T #This is whether you want a line plot of the total or the map

Pop.Groups <- seq(1,12)

#### SET UP FISHING EFFORT FOR EACH LEVEL OF F ####

q <- 0.005 # Value often used when we don't know what our value of q is meant to be

## Effort values 
F_values <- seq(0, 0.5, 0.05) # These are the values of F that we want to cycle through to see where our MSY is
E_values <- as.data.frame(array(0, dim = c(11, 13))) %>% 
  mutate(V1 = F_values/q) # These are my yearly effort values that come from catchability and the level F we want 

names(E_values)[1:13] <- c("Yearly_Total", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

## Allocation to months
Monthly_effort <- E_values[,2:13]

for(E in 1:11){
  for(M in 1:12){
    Monthly_effort[E,M] <- month.ave[M,2] * E_values[E,1] 
  }
}

Monthly_effort <- Monthly_effort %>% 
  mutate(FM = F_values)

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
BR <- st_sf(BR)

water <- water %>% 
  mutate(DistBR = 0) %>% 
  mutate(cell_area = st_area(Spatial))

centroids <- st_centroid_within_poly(water)

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
  mutate(Area = as.vector((water$cell_area)/100000)) %>% #Cells are now in km^2 but with no units
  rename("Bd_BR"=V1) %>% 
  rename("ExM_BR" = V2) %>% 
  rename("Tb_BR" = V3) %>% 
  rename("CrB_BR"=V4)

## Work out the utilities for each cell
BR_U <- as.data.frame(matrix(0, nrow=NCELL, ncol=4)) #Set up data frame to hold utilities of cells

BR_U <- BR_U %>% #make sure you give the columns good names so that you know what they are
  rename("Bd_U"=V1) %>% 
  rename("ExM_U" = V2) %>% 
  rename("Tb_U" = V3) %>% 
  rename("CrB_U"=V4)

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
         Area= Area*a) %>% 
  mutate(vj = Bd_BR+ExM_BR+Tb_BR+CrB_BR+Area)
  
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
Fishing_MSY <- array(0, dim=c(NCELL, 12,11)) #This array has a row for every cell, a column for every month, and a layer for every value of F
Months <- array(0, dim=c(NCELL, 12))
Ramps <- array(0, dim=c(NCELL, 4))
layer <- 1

for(fm in 1:11){
  
  for(MONTH in 1:12){
    
    for(RAMP in 1:4){
      
      temp <- Full_Effort %>% 
        filter(FM == as.numeric(F_values[fm])) %>% 
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

#### SET UP INITIAL POPULATION ####
YearlyTotal <- StartPop
survived.age <- array(0, dim=c(length(Ages),12))
YPR <- array(0, dim=c(length(Ages), 1))
Biomass <- array(0, dim=c(length(Ages), 1))
SB <- array(0, dim=c(length(Ages), 1))
YPR.F <- array(0, dim=c(11, 2))
YPR.F[,1] <- seq(0,0.5,0.05)
Biomass.F <- array(0, dim=c(11, 2))
Biomass.F[,1] <- seq(0,0.5,0.05)
SB.F <- array(0, dim=c(11, 2))
SB.F[,1] <- seq(0,0.5,0.05)

bio.catch <- array(0, dim=c(NCELL, length(Ages)))
monthly.catch <- array(0, dim=c(length(Ages), 12))
Flevel.catch <- array(0, dim=c(11, 1))

NCELL <- nrow(water)
Ages <- seq(1,30) #These are the ages you want to plot 

##### RUN MODEL #####
# Need a loop that iterates over the values of F (start with 0 - 0.5 in 0.05 increments)
# Need to get it to produce the plots we want as well 
# I think the only function I'll be using from the main model is the mortality function
for(FM in 1:11){
  for(MONTH in 1:12){
    for(A in 1:dim(StartPop)[3]){
      
      if(MONTH==12 & 2<=A & A<30){
        survived.catch <- mortality.func(Age=A, Nat.Mort=M, Effort=Fishing_MSY, Max.Cell = NCELL,
                                         Month=MONTH, Select=selectivity, Population=StartPop, Year=FM)
        
        YearlyTotal[ ,1, A+1] <- survived.catch[[1]]
        
        
        # Calculate catch
        n.catch <- survived.catch[[2]]
        
        bio.catch[ ,A] <- n.catch * weight[(A*12)+1]
      } else if (MONTH!=12 & A>0) {
        survived.catch <- mortality.func(Age=A, Nat.Mort=M, Effort=Fishing_MSY, Max.Cell = NCELL,
                                         Month=MONTH, Select=selectivity, Population=YearlyTotal, Year=FM) # Changing Year to be the level of fishing mortality
        
        YearlyTotal[ ,MONTH+1,A] <- survived.catch[[1]]
        
        # Calculate catch
        n.catch <- survived.catch[[2]]
        
        bio.catch[ ,A] <- n.catch * weight[(A*12)+1] # Catch of each age class in each month
      }
     
       survived.age[A,MONTH] <- sum(survived.catch[[1]]) # Gives us one value for all fish of that age group that survived
    
       }
    # This is where we just take the values last month of the year
    if(MONTH==12){
      for(A in 1:dim(YearlyTotal)[3])
        
      YPR <- ypr.func(survived.age, F_values, M, weight, Age=A, selectivity)
      #Biomass[A] <- bio.func(blah)
      #SB[A] <- SB.func(blah)
      
      
      YPR.F[FM,2] <- sum(YPR)
      #Biomass.F[FM,2] <- sum(Biomass[,12])
      #SB.F[FM,2] <- sum(SB[,12])
    } else { }
    
    monthly.catch[1:30,MONTH] <- colSums(bio.catch) # Total catch in each month
    
  }
  Flevel.catch[FM,1] <- sum(monthly.catch) # Total catch across the year for each level of F
  
  # Create plots at the very end of the loop with all of the different values of F
  if(A==30 & MONTH==12){
    MSY.plots <- msy.plot.func(YPR.F, Biomass.F, SB.F) # List with my three different plots in it
  } else { }
}
MSY.plots[[1]]



#### WORKING OUT BIOMASS ####
Years <- seq(5, 55, 5)
Years[12] <- 59

TotMatBio <- matrix(0, ncol=12, nrow=30)

for(Y in 1:12){
  
  Population <-  get(paste0("Year", Years[Y]))
  
  for (A in 1:30){
    adults <- Population[ ,10, ] %>% 
      colSums(.) # Gives us just females because they are the limiting factor for reproduction
    adults <- adults * 0.5
    
    MatBio<- lapply(1:dim(Population)[3], function(Age){
      SB <- adults[Age] * maturity[Age,1] * weight[(Age*12)+1] #Gives us spawning biomass in each age group at the end of the year, hence the x 12+1 as it starts at 1 not zero
      TotMatBio <- sum(SB) #Gives us total mature spawning biomass
    })
    MatBio <- do.call(rbind, MatBio)
    TotMatBio[ ,Y] <- MatBio
  }
}

TotMatBio <- as.data.frame(TotMatBio) %>% 
  `colnames<-`(c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018")) %>% 
  mutate(Age = seq(1, 30, 1)) %>% 
  pivot_longer(cols=-c(Age), names_to="Year", values_to="Number") %>% 
  mutate(Year = as.factor(Year)) %>% 
  group_by(Year) %>% 
  summarise(TotalBio=sum(Number)) %>% 
  mutate(TotalBio = TotalBio/1000)

YearlyFishing <- fishing %>% 
  colSums(fishing) %>% 
  colSums(.)

YearlyFishing <- YearlyFishing[Years, ]
YearlyFishing <- as.data.frame(YearlyFishing)

YearlyFishing <- YearlyFishing*q


#### PLOTTING ####
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

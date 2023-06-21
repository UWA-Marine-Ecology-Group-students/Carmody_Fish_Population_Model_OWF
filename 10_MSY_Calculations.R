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
Water <- readRDS(paste0(model.name, sep="_","water"))
YearlyTotal <- readRDS(paste0(model.name, sep="_", "BurnInPop"))
Selectivity <- readRDS("selret")
Mature <- readRDS("maturity")
Weight <- readRDS("weight")
month.ave <- readRDS("Average_Monthly_Effort")


NCELL <- nrow(Water)

#### SET UP FOR NTZ AND FISHED AREA ####
setwd(sp_dir)
bathy <- raster("ga_bathy_ningaloocrop.tif")
WHA <- st_read("2013_02_WorldHeritageMarineProgramme.shp") %>% 
  st_transform(4283)%>%
  st_make_valid %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) 


#* Create list of cells to restrict plots to shallow Water (<30m)
Water_points <- st_centroid_within_poly(Water) 

Water_bathy <- raster::extract(bathy, Water_points, fun=mean, df=TRUE)

Water_bathy <- Water_bathy %>% 
  mutate(ID = as.factor(ID))

model_WHA <- Water %>% 
  st_intersects(., WHA) %>% 
  as.data.frame()

Water_WHA <- Water[c(as.numeric(model_WHA$row.id)), ]

Water_shallow <- Water_WHA %>% 
  mutate(ID = as.factor(ID)) %>% 
  left_join(., Water_bathy, by="ID") %>% 
  rename(bathy = "ga_bathy_ningaloocrop") %>% 
  filter(bathy >= c(-30)) %>% 
  filter(!is.na(bathy))

shallow_cells_NTZ <- Water_shallow %>% 
  filter(Fished_2017 %in% c("N")) %>% 
  st_drop_geometry(.)

shallow_cells_F <- Water_shallow %>% 
  filter(Fished_2017 %in% c("Y")) %>% 
  st_drop_geometry(.)

shallow_NTZ_ID <- as.numeric(levels(shallow_cells_NTZ$ID))[as.integer(shallow_cells_NTZ$ID)]
shallow_F_ID <- as.numeric(levels(shallow_cells_F$ID))[as.integer(shallow_cells_F$ID)]

NCELL <- length(shallow_F_ID)

#### SET UP FISHING EFFORT FOR EACH LEVEL OF F ####

## Effort values 

# F_finite_values <- seq(0, 0.9, 0.05) # These are the values of F that we want to cycle through to see where our MSY is
F_finite_values <- seq(0, 0.15, 0.01)
E_values <- as.data.frame(array(0, dim = c(16, 13))) %>% 
  mutate(V1 = F_finite_values) 
  # mutate(V1 = V1/q)

names(E_values)[1:13] <- c("Yearly_Total", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

## Allocation to months
Monthly_effort <- E_values[,2:13]

for(E in 1:16){
  for(M in 1:12){
    Monthly_effort[E,M] <- month.ave[M,2] * E_values[E,1] 
  }
}

Monthly_effort <- Monthly_effort %>% 
  mutate(FM = F_finite_values)

# Allocate effort to cells
NCELL <- length(shallow_NTZ_ID) + length(shallow_F_ID)

Fishing_MSY <- array(0, dim=c(NCELL, 12,16)) #This array has a row for every cell, a column for every month, and a layer for every value of F

Effort.List <- list()

for(FM in 1:16){
  
  for(M in 1:12){
    
    Fishing_MSY[,M,FM] <- Monthly_effort[FM,M]
  }
  
  Effort <- abind(Fishing_MSY[,,FM],Fishing_MSY[,,FM], along=3)
  Effort <- abind(Effort, Effort, along=3)
  Effort <- abind(Effort, Effort, along=3)
  Effort <- abind(Effort, Effort, along=3)
  Effort <- abind(Effort, Effort, along=3)
  Effort <- abind(Effort, Effort, along=3)
  
  Effort.List[[FM]] <- Effort
  
}

#### PARAMETER VALUES ####
MaxAge = 30
MaxYear = 50
MaxCell = length(shallow_NTZ_ID) + length(shallow_F_ID)
PF = 0.5 # Proportion female

## Natural Mortality
# We have instantaneous mortality from Marriott et al 2011 and we need to convert that into monthly mortality
NatMort = 0.146
step = 1/12 # We're doing a monthly time step here

# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
BHa = 0.4344209 
BHb = 0.003759261 

# Collect catch data
age.catch <- array(0, dim=c(12, MaxAge, MaxYear))
catch.by.age <- array(0, dim=c(MaxAge, MaxYear))
FM.Age.Catches.All <- list()
FM.Age.Catches.NTZ <- list()
FM.Age.Catches.F <- list()

catch.by.weight <- array(0, dim=c(MaxCell, MaxYear))
FM.Weight.Catches.All <- list()
FM.Weight.Catches.NTZ <- list()
FM.Weight.Catches.F <- list()

# Collect Spawning biomass
Fem_SB <- array(0, dim=c(MaxAge, MaxYear))
FM.Spawning.Bio.All <- list()
FM.Spawning.Bio.NTZ <- list()
FM.Spawning.Bio.F <- list()

# Collect Total Biomass
survived.age <- array(0, dim=c(MaxAge, MaxYear))
survived.age.bio <- array(0, dim=c(MaxAge, MaxYear))
FM.Tot.Bio.All <- list()
FM.Tot.Bio.NTZ <- list()
FM.Tot.Bio.F <- list()

# Yield per Recruit
YPR.by.age <- array(0, dim=(c(MaxAge, MaxYear)))
FM.YPR.All <- list()
FM.YPR.NTZ <- list()
FM.YPR.F <- list()

# Inside outside plots
# Sp.Pop.F <- array(0, dim=c(length(shallow_F_ID), MaxAge, MaxYear))
# Sp.Pop.NTZ <- array(0, dim=c(length(shallow_NTZ_ID), MaxAge, MaxYear))
# 
# SIM.Sp.F <- list()
# SIM.SP.NTZ <- list()

Selectivity <- Selectivity[,,45:59] # Have to modify this to just be years where selectivity is the same as it is now, otherwise we get the selectivity where everything gets fished
Selectivity <- abind(Selectivity, Selectivity, along=3)
Selectivity <- abind(Selectivity, Selectivity, along=3)

##### RUN MSY MODEL #####
# Need to get it to produce the plots we want as well 
# Want it to run for a year and then get the values for the population at the end of the year 
setwd(sg_dir)
Cells <- c(shallow_F_ID, shallow_NTZ_ID)

for(FM in 1:16){
  
  print(FM)
  
  YearlyTotal <- readRDS(paste0(model.name, sep="_","BurnInPop"))
  
  YearlyTotal <- YearlyTotal[as.numeric(Cells),,]
  
  Effort <- Effort.List[[FM]]
  
  
  for (YEAR in 0:(MaxYear-1)){ # Start of model year loop - set it as 45 to make sure the selectivity is the same as present day
    
    ## Loop over all the Rcpp functions in the model
    
    ModelOutput <- RunModelfunc_cpp(YEAR, MaxAge, MaxYear, MaxCell, NatMort, BHa, BHb, PF, 
                                    AdultMove, Mature, Weight, Settlement, 
                                    YearlyTotal, Selectivity, Effort)
    
    survived.age[ ,YEAR+1] <- colSums(ModelOutput$YearlyTotal[,12,]) #No. fish that survived the year in each age class
    survived.age.bio[ ,YEAR+1] <- survived.age[,YEAR+1] * Weight[,12]
    
    ## For BMSY plots inside and outside NTZ
    # Sp.Pop.F[,,YEAR+1] <- ModelOutput$YearlyTotal[c(shallow_F_ID),12, ] # Saving the population at the end of the year in cells <30m depth for plots
    # Sp.Pop.NTZ[,,YEAR+1] <- ModelOutput$YearlyTotal[c(shallow_NTZ_ID),12, ] # Saving the population at the end of the year in cells <30m depth for plots
    
    ## Catch data
    monthly.catch.weight <- ModelOutput$month_catch_weight
    catch.by.weight[ ,YEAR+1] <- rowSums(monthly.catch.weight[,,3:30], dims=1)
    
    age.catch[,,YEAR+1] <- colSums(ModelOutput$month_catch) #This is the number of fish in each age class caught in each month
    catch.by.age[,YEAR+1] <- colSums(age.catch[,,YEAR+1]) # number of fish caught at by the end of the year in each age class
    
    ## Spawning Biomass
    Fem_SB[, YEAR+1] <- ModelOutput$Fem_SB
    
    ## Yield Per Recruit 
    # monthly.YPR.NTZ <- colSums(ModelOutput$Monthly_YPR_age[c(shallow_NTZ_ID),, ])
    # monthly.YPR.F <- colSums(ModelOutput$Monthly_YPR_age[c(shallow_F_ID),, ])

    monthly.YPR <- colSums(ModelOutput$Monthly_YPR_age)
    YPR.by.age[ ,YEAR+1] <- colSums(monthly.YPR[,1:30])
    
    # YPR.by.age.NTZ[ ,YEAR+1] <- colSums(monthly.YPR.NTZ[,1:30])
    # YPR.by.age.F[ ,YEAR+1] <- colSums(monthly.YPR.F[,1:30])
    

  } # End of model year loop
  
  # SIM.Sp.F[[FM]] <- Sp.Pop.F
  # SIM.SP.NTZ[[FM]] <- Sp.Pop.NTZ
  
  FM.Weight.Catches.All[[FM]] <- catch.by.weight
  FM.Age.Catches.All[[FM]] <- catch.by.age
  FM.Spawning.Bio.All[[FM]] <- Fem_SB
  FM.Tot.Bio.All[[FM]] <- survived.age.bio
  FM.YPR.All[[FM]] <- YPR.by.age

}

MSY.Plot.Data.All <- as.data.frame(array(0, dim=c(16,5))) %>%
  rename(Fishing.Mort = "V1",
         YPR = "V2",
         Total.Bio = "V3",
         Spawn.Bio = "V4",
         Zone = "V5") %>%
  mutate(Fishing.Mort = seq(0,0.15, 0.01))

MSY.Plot.Data.NTZ <- as.data.frame(array(0, dim=c(16,5))) %>%
  rename(Fishing.Mort = "V1",
         YPR = "V2",
         Total.Bio = "V3",
         Spawn.Bio = "V4",
         Zone = "V5") %>%
  mutate(Fishing.Mort = seq(0,0.15, 0.01))

MSY.Plot.Data.F <- as.data.frame(array(0, dim=c(16,5))) %>%
  rename(Fishing.Mort = "V1",
         YPR = "V2",
         Total.Bio = "V3",
         Spawn.Bio = "V4",
         Zone = "V5") %>%
  mutate(Fishing.Mort = seq(0,0.15, 0.01))

MSY.Plot.Data <- list()
MSY.Plot.Data[[1]] <- MSY.Plot.Data.All
MSY.Plot.Data[[2]] <- MSY.Plot.Data.NTZ
MSY.Plot.Data[[3]] <- MSY.Plot.Data.F

FM.YPR <- list()
FM.YPR[[1]] <- FM.YPR.All
FM.YPR[[2]] <- FM.YPR.NTZ
FM.YPR[[3]] <- FM.YPR.F


## Yield Per Recruit plot

for(A in 1:1){
  
  area <- FM.YPR[[A]]
  Plot.Data <-  MSY.Plot.Data[[A]]
  
  for(FM in 1:16){
    temp <- area[[FM]]
    temp <- sum(temp[,50])
    
    Plot.Data[FM,"YPR"] <- as.data.frame(temp)
  }
  MSY.Plot.Data[[A]] <- Plot.Data
  
}

Plots.YPR <- list()
for(A in 1:3){
  Plot.Data <- MSY.Plot.Data[[A]]
  
  plot <- ggplot(Plot.Data)+
    geom_line(aes(x=Fishing.Mort, y=YPR))+
    theme_classic()+
    ggplot2::annotate("text", x=0.0, y=0.4, label=A, size = 2.5, fontface=1)
  
  Plots.YPR[[A]] <- plot
}
Plots.YPR[[3]]


## Biomass Plot

FM.Tot.Bio <- list()
FM.Tot.Bio[[1]] <- FM.Tot.Bio.All
FM.Tot.Bio[[2]] <- FM.Tot.Bio.NTZ
FM.Tot.Bio[[3]] <- FM.Tot.Bio.F

for(A in 1:1){
  
  area <- FM.Tot.Bio[[A]]
  Plot.Data <-  MSY.Plot.Data[[A]]
  
  for(FM in 1:16){
    temp <- area[[FM]] * Weight[,12]
    temp <- sum(temp[,50])
    
    Plot.Data[FM,"Total.Bio"] <- temp
  }
 
   MSY.Plot.Data[[A]] <- Plot.Data

}

Plots.Bio <- list()
for(A in 1:3){
  Plot.Data <- MSY.Plot.Data[[A]]
  
  plot <- ggplot(Plot.Data)+
    geom_line(aes(x=Fishing.Mort, y=Total.Bio))+
    theme_classic()+
    ggplot2::annotate("text", x=0.0, y=0.4, label=A, size = 2.5, fontface=1)
  
  Plots.Bio[[A]] <- plot
}
Plots.Bio[[3]]





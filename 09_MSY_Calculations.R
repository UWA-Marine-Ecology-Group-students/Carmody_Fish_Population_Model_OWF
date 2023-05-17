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

#### SET UP FISHING EFFORT FOR EACH LEVEL OF F ####

## Effort values 
q <- 0.000005

F_finite_values <- seq(0, 0.9, 0.05) # These are the values of F that we want to cycle through to see where our MSY is
E_values <- as.data.frame(array(0, dim = c(19, 13))) %>% 
  mutate(V1 = F_finite_values) 
  # mutate(V1 = V1/q)

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

#### CALCULATE CATCHABILITY VARYING BY CELL SIZE ####
water <- water %>% 
  mutate(Area = as.vector((water$cell_area)/1000000))

water_area <- water %>% 
  dplyr::select(Fished_1960, Fished_1987, Fished_2005, Fished_2017, Area) %>% 
  mutate(Area_60 = ifelse(Fished_1960 %in% c("N"), 0, Area)) %>% 
  mutate(Area_87 = ifelse(Fished_1987 %in% c("N"), 0, Area)) %>% 
  mutate(Area_05 = ifelse(Fished_2005 %in% c("N"), 0, Area)) %>% 
  mutate(Area_17 = ifelse(Fished_2017 %in% c("N"), 0, Area)) %>% 
  dplyr::select(Area_60, Area_87, Area_05, Area_17) %>% 
  mutate(Sum_60 = sum(Area_60)) %>% 
  mutate(Sum_87 = sum(Area_87)) %>% 
  mutate(Sum_05 = sum(Area_05)) %>% 
  mutate(Sum_17 = sum(Area_17)) %>% 
  st_drop_geometry() %>% 
  mutate(q_60 = Area_60/Sum_60) %>% 
  mutate(q_87 = Area_87/Sum_87) %>% 
  mutate(q_05 = Area_05/Sum_05) %>% 
  mutate(q_17 = Area_17/Sum_17) %>% 
  mutate(ID = row_number())

spatial_q <- array(0.000005, dim=c(NCELL, 59))

for (y in 31:51){
  spatial_q[ ,y] <- spatial_q[ ,y-1] * 1.02
}

spatial_q[,52:59] <- spatial_q[,51]

for(COL in 1:27){
  for (ROW in 1:NCELL){
    
    spatial_q[ROW,COL] <- 0.3*(spatial_q[ROW,COL]/water_area[ROW,9])
    
  }
}

for(COL in 28:45){
  for (ROW in 1:NCELL){
    
    spatial_q[ROW,COL] <- 0.3*(spatial_q[ROW,COL]/water_area[ROW,10])
    
  }
}

for(COL in 46:57){
  for (ROW in 1:NCELL){
    
    spatial_q[ROW,COL] <- 0.3*(spatial_q[ROW,COL]/water_area[ROW,11])
    
  }
}

for(COL in 58:59){
  for (ROW in 1:NCELL){
    
    spatial_q[ROW,COL] <- 0.3*(spatial_q[ROW,COL]/water_area[ROW,12])
    
  }
}

spatial_q[spatial_q == Inf] <- 0

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
  mutate(Area = as.vector((water$cell_area)/1000000)) #Cells are now in km^2 but with no units

## Work out the utilities for each cell
BR_U <- as.data.frame(matrix(0, nrow=NCELL, ncol=4)) #Set up data frame to hold utilities of cells

BR_U <- BR_U %>% #make sure you give the columns good names so that you know what they are
  rename("Bd_U"=V1) %>% 
  rename("ExM_U" = V2) %>% 
  rename("Tb_U" = V3) %>% 
  rename("CrB_U"=V4)

# Add coefficients to each variable - all the BRs are the same but make them negative because the further away they are the less likely people are to visit
a = 0.0005
b = -0.001
c = -0.001
d = -0.001
e = -0.001

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

# Fishing_MSY <- Fishing_MSY*q



Effort.List <- list()

for(FM in 1:19){
  
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
NCELL = nrow(water)
MaxAge = 30
MaxYear = 50
MaxCell = NCELL
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
FM.Age.Catches <- list()

catch.by.weight <- array(0, dim=c(MaxCell, MaxYear))
FM.Weight.Catches <- list()

# Collect Spawning biomass
Fem_SB <- array(0, dim=c(MaxAge, MaxYear))
FM.Spawning.Bio <- list()

# Collect Total Biomass
survived.age <- array(0, dim=c(MaxAge, MaxYear))
survived.age.bio <- array(0, dim=c(MaxAge, MaxYear))
FM.Tot.Bio <- list()

# Yield per Recruit
YPR.by.age <- array(0, dim=(c(MaxAge, MaxYear)))
YPR.by.age.NTZ <- array(0, dim=(c(MaxAge, MaxYear)))
YPR.by.age.F <- array(0, dim=(c(MaxAge, MaxYear)))
FM.YPR <- list()
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

for(FM in 1:10){
  
  print(FM)
  
  YearlyTotal <- readRDS(paste0(model.name, sep="_","BurnInPop"))
  
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
  
  FM.Weight.Catches[[FM]] <- catch.by.weight
  FM.Age.Catches[[FM]] <- catch.by.age
  FM.Spawning.Bio[[FM]] <- Fem_SB
  FM.Tot.Bio[[FM]] <- survived.age.bio
  FM.YPR[[FM]] <- YPR.by.age
  # FM.YPR.NTZ[[FM]] <- YPR.by.age.NTZ
  # FM.YPR.F[[FM]] <- YPR.by.age.F

}

MSY.Plot.Data <- as.data.frame(array(0, dim=c(10,5))) %>% 
  rename(Fishing.Mort = "V1",
         YPR = "V2",
         Total.Bio = "V3",
         Spawn.Bio = "V4",
         Zone = "V5") %>% 
  mutate(Fishing.Mort = seq(0,0.45, 0.05))


## Yield Per Recruit plot
for(FM in 1:10){
  temp <- FM.YPR[[FM]]
  temp <- sum(temp[,50])

  MSY.Plot.Data[FM,"YPR"] <- as.data.frame(temp)
}

# for(FM in 1:10){
#   temp <- FM.YPR.NTZ[[FM]]
#   temp <- sum(temp[,50]) %>% 
#     as.data.frame() %>% 
#     mutate(Zone="NTZ")
#   
#   MSY.Plot.Data[FM,c("YPR", "Zone")] <- as.data.frame(temp) 
#   
#   temp <- FM.YPR.F[[FM]]
#   temp <- sum(temp[,50]) %>% 
#     as.data.frame() %>% 
#     mutate(Zone="F")
#   
#   MSY.Plot.Data[FM+10,c("YPR", "Zone")] <- as.data.frame(temp) 
# }
ggplot(MSY.Plot.Data)+
  geom_line(aes(x=Fishing.Mort, y=YPR))+
  theme_classic()


## Biomass Plot
for(FM in 1:10){
  temp <- FM.Tot.Bio[[FM]]
  temp <- sum(temp[,50])

  MSY.Plot.Data[FM,"Total.Bio"] <- temp
}

# for(FM in 1:10){
#   temp <- colSums(SIM.SP.NTZ[[FM]]) * Weight[,12]
#   temp <- sum(temp[,50])
#   
#   MSY.Plot.Data[FM,"Total.Bio"] <- temp
#   
#   temp <- colSums(SIM.Sp.F[[FM]]) * Weight[,12]
#   temp <- sum(temp[,50])
#   
#   MSY.Plot.Data[FM+10,"Total.Bio"] <- temp
# }

ggplot(MSY.Plot.Data)+
  geom_line(aes(x=Fishing.Mort, y=Total.Bio))+
  theme_classic()

## Spawning Biomass Plot
for(FM in 1:10){
  temp <- FM.Spawning.Bio[[FM]]
  temp <- sum(temp[,50])
  
  MSY.Plot.Data[FM,"Spawn.Bio"] <- temp
}

ggplot(MSY.Plot.Data)+
  geom_line(aes(x=Fishing.Mort, y=Spawn.Bio))+
  theme_classic()


#### CALCULATE F FROM MY MODEL DATA ####
setwd(pop_dir)

Weight <- as.data.frame(Weight)
Mature <- as.data.frame(Mature)

Age.Dist <- list()
Age.Catch <- list()

# Each layer is a simulation, rows are cells and columns are years
Age.Dist[[1]] <- readRDS(paste0(model.name, sep="_", "Age_Distribution", sep="_", "S00"))
Age.Dist[[2]] <- readRDS(paste0(model.name, sep="_", "Age_Distribution", sep="_", "S01"))
Age.Dist[[3]] <- readRDS(paste0(model.name, sep="_", "Age_Distribution", sep="_", "S02"))
Age.Dist[[4]] <- readRDS(paste0(model.name, sep="_", "Age_Distribution", sep="_", "S03"))

Age.Catch[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Age", sep="_", "S00"))
Age.Catch[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Age", sep="_", "S01"))
Age.Catch[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Age", sep="_", "S02"))
Age.Catch[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Age", sep="_", "S03"))

Names <- c("Historical and Current Management", "No Spatial Management", 
           "Temporal and Spatial Management","Temporal Management Only" )

## Calculate mean number of fish and biomass in each year of the model
Ages <- list()
Catches <- list()
catch <- array(0, dim=c(30,59,100))
Biomass <- list()
Spawning.Biomass <- list()


for(S in 1:4){
  
  # Ages
  temp <- Age.Dist[[S]]
  
  temp2 <- temp %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) 
  
  Ages[[S]] <- temp2
  
  # Catches
  temp.c <- Age.Catch[[S]]
  for(SIM in 1:100){
    
    catch[,,SIM] <- temp.c[[SIM]]
  
  }
  
  temp2.c <- catch %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) 
  
  Catches[[S]] <- temp2.c
  
  #Biomass
  temp3 <- temp %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) %>% 
    mutate_all(.,function(col){Weight$V12*col}) 
    
  Biomass[[S]] <- temp3
  
  #Spawning biomass
  temp4 <- temp %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) %>% 
    mutate_all(.,function(col){Mature$V12*col}) %>% 
    mutate_all(.,function(col){Weight$V12*col}) 
  
  Spawning.Biomass[[S]] <- temp4*0.5
  
}


## Calculate fishing mortality
Mortality <- array(0, dim=c(59, 3))
FM.Scenario <- list()

for(S in 1:4){
  
  temp <- Ages[[S]]
  temp <- temp %>% 
    colSums(.) %>% 
    as.data.frame(.) %>% 
    rename(Pop = ".")
  
  temp.c <- Catches[[S]]
  temp.c <- temp.c %>% 
    colSums(.) %>% 
    as.data.frame(.) 

  temp2 <- temp %>% 
    mutate(Caught = temp.c) %>% 
    mutate(FM = Caught/Pop)
      
  FM.Scenario[[S]] <- temp2$FM
  
}

## Biomass and Spawning Biomass
Kobe.Data <- NULL

for(S in 1:4){
  
  temp <- FM.Scenario[[S]]
  temp2 <- Spawning.Biomass[[S]]
  temp3 <- Biomass[[S]]
  
  temp2 <- temp2 %>% 
    colSums()
  temp3 <- temp3 %>% 
    colSums()
  
  temp <- temp %>% 
    as.data.frame() %>% 
    rename(FM = ".") %>% 
    mutate(Biomass = temp3) %>% 
    mutate(Spawning.Biomass = temp2) %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Year = seq(1960,2018, 1))

  
  Kobe.Data <- rbind(Kobe.Data, temp)
  
}

#### Make Kobe Plots ####
F.MSY = 0.1
Bio.MSY = 30449.206
SB.MSY = 12818.6052

Kobe.Data <- Kobe.Data %>% 
  mutate(Rel.FM = as.numeric(FM/F.MSY)) %>% 
  mutate(Rel.Bio = as.numeric(Biomass/Bio.MSY))%>% 
  mutate(Rel.SB = as.numeric(Spawning.Biomass/SB.MSY)) %>% 
  rename(Mod.Year = "Year") %>% 
  mutate(Mod.Year = as.character(Mod.Year))

Bio.Plot.S00 <- Kobe.Data %>% 
  filter(Scenario %in% c("Historical and Current Management")) %>% 
  ggplot(., aes(x=Rel.Bio, y=Rel.FM, label=Mod.Year))+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 3, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 3.5,fill = "red")+
  annotate("rect", xmin = 1, xmax = 3, ymin = 1, ymax = 3.5,fill = "orange")+
  geom_point()+
  geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_path()+
  scale_y_continuous(breaks=seq(0,3.5, 0.5), limits=c(0,3.5)) +
  scale_x_continuous(breaks=seq(0,3, 0.5), limits=c(0,3)) +
  theme_classic()+
  ylab("F/Fmsy")+
  xlab("B/Bmsy")


Bio.Plot.S01 <- Kobe.Data %>% 
  filter(Scenario %in% c("No Spatial Management")) %>% 
  ggplot(., aes(x=Rel.Bio, y=Rel.FM, label=Mod.Year))+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 3, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 3.5,fill = "red")+
  annotate("rect", xmin = 1, xmax = 3, ymin = 1, ymax = 3.5,fill = "orange")+
  geom_point()+
  geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_path()+
  scale_y_continuous(breaks=seq(0,3.5, 0.5), limits=c(0,3.5)) +
  scale_x_continuous(breaks=seq(0,3, 0.5), limits=c(0,3)) +
  theme_classic()+
  ylab("F/Fmsy")+
  xlab("B/Bmsy")

Bio.Plot.S02 <- Kobe.Data %>% 
  filter(Scenario %in% c("Temporal Management Only")) %>% 
  ggplot(., aes(x=Rel.Bio, y=Rel.FM, label=Mod.Year))+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 3, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 3.5,fill = "red")+
  annotate("rect", xmin = 1, xmax = 3, ymin = 1, ymax = 3.5,fill = "orange")+
  geom_point()+
  geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_path()+
  scale_y_continuous(breaks=seq(0,3.5, 0.5), limits=c(0,3.5)) +
  scale_x_continuous(breaks=seq(0,3, 0.5), limits=c(0,3)) +
  theme_classic()+
  ylab("F/Fmsy")+
  xlab("B/Bmsy")

Bio.Plot.S03 <- Kobe.Data %>% 
  filter(Scenario %in% c("Temporal and Spatial Management")) %>% 
  ggplot(., aes(x=Rel.Bio, y=Rel.FM, label=Mod.Year))+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 3, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 3.5,fill = "red")+
  annotate("rect", xmin = 1, xmax = 3, ymin = 1, ymax = 3.5,fill = "orange")+
  geom_point()+
  geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_path()+
  scale_y_continuous(breaks=seq(0,3.5, 0.5), limits=c(0,3.5)) +
  scale_x_continuous(breaks=seq(0,3, 0.5), limits=c(0,3)) +
  theme_classic()+
  ylab("F/Fmsy")+
  xlab("B/Bmsy")


#### CHECKING KOBE PLOT CORRESPONDS TO CATCH ####
yearly.catch <- NULL

for(S in 1:4){
  
  temp <- Catches[[S]]
  
  temp2 <- temp %>% 
    colSums(.) %>% 
    as.data.frame() %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Mod.Year = seq(1960, 2018, 1)) %>% 
    rename(Catch = ".")
  
  yearly.catch <- rbind(yearly.catch, temp2)
  
}

catch.plot <- yearly.catch %>% 
  filter(Scenario %in% c("Historical and Current Management")) %>% 
  ggplot()+
  geom_line(aes(x=Mod.Year, y=Catch))+
  theme_classic()+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")
  

#### BIOMASS OVER TIME INSIDE AND OUTSIDE NTZ ####
setwd(pop_dir)

NTZ.pop <- list()

NTZ.pop[[1]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S00"))
NTZ.pop[[2]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S01"))
NTZ.pop[[3]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S02"))
NTZ.pop[[4]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S03"))

F.pop <- list()

F.pop[[1]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S00"))
F.pop[[2]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S01"))
F.pop[[3]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S02"))
F.pop[[4]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S03"))

bio.ntz <- array(0, dim=c(30, 59, 100))
bio.F <- array(0, dim=c(30, 59, 100))

bio.mean.NTZ <- NULL
bio.sd.NTZ <- NULL
bio.mean.F <- NULL
bio.sd.F <- NULL

for(S in 1:4){
  
  temp.NTZ <- NTZ.pop[[S]]
  temp.F <- F.pop[[S]]
  
  for(SIM in 1:100){
    
    temp2.NTZ <- temp.NTZ[[SIM]]
    temp3.NTZ <- colSums(temp2.NTZ, dim=1)
    temp4.NTZ <- temp3.NTZ*Weight[,12]
    
    bio.ntz[,,SIM] <- temp4.NTZ
    
    temp2.F <- temp.F[[SIM]]
    temp3.F <- colSums(temp2.F, dim=1)
    temp4.F <- temp3.F*Weight[,12]
    
    bio.F[,,SIM] <- temp4.F
  }
  
  temp.mean.NTZ <- bio.ntz %>% 
    colSums(., dim=1) %>% 
    rowMeans(.) %>% 
    as.data.frame(.) %>% 
    mutate(Zone = "NTZ") %>% 
    rename(Biomass = ".") %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Mod_Year = seq(1960,2018,1))
  
  bio.mean.NTZ <- rbind(bio.mean.NTZ, temp.mean.NTZ)
  
  temp.mean.F <- bio.F %>% 
    colSums(., dim=1) %>% 
    rowMeans(.) %>% 
    as.data.frame(.) %>% 
    mutate(Zone = "F") %>% 
    rename(Biomass = ".") %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Mod_Year = seq(1960,2018,1))
  
  bio.mean.F <- rbind(bio.mean.F, temp.mean.F)
  
  temp.SD.NTZ <- bio.ntz %>% 
    colSums(., dim=1) %>% 
    rowSds(.) %>% 
    as.data.frame(.) %>% 
    mutate(Zone = "NTZ") %>% 
    rename(Biomass = ".") %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Mod_Year = seq(1960,2018,1))
  
  bio.sd.NTZ <- rbind(bio.sd.NTZ, temp.SD.NTZ)
  
  temp.SD.F <- bio.F %>% 
    colSums(., dim=1) %>% 
    rowSds(.) %>% 
    as.data.frame(.) %>% 
    mutate(Zone = "F") %>% 
    rename(Biomass = ".") %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Mod_Year = seq(1960,2018,1))
  
  bio.sd.F <- rbind(bio.sd.F, temp.SD.F)
  
}

bio.mean <- rbind(bio.mean.NTZ, bio.mean.F)
bio.SD <- rbind(bio.sd.NTZ, bio.sd.F)

bio.mean <- bio.mean %>% 
  mutate(SD = bio.SD$Biomass)

#* Plot by scenario and inside outside ####

Pre_1987_NTZ <- bio.mean %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) 

Pre_1987_F <- bio.mean %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F"))

Biomass_F <- bio.mean %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current Management") & Mod_Year>1985, "Historical and\ncurrent management", 
                                                                 ifelse(Scenario %in% c("Temporal Management Only") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("No Spatial Management"), "No NTZs or\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.biomass <- bio.mean %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current Management") & Mod_Year>1985, "Historical and\ncurrent management", 
                                                                 ifelse(Scenario %in% c("Temporal Management Only") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("No Spatial Management"), "No NTZs or\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup))%>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Biomass, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Biomass, ymin=Biomass-SD, ymax=Biomass+SD, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Biomass_F, aes(x=Mod_Year, y=Biomass, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Biomass_F, aes(x=Mod_Year, y=Biomass, ymin=Biomass-SD, ymax=Biomass+SD, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent management"="#36753B", "No NTZs or\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("Pre-1987"="grey20", "Historical and\ncurrent management"="#36753B", "No NTZs or\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                                 "Temporal\nmanagement only"="#BBCC33"), name= "Spatial and Temporal\nManagement Scenario")+ 
  geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Biomass, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Biomass, ymin=Biomass-SD, ymax=Biomass+SD, group=Scenario), fill="grey20",alpha=0.2)+
  geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Biomass, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Biomass, ymin=Biomass-SD, ymax=Biomass+SD, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  ylab("Biomass (Kg)")+
  xlab("Year")+
  xlim(1960,2020)+
  scale_linetype_manual(values = c("longdash", "solid" ), labels=c("Always Fished", "NTZ Area"), name="Model Area")+
  theme(legend.title = element_text(size=9, face="bold"), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")
line.biomass

## Biomass Relative to B MSY plot

F.MSY = 0.1
Bio.MSY.NTZ = 11250.2377
Bio.MSY.F = 2886.8597


bio.mean.NTZ <- bio.mean %>% 
  filter(Zone %in% c("NTZ")) %>% 
  mutate(Rel.Bio = as.numeric(Biomass/Bio.MSY.NTZ))

bio.mean.F <- bio.mean %>% 
  filter(Zone %in% c("F")) %>% 
  mutate(Rel.Bio = as.numeric(Biomass/Bio.MSY.F))

NTZ.bio <- bio.mean.NTZ %>% 
  filter(Scenario %in% c("Historical and Current Management")) %>% 
  ggplot()+
  annotate("rect", xmin = 1960, xmax = 2018, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1960, xmax = 2018, ymin = 1, ymax = 3.5,fill = "green")+
  geom_line(aes(x=Mod_Year, y=Rel.Bio))+
  theme_classic()+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  geom_hline(yintercept=1.0, colour="grey20") +
  ylab("B/Bmsy")+
  xlab("Year")

F.bio <- bio.mean.F %>% 
  filter(Scenario %in% c("Historical and Current Management")) %>% 
  ggplot()+
  annotate("rect", xmin = 1960, xmax = 2018, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1960, xmax = 2018, ymin = 1, ymax = 3.5,fill = "green")+
  geom_line(aes(x=Mod_Year, y=Rel.Bio))+
  theme_classic()+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  geom_hline(yintercept=1.0, colour="grey20") +
  ylab("B/Bmsy")+
  xlab("Year")
  




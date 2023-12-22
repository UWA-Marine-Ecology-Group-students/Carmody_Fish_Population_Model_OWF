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
msy_dir <- paste(working.dir, "MSY_Outputs", sep="/")

model.name <- "ningaloo"

## Read in functions
setwd(working.dir)
sourceCpp("X_Model_RccpArm.cpp")
source("X_Functions.R")

#### READ IN DATA ####
setwd(sg_dir)
AdultMove <- readRDS(paste0(model.name, sep="_", "movement_medium"))
Settlement <- readRDS(paste0(model.name, sep="_","recruitment")) 
Effort <- readRDS(paste0(model.name, sep="_", "fishing"))
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
Water <- readRDS(paste0(model.name, sep="_","water")) %>% 
  st_make_valid()
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
  st_make_valid() %>% 
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
  #filter(bathy >= c(-30)) %>% 
  filter(!is.na(bathy))

shallow_cells_NTZ <- Water_shallow %>% 
  filter(Fished_2017 %in% c("N")) %>% 
  st_drop_geometry(.)

shallow_cells_F <- Water_shallow %>% 
  filter(Fished_2017 %in% c("Y")) %>% 
  st_drop_geometry(.)

shallow_NTZ_ID <- as.numeric(levels(shallow_cells_NTZ$ID))[as.integer(shallow_cells_NTZ$ID)]
shallow_F_ID <- as.numeric(levels(shallow_cells_F$ID))[as.integer(shallow_cells_F$ID)]

# NCELL <- length(shallow_NTZ_ID)

#### SET UP FISHING EFFORT FOR EACH LEVEL OF F ####

## Effort values 

# F_finite_values <- seq(0, 0.9, 0.05) # These are the values of F that we want to cycle through to see where our MSY is
F_finite_values <- 0.1
E_values <- as.data.frame(array(0, dim = c(9, 13))) %>% 
  mutate(V1 = F_finite_values) 
  # mutate(V1 = V1/q)

names(E_values)[1:13] <- c("Yearly_Total", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

## Allocation to months
Monthly_effort <- E_values[,2:13]

for(E in 1:9){
  for(M in 1:12){
    Monthly_effort[E,M] <- month.ave[M,2] * E_values[E,1] 
  }
}

Monthly_effort <- Monthly_effort %>% 
  mutate(FM = F_finite_values)

# Allocate effort to cells
NCELL <- nrow(Water)

Fishing_MSY <- array(0, dim=c(NCELL, 12,9)) #This array has a row for every cell, a column for every month, and a layer for every value of F

Effort.List <- list()

for(FM in 1:9){
  
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
MaxCell = NCELL
PF = 0.5 # Proportion female

## Natural Mortality
# We have instantaneous mortality from Marriott et al 2011 and we need to convert that into monthly mortality
NatMort = 0.146
step = 1/12 # We're doing a monthly time step here

# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
BHa = 0.4344209 #0.4344209
BHb = 0.009398152 #0.003759261 #0.0009398152 #0.01889882

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
Sp.Pop.F <- array(0, dim=c(length(shallow_F_ID), MaxAge, MaxYear))
Sp.Pop.NTZ <- array(0, dim=c(length(shallow_NTZ_ID), MaxAge, MaxYear))
# 
SIM.Sp.F <- list()
SIM.SP.NTZ <- list()

Selectivity <- Selectivity[,,45:59] # Have to modify this to just be years where selectivity is the same as it is now, otherwise we get the selectivity where everything gets fished
Selectivity <- abind(Selectivity, Selectivity, along=3)
Selectivity <- abind(Selectivity, Selectivity, along=3)

##### RUN MSY MODEL #####
# Need to get it to produce the plots we want as well 
# Want it to run for a year and then get the values for the population at the end of the year 
setwd(sg_dir)

for(FM in 1:1){
  
  print(FM)
  
  YearlyTotal <- readRDS(paste0(model.name, sep="_", "BurnInPop"))
  
  #YearlyTotal <- YearlyTotal[as.numeric(Cells),,]

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
    Sp.Pop.NTZ[,,YEAR+1] <- ModelOutput$YearlyTotal[c(shallow_NTZ_ID),12, ] # Saving the population at the end of the year in cells <30m depth for plots
    
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
  
  SIM.Sp.F[[FM]] <- Sp.Pop.F
  SIM.SP.NTZ[[FM]] <- Sp.Pop.NTZ
  
  FM.Weight.Catches.All[[FM]] <- catch.by.weight
  FM.Age.Catches.All[[FM]] <- catch.by.age
  FM.Spawning.Bio.All[[FM]] <- Fem_SB
  FM.Tot.Bio.All[[FM]] <- survived.age.bio
  FM.YPR.All[[FM]] <- YPR.by.age

}

MSY.Plot.Data.All <- as.data.frame(array(0, dim=c(9,5))) %>%
  rename(Fishing.Mort = "V1",
         YPR = "V2",
         Total.Bio = "V3",
         Spawn.Bio = "V4",
         Zone = "V5") %>%
  mutate(Fishing.Mort = seq(0, 0.2, 0.025))


MSY.Plot.Data <- list()
MSY.Plot.Data[[1]] <- MSY.Plot.Data.All
MSY.Plot.Data[[2]] <- MSY.Plot.Data.All
MSY.Plot.Data[[3]] <- MSY.Plot.Data.F

FM.YPR <- list()
FM.YPR[[1]] <- FM.YPR.All
FM.YPR[[2]] <- FM.YPR.NTZ
FM.YPR[[3]] <- FM.YPR.F


## Yield Per Recruit plot

for(A in 1:1){
  
  area <- FM.YPR[[A]]
  Plot.Data <-  MSY.Plot.Data[[A]]
  
  for(FM in 1:9){
    temp <- area[[FM]]
    temp <- sum(temp[,50])
    
    Plot.Data[FM,"YPR"] <- as.data.frame(temp)
  }
  MSY.Plot.Data[[A]] <- Plot.Data
  
}

Plots.YPR <- list()
for(A in 1:1){
  Plot.Data <- MSY.Plot.Data[[A]]
  
  plot <- ggplot(Plot.Data)+
    geom_line(aes(x=Fishing.Mort, y=YPR))+
    theme_classic()
  
  Plots.YPR[[A]] <- plot
}
Plots.YPR[[1]]


## Biomass Plot

FM.Tot.Bio <- list()
FM.Tot.Bio[[1]] <- FM.Tot.Bio.All
FM.Tot.Bio[[2]] <- FM.Tot.Bio.NTZ
FM.Tot.Bio[[3]] <- FM.Tot.Bio.F

for(A in 1:2){
  
  area <- FM.Tot.Bio[[A]]
  Plot.Data <-  MSY.Plot.Data[[A]]
  
  for(FM in 1:21){
    temp <- area[[FM]] * Weight[,12]
    temp <- sum(temp[,50])
    
    Plot.Data[FM,"Total.Bio"] <- temp
  }
 
   MSY.Plot.Data[[A]] <- Plot.Data

}

Plots.Bio <- list()
for(A in 1:1){
  Plot.Data <- MSY.Plot.Data[[A]]
  
  plot <- ggplot(Plot.Data)+
    geom_line(aes(x=Fishing.Mort, y=Total.Bio))+
    theme_classic()+
    ggplot2::annotate("text", x=0.0, y=0.4, label=A, size = 2.5, fontface=1)
  
  Plots.Bio[[A]] <- plot
}
Plots.Bio[[1]]


## Spawning biomass

for(A in 1:1){
  
  area <- FM.Tot.Bio[[A]]
  Plot.Data <-  MSY.Plot.Data[[A]]
  
  for(FM in 1:9){
    temp <- area[[FM]] * Mature[,12]
    temp <- temp * Weight[,12]
    temp <- sum(temp[,50])
    
    Plot.Data[FM,"Spawn.Bio"] <- temp
  }
  
  MSY.Plot.Data[[A]] <- Plot.Data
  
}
MSY.Plot.Data[[1]]

## Spawning biomass inside NTZ

for(A in 2:2){
  

  Plot.Data <-  MSY.Plot.Data[[A]]
  
  for(FM in 1:9){
    area <- SIM.SP.NTZ[[FM]]
    temp <- area[,,50] * Mature[,12]
    temp <- temp * Weight[,12]
    temp <- sum(temp)
    
    Plot.Data[FM,"Spawn.Bio"] <- temp
  }
  
  MSY.Plot.Data[[A]] <- Plot.Data
  
}
MSY.Plot.Data[[2]]

#### SETTING F TO DIFFERENT LEVELS FOR KOBE PLOT ####
## Set Effort to be whatever you want for the Kobe plot
Fishing_mort <- "1-5xNatMort"
#F_finite_values <- seq(0, 0.9, 0.05) # These are the values of F that we want to cycle through to see where our MSY is
F_finite_values <- 1.5*NatMort
E_values <- as.data.frame(array(0, dim = c(9, 13))) %>% 
  mutate(V1 = F_finite_values) 
# mutate(V1 = V1/q)

names(E_values)[1:13] <- c("Yearly_Total", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

## Allocation to months
Monthly_effort <- E_values[,2:13]

for(E in 1:9){
  for(M in 1:12){
    Monthly_effort[E,M] <- month.ave[M,2] * E_values[E,1] 
  }
}

Monthly_effort <- Monthly_effort %>% 
  mutate(FM = F_finite_values)

# Allocate effort to cells
NCELL <- nrow(Water)

Fishing_MSY <- array(0, dim=c(NCELL, 12,9)) #This array has a row for every cell, a column for every month, and a layer for every value of F

Effort.List <- list()

for(FM in 1:9){
  
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

## Set selectivity to be what it is now
Selectivity <- Selectivity[,,45:59] # Have to modify this to just be years where selectivity is the same as it is now, otherwise we get the selectivity where everything gets fished
Selectivity <- abind(Selectivity, Selectivity, along=3)
Selectivity <- abind(Selectivity, Selectivity, along=3)

#### PARAMETER VALUES ####
## Natural Mortality
# We have instantaneous mortality from Marriott et al 2011 and we need to convert that into monthly mortality
NatMort = 0.146

## Potentially 18% mortality rate per day for larvae https://www.nature.com/articles/s41598-022-25629-w
## Larval duration is thought to be ~45 days https://www.researchgate.net/publication/221894788_Understanding_age-specific_dispersal_in_fishes_through_hydrodynamic_modelling_genetic_simulations_and_microsatellite_DNA_analysis#pf9


# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
BHa = 0.4344209 #0.4344209
BHb = 0.009398152 #0.003759261 #0.0009398152 #0.01889882
PF = 0.5

# Model settings
MaxCell = nrow(Water) # Number of cells in the model
MaxAge = 30 # The max age of the fish in the model -1 to account for the fact that Rcpp functions start from 0
MaxYear = 59 # Number of years the model should run for -1 to account for the fact that Rcpp functions start from 0
# PlotTotal <- T #This is whether you want a line plot of the total or the map

Pop.Groups <- seq(1,12)

#### SET UP LISTS TO HOLD THE PLOTS ####
SpatialPlots <- list()
LengthPlots <- list()
AgePlots <- list()
TimesPlotted = 0

#### SET UP INITIAL POPULATION ####

## To save population information
PopTotal <- array(0, dim=c(MaxCell, 12, MaxYear)) # This is our total population, all ages are summed and each column is a month (each layer is a year)
Total <- array(NA, dim=c(MaxYear,1)) # For plotting
Pop.Total.Dist <- array(0, dim=c(MaxCell, MaxYear))

## To save catch information
age.catch <- array(0, dim=c(12, MaxAge, MaxYear))
catch.by.cell <- array(0, dim=c(MaxCell, MaxYear))
catch.by.age <- array(0, dim=(c(MaxAge, MaxYear)))
catch.by.weight <- array(0, dim=(c(MaxCell, MaxYear)))

## Save all information by simulation
Sim.Pop <- array(0, dim=c(MaxYear, 100))
Sim.Ages <- array(0, dim=c(MaxAge, MaxYear, 100))

Sp.Pop.F <- array(0, dim=c(length(shallow_F_ID), MaxAge, MaxYear))
Sp.Pop.NTZ <- array(0, dim=c(length(shallow_NTZ_ID), MaxAge, MaxYear))

SIM.Sp.F <- list()
SIM.Sp.NTZ <- list()
SIM.N.Dist <- list()
SIM.N.Catches <- list()
SIM.Age.Catches <- list()
SIM.Weight.Catches <- list()

#### RUN MODEL ####
Start=Sys.time()
for (SIM in 1:50){ # Simulation loop
  
  #### SET UP LISTS TO HOLD THE PLOTS ####
  SpatialPlots <- list()
  LengthPlots <- list()
  AgePlots <- list()
  TimesPlotted = 0
  
  #### SET UP INITIAL POPULATION ####
  PopTotal <- array(0, dim=c(MaxCell, 12, MaxYear)) # This is our total population, all ages are summed and each column is a month (each layer is a year)
  Total <- array(NA, dim=c(MaxYear,1)) # For plotting
  
  print(paste0("SIM", sep=" ", SIM))
  
  setwd(sg_dir)
  YearlyTotal <- readRDS(paste0(model.name, sep="_", "BurnInPop_High_M"))
  
  for (YEAR in 25:(MaxYear-1)){ # Start of model year loop
    
    print(YEAR)
    
    ## Loop over all the Rcpp functions in the model
    ModelOutput <- RunModelfunc_cpp(YEAR, MaxAge, MaxYear, MaxCell, NatMort, BHa, BHb, PF, AdultMove, Mature, Weight, Settlement, 
                                    YearlyTotal, Selectivity, Effort)
    
    ## Get outputs from the model
    # Have to add 1 to all YEAR because the loop is now starting at 0
    ## Abundance in different areas
    
    PopTotal[ , ,YEAR+1] <- rowSums(ModelOutput$YearlyTotal[,,1:30], dim=2) # This flattens the matrix to give you the number of fish present in the population each month in each cell, with layers representing the year
    
    # Whole area
    Water$pop <- PopTotal[ , 12, YEAR+1] # We just want the population at the end of the year
    Total[YEAR+1,1] <- sum(Water$pop) # Add this to a dataframe we can then use later 
    
    # By zone
    Sp.Pop.F[,,YEAR+1] <- ModelOutput$YearlyTotal[c(shallow_F_ID),12, ] # Saving the population at the end of the year in cells <30m depth for plots
    Sp.Pop.NTZ[,,YEAR+1] <- ModelOutput$YearlyTotal[c(shallow_NTZ_ID),12, ] # Saving the population at the end of the year in cells <30m depth for plots
    
    # By cell so we can get distances to boat ramps
    Pop.Total.Dist[ ,YEAR+1] <- PopTotal[,12,YEAR+1]
    
    ## Catch data
    monthly.catch <- ModelOutput$month_catch
    age.catch[,,YEAR+1] <- colSums(ModelOutput$month_catch) #This is the number of fish in each age class caught in each month
    
    catch.by.cell[,YEAR+1] <- rowSums(monthly.catch[,,3:30], dims=1) # Number of legal size fish caught in each cell 
    catch.by.age[,YEAR+1] <- colSums(age.catch[,,YEAR+1]) # number of fish caught at by the end of the year in each age class
    
    monthly.catch.weight <- ModelOutput$month_catch_weight
    catch.by.weight[ ,YEAR+1] <- rowSums(monthly.catch.weight[,,3:30], dims=1)
    
    # if(YEAR==58){ # Have to subtract one because the loop now starts at 0
    #   TotalPop <- as.data.frame(Total) %>%
    #     rename(Tot.Pop="V1")
    #   TotalPop$Year <- seq(1960, 2018, by=1)
    #   TotalPlot <- total.plot.func(pop=TotalPop)
    #   print(TotalPlot)
    # } else { }
    # 
    # if(YEAR %%5==0|YEAR==58){
    #   TimesPlotted <- TimesPlotted+1
    #   #SpatialPlots[[TimesPlotted]] <- spatial.plot.func(area=Water, pop=Total, pop.breaks=pop.groups, colours="PuBu")
    #   #AgePlots[[TimesPlotted]] <- age.plot.func(pop=YearlyTotal, NTZs=NoTake)
    #   
    # } else { }
    
    Sim.Pop[YEAR+1, SIM] <- Total[YEAR+1,1] 
    Sim.Ages[ ,YEAR+1,SIM] <- colSums(ModelOutput$YearlyTotal[,12,1:30]) #Number of fish present in age age group at the end of the year
  } # End of model year loop
  
  ## Population in different zones
  SIM.Sp.F[[SIM]] <- Sp.Pop.F
  SIM.Sp.NTZ[[SIM]] <- Sp.Pop.NTZ
  SIM.N.Dist[[SIM]] <- Pop.Total.Dist
  
  ## Catches
  SIM.N.Catches[[SIM]] <- catch.by.cell # Catches in each cell
  SIM.Age.Catches[[SIM]] <- catch.by.age # Catches by age in each month of the year in each year
  SIM.Weight.Catches[[SIM]] <- catch.by.weight
  
  if(SIM==100){ # Saving if statement
    setwd(msy_dir)
    
    ## Population
    # Total population
    filename <- paste0(model.name, sep="_", "Total_Population", sep="_", Scenario, sep="_", Fishing_mort)
    saveRDS(Sim.Pop, file=filename)
    
    # Numbers of each age that make it to the end of each year
    filename <- paste0(model.name, sep="_", "Age_Distribution", sep="_", Scenario, sep="_", Fishing_mort)
    saveRDS(Sim.Ages, file=filename)
    
    # Numbers of fish of each age, inside and outside sanctuary zones
    filename <- paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", Scenario, sep="_", Fishing_mort)
    saveRDS(SIM.Sp.NTZ, file=filename)
    
    filename <- paste0(model.name, sep="_", "Sp_Population_F", sep="_", Scenario, sep="_", Fishing_mort)
    saveRDS(SIM.Sp.F, file=filename)
    
    # Number of fish in each cell at the end of each year
    filename <- paste0(model.name, sep="_", "Cell_Population", sep="_", Scenario, sep="_", Fishing_mort)
    saveRDS(SIM.N.Dist, file=filename)
    
    ## Catches
    # Catch in each year by age
    filename <- paste0(model.name, sep="_", "Catch_by_Age", sep="_", Scenario, sep="_", Fishing_mort)
    saveRDS(SIM.Age.Catches, file=filename)
    
    # catch in each cell across the year
    filename <- paste0(model.name, sep="_", "Catch_by_Cell", sep="_", Scenario, sep="_", Fishing_mort)
    saveRDS(SIM.N.Catches, file=filename)
    
    # Catch in each cell by weight across the year
    filename <- paste0(model.name, sep="_", "Catch_by_Weight", sep="_", Scenario, sep="_", Fishing_mort)
    saveRDS(SIM.Weight.Catches, file=filename)
    
  } else { }# End saving if statement
  
} # End simulation loop
End = Sys.time() 
Runtime = End - Start
Runtime

finished_email <- gm_mime() %>% 
  gm_to("charlotte.aston@research.uwa.edu.au") %>% 
  gm_from("charlotte.aston@marineecology.io") %>% 
  gm_subject("MSY code has finished running") %>% 
  gm_text_body("F at M scenario done, Yew!")

d <- gm_create_draft(finished_email)
gm_send_draft(d)


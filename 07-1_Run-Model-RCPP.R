###################################################

# Running the full model with whole area and fishing
# Produces plots every five years
# Get spatial plots as well as age based histograms
# See script six for more complex plots
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

rm(list = ls())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")
pop_dir <-  paste(working.dir, "Output_Population", sep="/")
sim_dir <- paste(working.dir, "Simulations", sep="/")

## Read in functions
setwd(working.dir)
sourceCpp("Model_RccpArm.cpp")
source("X_Functions.R")

#### PRE-SETS ####

## Create colours for the plot
pop.groups <- c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150)
my.colours <- "PuBu"

model.name <- "ningaloo"

#### LOAD FILES ####

## Normal Model Files
setwd(sg_dir)
AdultMove <- readRDS(paste0(model.name, sep="_", "movement"))
Effort <- readRDS(paste0(model.name, sep="_", "fishing"))
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
Water <- readRDS(paste0(model.name, sep="_","water"))
YearlyTotal <- readRDS(paste0(model.name, sep="_", "BurnInPop"))
Selectivity <- readRDS("selret")
Mature <- readRDS("maturity")
Weight <- readRDS("weight")

BurnInPop <- YearlyTotal

Settlement <- readRDS(paste0(model.name, sep="_","recruitment")) 
Settlement <- as.vector(Settlement[,2]) # At some point should take this out and put it into the movement file and save it correctly

## Simulation Files
## Need to set different seeds for each scenario so that they are different but run the same every time 
# setwd(sim_dir)
# Effort <- readRDS(paste0(model.name, sep="_", "S03_fishing"))

#### PARAMETER VALUES ####
## Natural Mortality
# We have instantaneous mortality from Marriott et al 2011 and we need to convert that into monthly mortality
NatMort = 0.146

## Potentially 18% mortality rate per day for larvae https://www.nature.com/articles/s41598-022-25629-w
## Larval duration is thought to be ~45 days https://www.researchgate.net/publication/221894788_Understanding_age-specific_dispersal_in_fishes_through_hydrodynamic_modelling_genetic_simulations_and_microsatellite_DNA_analysis#pf9


# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
BHa = 0.4344209 #0.4344209
BHb =	0.01889882
PF = 0.5

# Model settings
MaxCell = nrow(Water) # Number of cells in the model
MaxAge = 29 # The max age of the fish in the model -1 to account for the fact that Rcpp functions start from 0
MaxYear = 58 # Number of years the model should run for -1 to account for the fact that Rcpp functions start from 0
# PlotTotal <- T #This is whether you want a line plot of the total or the map

Pop.Groups <- seq(1,12)

#### SET UP LISTS TO HOLD THE PLOTS ####
SpatialPlots <- list()
LengthPlots <- list()
AgePlots <- list()
TimesPlotted = 0

#### SET UP INITIAL POPULATION ####
PopTotal <- array(0, dim=c(MaxCell, 12, MaxYear+1)) # This is our total population, all ages are summed and each column is a month (each layer is a year)
Total <- array(NA, dim=c(MaxYear+1,1)) # For plotting
bio.catch <- array(0, dim=c(MaxCell, MaxAge+1))
yearly.catch <- array(0, dim=(c(MaxYear+1, 1)))

Sim_Pop <- array(0, dim=c(MaxYear+1, 100))
Sim_Catches <- array(0, dim=c(MaxYear+1, 100))
Sim_Ages <- array(0, dim=c(MaxYear+1, MaxAge+1, 100))

Sim_YearlyTotals <- list()

#### RUN MODEL ####
BurnIn = F #This is to swap the model between burn in and running the model properly
setwd(pop_dir)

RUN <- 0

Start=Sys.time()
for (SIM in 1:100){
  
  #### SET UP LISTS TO HOLD THE PLOTS ####
  SpatialPlots <- list()
  LengthPlots <- list()
  AgePlots <- list()
  TimesPlotted = 0
  
  #### SET UP INITIAL POPULATION ####
  YearlyTotal <- BurnInPop
  PopTotal <- array(0, dim=c(MaxCell, 12, MaxYear+1)) # This is our total population, all ages are summed and each column is a month (each layer is a year)
  Total <- array(NA, dim=c(MaxYear+1,1)) # For plotting
  yearly.catch <- array(0, dim=(c(MaxYear+1, 1)))
  
  print(paste0("SIM", sep=" ", SIM))
  
  for (YEAR in 0:MaxYear){
    
    print(YEAR)
    RUN <- RUN+1
    
    ## Loop over all the Rcpp functions in the model
    ModelOutput <- RunModelfunc_cpp(YEAR, MaxAge, MaxYear, MaxCell, NatMort, BHa, BHb, PF, AdultMove, Mature, Weight, Settlement, 
                                    YearlyTotal, Selectivity, Effort)
    
    ## Save some outputs from the model 
    # Have to add 1 to all YEAR because the loop is now starting at 0
    PopTotal[ , , YEAR+1] <- rowSums(ModelOutput$YearlyTotal[,,1:30], dim=2) # This flattens the matrix to give you the number of fish present in the population each month in each cell, with layers representing the year
    
    total.catch <- ModelOutput$Month_catch # This gives you catch by weight in each cell for each month, which each layer representing an age class
    bio.catch <- colSums(total.catch[,,1:MaxAge], dim=2) # Biomass of fish caught in each age group (there are no fish caught in age 30 because at this point they would be dead)

    yearly.catch[YEAR+1,1] <- sum(bio.catch)
    
    Water$pop <- PopTotal[ , 12, YEAR+1] # We just want the population at the end of the year
    
    ## Plotting ##
    Total[YEAR+1,1] <- sum(Water$pop)
    #print(Total[YEAR+1,1])
    
    if(BurnIn==F & YEAR==58|BurnIn==T & YEAR==58){ # Have to subtract one because the loop now starts at 0
      TotalPop <- as.data.frame(Total) %>%
        rename(Tot.Pop="V1")
      TotalPop$Year <- seq(1960, 2018, by=1)
      TotalPlot <- total.plot.func(pop=TotalPop)
      print(TotalPlot)
    } else { }
    
    if(BurnIn==F & YEAR %%5==0|BurnIn==F & YEAR==58){
      TimesPlotted <- TimesPlotted+1
      SpatialPlots[[TimesPlotted]] <- spatial.plot.func(area=Water, pop=Total, pop.breaks=pop.groups, colours="PuBu")
      AgePlots[[TimesPlotted]] <- age.plot.func(pop=YearlyTotal, NTZs=NoTake)
      #LengthPlots[[TimesPlotted]] <- length.plot.func()
      
    } else { }
    
    Sim_Pop[YEAR+1, SIM] <- Total[YEAR+1,1]
    Sim_Ages[YEAR+1, ,SIM] <- colSums(ModelOutput$YearlyTotal[,12,1:30])
    Sim_Catches[YEAR+1,SIM] <- yearly.catch[YEAR+1,1]
    
    Sim_YearlyTotals[[RUN]] <- ModelOutput["YearlyTotal"]
    
    # filename <- paste0(model.name, sep="_", "Rcpp_YearlyTotal", sep="_", YEAR)
    # saveRDS(ModelOutput["YearlyTotal"], file=filename)

    # filename <- paste0(model.name, sep="_", "Rcpp_BHRecruits", sep="_", YEAR)
    # saveRDS(ModelOutput["BH_recs"], file=filename)
    
  }
  
  ## NEED TO FIGURE OUT A WAY TO SAVE THE CELL REFERENCES AS WELL ##
  if(SIM==100){
    setwd(pop_dir)
    
    filename <- paste0(model.name, sep="_", "All_Yearly_Population", sep="_", "normal")
    saveRDS(Sim_YearlyTotals, file=filename)
    
    filename <- paste0(model.name, sep="_", "Total_Population", sep="_", "normal")
    saveRDS(Sim_Pop, file=filename)
    
    filename <- paste0(model.name, sep="_", "Age_Distribution", sep="_", "normal")
    saveRDS(Sim_Ages, file=filename)
    
    filename <- paste0(model.name, sep="_", "Yearly_Catch", sep="_", "normal")
    saveRDS(Sim_Catches, file=filename)
  }

}
End = Sys.time() 
Runtime = End - Start
Runtime

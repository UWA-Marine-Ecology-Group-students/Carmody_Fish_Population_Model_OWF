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

rm(list = ls())
#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")
pop_dir <-  paste(working.dir, "Output_Population", sep="/")
sim_dir <- paste(working.dir, "Simulations", sep="/")

#### PRE-SETS ####

model.name <- "small"

#### LOAD FILES ####

## Normal Model Files
setwd(sg_dir)
movement <- readRDS(paste0(model.name, sep="_", "movement"))
juv_movement <- readRDS(paste0(model.name, sep="_","movement")) # Have taken out having juvenile movement and now it's the same as adult but I don't want to change the code for the functions  
recruitment <- readRDS(paste0(model.name, sep="_","recruitment"))
fishing <- readRDS(paste0(model.name, sep="_", "fishing"))
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
water <- readRDS(paste0(model.name, sep="_","water"))
YearlyTotal <- readRDS(paste0(model.name, sep="_", "BurnInPop"))
selectivity <- readRDS("selret")
maturity <- readRDS("maturity")
weight <- readRDS("weight")

## Simulation Files
setwd(sim_dir)
fishing <- readRDS(paste0(model.name, sep="_", "S03_fishing"))


## Read in functions
setwd(working.dir)
source("X_Functions.R")

#### PARAMETER VALUES ####
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

#### SET UP LISTS TO HOLD THE PLOTS ####
SpatialPlots <- list()
LengthPlots <- list()
AgePlots <- list()
TimesPlotted <- 0

#### SET UP INITIAL POPULATION ####
PopTotal <- array(0, dim=c(NCELL, 12, length(Time))) # This is our total population, all ages are summed and each column is a month (each layer is a year)
Total <- array(NA, dim=c(length(Time),1)) # For plotting
bio.catch <- array(0, dim=c(NCELL, length(Ages)))
monthly.catch <- array(0, dim=c(length(Ages), 12))
yearly.catch <- array(0, dim=(c(length(Time), 3)))

#### RUN MODEL ####
BurnIn = F #This is to swap the model between burn in and running the model properly
setwd(sim_dir)

for(YEAR in 1:length(Time)){
  
  for(MONTH in 1:12){
    
    ## Movement
    for(A in 1:dim(YearlyTotal)[3]){

      YearlyTotal[ , MONTH, A] <- movement.func(Age=A, Month=MONTH, Population=YearlyTotal, Max.Cell=NCELL, Adult.Move= movement,
                                                Juv.Move=juv_movement)

    }  # End bracket for movement

    ## Mortality
    
    for(A in 1:dim(YearlyTotal)[3]){
      
      if(MONTH==12 & 2<=A & A<30){
        survived.catch <- mortality.func(Age=A, Nat.Mort=M, Effort=fishing, Max.Cell = NCELL,
                                                Month=MONTH, Select=selectivity, Population=YearlyTotal, Year=YEAR)
      
        YearlyTotal[ ,1, A+1] <- survived.catch[[1]]
       
         # Calculate catch
        n.catch <- survived.catch[[2]]
        
        bio.catch[ ,A] <- n.catch * weight[(A*12)+1]
        
      } else if (MONTH!=12) {
        survived.catch <- mortality.func(Age=A, Nat.Mort=M, Effort=fishing, Max.Cell = NCELL,
                                         Month=MONTH, Select=selectivity, Population=YearlyTotal, Year=YEAR)
        
        YearlyTotal[ ,MONTH+1, A] <- survived.catch[[1]]
        
        # Calculate catch
        n.catch <- survived.catch[[2]] # Fish caught in each cell of one age in one month
        
        bio.catch[,A] <- n.catch * weight[(A*12)+1]
      
      } else { }
      
      
    } # End Mortality
    
    ## Recruitment
    if(MONTH==10){
      Recruits <- recruitment.func(Population=YearlyTotal, settlement=recruitment, Max.Cell=NCELL, 
                                   BHa=alpha, BHb=beta, Mature=maturity, Weight=weight, PF=0.5)
      
      YearlyTotal[ ,1,1] <- Recruits
      
    } else { }
    # End Recruitment
    
    monthly.catch[1:30,MONTH] <- colSums(bio.catch)
    
  } #End bracket for months
  
  PopTotal[ , , YEAR] <- rowSums(YearlyTotal[,,Ages], dim=2) # This flattens the matrix to give you the number of fish present in the population each month, with layers representing the years
  
  yearly.catch[YEAR,3] <- sum(monthly.catch)
  
  print(YEAR)
  water$pop <- PopTotal[ , 12, YEAR] # We just want the population at the end of the year
  
  ## Plotting ##
  Total[YEAR,1] <- sum(water$pop)
  print(Total[YEAR,1])
  
  if(BurnIn==F & YEAR==59|BurnIn==T & YEAR==50){
    TotalPop <- as.data.frame(Total) %>% 
      rename(Tot.Pop="V1")
    TotalPop$Year <- seq(1960, 2018, by=1)
    TotalPlot <- total.plot.func(pop=TotalPop) 
    print(TotalPlot)
  } else { }
  
  if(BurnIn==F & YEAR %%5==0|BurnIn==F & YEAR==59){
    TimesPlotted <- TimesPlotted+1
    SpatialPlots[[TimesPlotted]] <- spatial.plot.func(area=water, pop=Total, pop.breaks=pop.groups, colours="PuBu")
    AgePlots[[TimesPlotted]] <- age.plot.func(pop=YearlyTotal, NTZs=NoTake)
    #LengthPlots[[TimesPlotted]] <- length.plot.func()
    
  } else { }
  
  filename <- paste0(model.name, sep="_", "S03_YearlyTotal", sep="_", YEAR)
  saveRDS(YearlyTotal, file=filename)
  
  Sys.sleep(0.3)
}


## SAVE THE CATCHES AND GET SPATIAL PLOT ##

saveRDS(yearly.catch, file=paste0(model.name, sep="_", "S03_YearlyCatch"))

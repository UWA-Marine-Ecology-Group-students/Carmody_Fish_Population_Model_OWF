###################################################

# Script just for running the full model 
# Requires files that are made in the first script
# Will run the full model with fishing
# Make sure you reset things like NCELL if you've
# used the test script

###################################################
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(stringr)
library(forcats)
library(RColorBrewer)
library(geosphere)


#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
m_dir <- paste(working.dir, "Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")

#### PRE-SETS ####

## Create colours for the plot
pop.groups <- c(0,50,100,500,1000,2000,3000,4000,5000,6000)
my.colours <- "RdBu"

## Read in functions
setwd(working.dir)
source("X_Functions.R")
source("03_Population-Set-Up.R")

#### LOAD FILES ####
setwd(sg_dir)
movement <- readRDS("movement")
juv_movement <- readRDS("juvmove")
recruitment <- readRDS("recruitment")
fishing <- readRDS("Fishing")
NoTake <- readRDS("NoTake")
water <- readRDS("water")
selectivity <- readRDS("selectivity")
maturity <- readRDS("maturity")
weight <- readRDS("weight")
start.pop <- readRDS("Starting_Pop")

#### PARAMETER VALUES ####

## Natural Mortality
# We have instantaneous mortality from Marriot et al 2011 and we need to convert that into monthly mortality
M <- 0.146
step <- 1/12 # We're doing a monthly timestep here

# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
alpha <- 0.3245958
beta <- 0.0001910434

#Fish movement parameters
SwimSpeed <- 1.0 # Swim 1km in a day - this is completely made up 

NCELL <- nrow(water)
Ages <- seq(1,30) #These are the ages you want to plot 
Time <- seq(1,100) #This is how long you want the model to run for
PlotTotal <- T #This is whether you want a line plot of the total or the map

Pop.Groups <- seq(1,12)

#### SET UP INITIAL POPULATION ####

Total <- array(NA, dim=c(length(Time),1))

YearlyTotal <- array(0, dim = c(NCELL,12,30)) #This is our yearly population split by age category (every layer is an age group)

start.pop.year <- start.pop %>% 
  slice(which(row_number() %% 12 == 1)) # Gives you the total in each age group at the end of the year

for(d in 1:dim(YearlyTotal)[3]){ # This allocates fish to cells randomly with the fish in age group summing to the total we caluclated above - beware the numbers change slightly due to rounding errors
  for(N in 1:start.pop.year[d,1]){
    cellID <- ceiling((runif(n=1, min = 0, max = 1))*NCELL)
    YearlyTotal[cellID,1,d] <- YearlyTotal[cellID,1,d]+1
  }
}

PopTotal <- array(0, dim=c(NCELL, 12, length(Time))) # This is our total population, all ages are summed and each column is a month (each layer is a year)

#### RUN MODEL ####
BurnIn = T #This is to swap the model between burn in and running the model properly

for(YEAR in 1:length(Time)){
  
  for(MONTH in 1:12){
    
    ## Movement - this is where you'd change things to suit the months
    for(A in 1:dim(YearlyTotal)[3]){
      
      YearlyTotal[ , MONTH, A] <- movement.func(Age=A, Month=MONTH, Population=YearlyTotal, Max.Cell=NCELL, Adult.Move= movement,
                                                Juv.Move=juv_movement)
      # }
      
    }  # End bracket for movement
    
    ## Mortality
    
    for(A in 1:dim(YearlyTotal)[3]){
      
      if(MONTH==12 & 2<=A & A<30){
        YearlyTotal[ ,1, A+1] <- mortality.func(Age=A, mort.50=A50, mort.95=A95, Nat.Mort=M, NTZ=NoTake, Effort=fishing, Max.Cell = NCELL,
                                                Month=MONTH, Select=Selectivity, Population=YearlyTotal, Year=YEAR)
      } else if (MONTH!=12) {
        YearlyTotal[ ,MONTH+1, A] <-mortality.func(Age=A, mort.50=A50, mort.95=A95, Nat.Mort=M, NTZ=NoTake, Effort=fishing, Max.Cell = NCELL,
                                                   Month=MONTH, Select=Selectivity, Population=YearlyTotal, Year=YEAR)
      } else { }
      
    } # End Mortality
    
    ## Recruitment
    if(MONTH==10){
      Recruits <- recruitment.func(Population=YearlyTotal, mat.95=M95, mat.50=M50, settlement=recruitment, 
                                   Max.Cell=NCELL, BHa=alpha, BHb=beta, Mature=maturity, Weight=weight, PF=0.5)
      
      YearlyTotal[ ,1,1] <- Recruits
      
    } else { }
    # End Recruitment
  } #End bracket for months
  
  PopTotal[ , , YEAR] <- rowSums(YearlyTotal[,,Ages]) # This flattens the matrix to give you the number of fish present in the population each month, with layers representing the ages
  
  
  print(YEAR)
  water$pop <- PopTotal[ , 12, YEAR] # We just want the population at the end of the year
  
  ## Plotting ##
  Total[YEAR,1] <- sum(water$pop)
  
  plots <- plotting.func(area=water, pop=Total, pop.breaks=pop.groups, colours="RdBu")
  print(plots)
  
  Sys.sleep(3)
}

## Save burn in population for use in the actual model
setwd(sg_dir)
saveRDS(YearlyTotal, file="BurnInPop")

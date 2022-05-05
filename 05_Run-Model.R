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

#### LOAD FILES ####
setwd(sg_dir)
movement <- readRDS("movement")
juv_movement <- readRDS("juvmove")
recruitment <- readRDS("recruitment")
fishing <- readRDS("fishing")
NoTake <- readRDS("NoTake")
water <- readRDS("water")

#### PARAMETER VALUES ####

## Natural Mortality
# We have instantaneous mortality from Marriot et al 2011 and we need to convert that into monthly mortality
yearly_surv=exp(-0.146)
monthly_mort=1-(yearly_surv^(1/12))

M <- monthly_mort # Natural mortality rate per month

#Ricker recruitment model parameters (these are currently just made up values)
a <- 7
b <- 0.0017
M50 <- 2 # From Grandcourt et al. 2010
M95 <- 5 # From Grandcourt et al. 2010 (technically M100)
relationship <- 0.76

#Fish movement parameters
SwimSpeed <- 1.0 # Swim 1km in a day - this is completely made up 

## Fishing mortality parameters
A50 <- 4 # For L. miniatus from Williams et al. 2010 # L. miniatus becomes vulnerable to fishing at about age two
A95 <- 6 # For L. miniatus from Williams et al. 2012
q <- 0.5 # Apparently this is what lots of stock assessments set q to be...

NCELL <- nrow(water)
Ages <- seq(4,8) #These are the ages you want to plot 
Time <- seq(1,500) #This is how long you want the model to run for
PlotTotal <- F #This is whether you want a line plot of the total or the map

Pop.Groups <- seq(1,12)

#### RUN MODEL ####

BurnIn = F #This is to swap the model between burn in and running the model properly
Total <- array(NA, dim=c(length(Time),1))

YearlyTotal <- array(0, dim = c(NCELL,12,30)) #This is our yearly population split by age category (every layer is an age group)
# If you change age you have to change it in the fish mortality function too
for(d in 1:dim(YearlyTotal)[3]){
  YearlyTotal[,1,d] <- matrix(floor(runif(NCELL, 1, 100))) #50 is too few
}

PopTotal <- array(0, dim=c(NCELL, 12, length(Time))) # This is our total population, all ages are summed and each column is a month (each layer is a year)

for(YEAR in 1:length(Time)){
  
  for(MONTH in 2:13){
    
    ## Movement - this is where you'd change things to suit the months
    for(A in 2:dim(YearlyTotal)[3]){
      
      if(MONTH==13){
        YearlyTotal[ , 12, A-1] <- movement.func(Age=A, Month=MONTH, Population=YearlyTotal, Max.Cell=NCELL, Adult.Move= movement,
                                                 Juv.Move=juv_movement)
      } else {
        YearlyTotal[ , MONTH, A-1] <- movement.func(Age=A, Month=MONTH, Population=YearlyTotal, Max.Cell=NCELL, Adult.Move= movement,
                                                    Juv.Move=juv_movement)
      } 
      
    }  # End bracket for movement
    
    ## Fishing Mortality
    
    for(A in 1:dim(YearlyTotal)[3]){
      
      #YearPop <- YearlyTotal[ , ,A]
      
      YearlyTotal[ ,MONTH-1 ,A] <- mortality.func(Age=A, mort.50=A50, mort.95=A95, Nat.Mort=M, NTZ=NoTake, Effort=fishing, Cell=CELL, Max.Cell = NCELL,
                                                  Month=MONTH, Year=YEAR, Population=YearlyTotal)
      
    } # End Mortality
    
    ## Recruitment
    
    if(MONTH==11){
      YearlyTotal[,1,1] <- recruitment.func(Population=YearlyTotal, Age=A, mat.95=M95, mat.50=M50, settlement=recruitment, 
                                            Max.Cell=NCELL, relationship=0.76)
    } else { } 
    # End Recruitment
  } #End bracket for months
  
  PopTotal[ , , YEAR] <- rowSums(YearlyTotal[,,Ages]) # This flattens the matrix to give you the number of fish present in the population each month, with layers representing the years
  # PopTotal[ , , YEAR] <- YearlyTotal[ , , 1]  # To look at specific parts of the matrix for the plot
  
  
  print(YEAR)
  water$pop <- PopTotal[ , 12, YEAR] # We just want the population at the end of the year
  
  ## Plotting ##
  Total[YEAR,1] <- sum(water$pop)
  
  plots <- plotting.func(area=water, pop=Total, pop.breaks=pop.groups, colours="RdBu")
  print(plots)
  
  Sys.sleep(3)
}

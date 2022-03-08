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

## Create colours for the plots
pop.groups <- c(0,500,1000,5000,10000,20000,30000,40000,50000,60000,
                70000,80000,90000,100000,200000,300000,400000,500000,600000)
my.colours <- "RdBu"

## Read in functions
setwd(working.dir)
source("X_Functions.R")

#### LOAD FILES ####
setwd(sg_dir)
movement <- readRDS("movement")
juv_movement <- readRDS("juvmove")
fishing <- readRDS("Fishing")
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

#Fish movement parameters
SwimSpeed <- 1.0 # Swim 1km in a day - this is completely made up 

## Fishing mortality parameters
A50 <- 4 # For L. miniatus from Williams et al. 2010 # L. miniatus becomes vulnerable to fishing at about age two
A95 <- 6 # For L. miniatus from Williams et al. 2012
q <- 0.5 # Apparently this is what lots of stock assessments set q to be...

#### RUN MODEL ####

# Fill the array with random numbers of fish for now
# Need to change the loop so that rather than filling in a column it adds a new column each time 

YearlyTotal <- array(0, dim = c(NCELL,12,100)) #This is our yearly population split by age category (every layer is an age group)
for(d in 1:100){
  YearlyTotal[,1,] <- matrix(floor(runif(NCELL, 1, 2000)))
}

PopTotal <- array(0, dim=c(NCELL, 12, 59)) # This is our total population, all ages are summed and each column is a month (each layer is a year)

Movement <- NULL

NCELL <- nrow(water)

# The time steps now represent 1 month (need to make sure you change the distance your fish can swim to match the timestep)

for(YEAR in 2:59){
  
  for(MONTH in 2:13){
    
      for(A in 2:dim(YearlyTotal)[3]){
        
        movement.func(Age=A, Month=MONTH, Population=YearlyTotal, Max.Cell=NCELL, Adult.Move= movement, Juv.Move=juv_movement)

      } #End bracket for movement 
    
      ## Fishing Mortality
      for(A in 1:dim(YearlyTotal)[3]){
        
        mortality.func(Age=A, mort.50=A50, mort.95=A95, Nat.Mort=M, NTZ=NoTake, Effort=Fishing, Cell=CELL, Max.Cell = NCELL,
                       Month=MONTH, Year=YEAR, Population=YearlyTotal)
      } # End Mortality
    
    # ## Recruitment
    # Recs <- matrix(0, nrow=1, ncol=1) #Blank matrix to add the recruits to
    # Recruitment <- as.matrix(colSums(Pop[,t,], dims=1)) 
    # for(A in 1:dim(Recruitment)[1]){
    #   Mature <- 1/(1+(exp(-log(19)*((A-M95)/(M95-M50))))) #Number of mature individuals in each age class
    #   S <- colSums(Recruitment)*Mature #Spawning stock
    #   Rec <- a*S*exp(-b*S) #Number of recruits from that age class
    #   Recs <- rbind(Recs, Rec) #Combine recruits from each age class into one dataframe
    # }
    
    PopTotal[ , , YEAR] <- rowSums(YearlyTotal) 
    # Collapses the 3D matrix and gives you the number of fish of all ages in each cell in each month
    
  }
  print(YEAR)
  water$pop <- PopTotal[ , 12, YEAR] # We just want the population at the end of the year
  
  ## Plotting ##
  map <- plotting.func(area=water, pop.breaks=pop.groups, colours="RdBu")
  print(map)
  
  Sys.sleep(3)
}


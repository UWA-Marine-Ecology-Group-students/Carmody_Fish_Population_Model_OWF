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
library(forcats)
library(RColorBrewer)
library(MQMF)
library(Rcpp)
library(RcppArmadillo)
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

## Read in functions
setwd(working.dir)
sourceCpp("X_Model_RccpArm.cpp")
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
Effort <- readRDS(paste0(model.name, sep="_", "burn_in_fishing"))
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
Water <- readRDS(paste0(model.name, sep="_","water"))
start.pop <- readRDS(paste0(model.name, sep="_","Starting_Pop"))
Selectivity <- readRDS("selret")
Mature <- readRDS("maturity")
maturity <- readRDS("maturity")
Weight <- readRDS("weight")

Settlement <- readRDS(paste0(model.name, sep="_","recruitment")) 

## Fix selectivity so you have the most up to date selectivity and you have enough layers for the model to run
Selectivity <- Selectivity[,,45:59]

for(i in 1:3){
  Selectivity <- abind(Selectivity, Selectivity, along=3)
}


setwd(sg_dir)
water <- readRDS(paste0(model.name, sep="_","water"))

#### PARAMETER VALUES ####
## Natural Mortality
# We have instantaneous mortality from Marriott et al 2011 and we need to convert that into monthly mortality
NatMort = 0.146

# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
BHa = 0.4344209 
BHb = 0.003132717
PF = 0.5

# Model settings
MaxCell = nrow(Water) # Number of cells in the model
MaxAge = 30 # The max age of the fish in the model -1 to account for the fact that Rcpp functions start from 0
MaxYear = 50 # Number of years the model should run for 
# PlotTotal <- T #This is whether you want a line plot of the total or the map

Pop.Groups <- seq(1,12)

#### SET UP LISTS TO HOLD THE PLOTS ####
SpatialPlots <- list()
LengthPlots <- list()
AgePlots <- list()
TimesPlotted = 0

#### SET UP INITIAL POPULATION ####
Total <- array(0, dim=c(MaxYear+1,1))

YearlyTotal <- array(0, dim = c(MaxCell,12,30)) #This is our yearly population split by age category (every layer is an age group)

start.pop.year <- start.pop %>% 
  slice(which(row_number() %% 12 == 1)) # Gives you the total in each age group at the end of the year

for(d in 1:dim(YearlyTotal)[3]){ # This allocates fish to cells randomly with the fish in age group summing to the total we caluclated above - beware the numbers change slightly due to rounding errors
  for(N in 1:start.pop.year[d,1]){
    cellID <- ceiling((runif(n=1, min = 0, max = 1))*MaxCell)
    YearlyTotal[cellID,1,d] <- YearlyTotal[cellID,1,d]+1
  }
}

PopTotal <- array(0, dim=c(MaxCell, 12, MaxYear)) # This is our total population, all ages are summed and each column is a month (each layer is a year)

#### RUN MODEL ####

Start=Sys.time()
for (YEAR in 0:49){
  
  print(YEAR)
  
  ## Loop over all the Rcpp functions in the model
  ModelOutput <- RunModelfunc_cpp(YEAR, MaxAge, MaxYear, MaxCell, NatMort, BHa, BHb, PF, AdultMove, Mature, Weight, Settlement, 
                                  YearlyTotal, Selectivity, Effort)
  
  ## Save some outputs from the model 
  # Have to add 1 to all YEAR because the loop is now starting at 0
  PopTotal[ , , YEAR+1] <- rowSums(ModelOutput$YearlyTotal[,,1:MaxAge], dim=2) # This flattens the matrix to give you the number of fish present in the population each month, with layers representing the years
  
  water$pop <- PopTotal[ , 12, YEAR+1] # We just want the population at the end of the year
  
  Total[YEAR+1,1] <- sum(water$pop)
  #print(Total[YEAR+1,1])
 

}

End=Sys.time()
Runtime = End - Start
Runtime

## PLOTTING
Total <- as.data.frame(Total)
plot(x=seq(1,51,1), y=Total$V1)

## Save burn in population for use in the actual model
setwd(sg_dir)
saveRDS(YearlyTotal, file=paste0(model.name, sep="_", "BurnInPop"))


#### RUN BURN IN AGAIN WITH HIGH LEVEL OF FISHING MORTALITY ####
#### SET UP LISTS TO HOLD THE PLOTS ####
SpatialPlots <- list()
LengthPlots <- list()
AgePlots <- list()
TimesPlotted = 0

#### SET UP INITIAL POPULATION ####
Total <- array(0, dim=c(MaxYear+1,1))

setwd(sg_dir)
YearlyTotal <- readRDS("ningaloo_BurnInPop") #This is our yearly population split by age category (every layer is an age group)
Effort <- readRDS("ningaloo_burn_in_fishing_High_M")

PopTotal <- array(0, dim=c(MaxCell, 12, MaxYear)) # This is our total population, all ages are summed and each column is a month (each layer is a year)

#### RUN MODEL ####

Start=Sys.time()
for (YEAR in 0:30){
  
  print(YEAR)
  
  ## Loop over all the Rcpp functions in the model
  ModelOutput <- RunModelfunc_cpp(YEAR, MaxAge, MaxYear, MaxCell, NatMort, BHa, BHb, PF, AdultMove, Mature, Weight, Settlement, 
                                  YearlyTotal, Selectivity, Effort)
  
  ## Save some outputs from the model 
  # Have to add 1 to all YEAR because the loop is now starting at 0
  PopTotal[ , , YEAR+1] <- rowSums(ModelOutput$YearlyTotal[,,1:MaxAge], dim=2) # This flattens the matrix to give you the number of fish present in the population each month, with layers representing the years
  
  water$pop <- PopTotal[ , 12, YEAR+1] # We just want the population at the end of the year
  
  Total[YEAR+1,1] <- sum(water$pop)
  #print(Total[YEAR+1,1])
  
  
}

End=Sys.time()
Runtime = End - Start
Runtime

## PLOTTING
Total <- as.data.frame(Total)
plot(x=seq(1,51,1), y=Total$V1)

setwd(sg_dir)
saveRDS(YearlyTotal, file=paste0(model.name, sep="_", "BurnInPop_High_M"))


setwd(sg_dir)
unfished.bio <- readRDS(paste0(model.name, sep="_", "BurnInPop"))

temp <- unfished.bio[,12,]
temp2 <- colSums(temp)
temp3 <- temp2 * maturity[,12]
Unf.Mat.Bio <- sum(temp3 * Weight[,12])


fished.bio <- YearlyTotal
temp <- fished.bio[,12,]
temp2 <- colSums(temp)
temp3 <- temp2 * maturity[,12]
Mat.Bio <- sum(temp3 * Weight[,12])


Mat.Bio/Unf.Mat.Bio

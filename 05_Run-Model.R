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
pop_dir <-  paste(working.dir, "Output_Population", sep="/")

#### PRE-SETS ####

## Create colours for the plot
pop.groups <- c(0,10,20,30,40,50,60,70,80,90,100,110,120,130)
my.colours <- "PuBu"

## Read in functions
setwd(working.dir)
source("X_Functions.R")

#### LOAD FILES ####
setwd(sg_dir)
movement <- readRDS("movement")
juv_movement <- readRDS("juvmove")
recruitment <- readRDS("recruitment")
fishing <- readRDS("fishing")
NoTake <- readRDS("NoTakeList")
water <- readRDS("water")
selectivity <- readRDS("selret")
maturity <- readRDS("maturity")
weight <- readRDS("weight")
YearlyTotal <- readRDS("BurnInPop")

#### PARAMETER VALUES ####
## Natural Mortality
# We have instantaneous mortality from Marriott et al 2011 and we need to convert that into monthly mortality
M <- 0.146
step <- 1/12 # We're doing a monthly time step here

# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
alpha <- 0.3245958
beta <- 0.0001910434

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

#### RUN MODEL ####
BurnIn = F #This is to swap the model between burn in and running the model properly
setwd(pop_dir)

for(YEAR in 1:length(Time)){
  
  for(MONTH in 1:12){
    
    ## Movement - this is where you'd change things to suit the months
    for(A in 1:dim(YearlyTotal)[3]){
      
      YearlyTotal[ , MONTH, A] <- movement.func(Age=A, Month=MONTH, Population=YearlyTotal, Max.Cell=NCELL, Adult.Move= movement,
                                                Juv.Move=juv_movement)
      
    }  # End bracket for movement
    
    ## Mortality
    
    for(A in 1:dim(YearlyTotal)[3]){
      
      if(MONTH==12 & 2<=A & A<30){
        YearlyTotal[ ,1, A+1] <- mortality.func(Age=A, Nat.Mort=M, Effort=fishing, Max.Cell = NCELL,
                                                Month=MONTH, Select=selectivity, Population=YearlyTotal, Year=YEAR)
      } else if (MONTH!=12) {
        YearlyTotal[ ,MONTH+1, A] <- mortality.func(Age=A, Nat.Mort=M, Effort=fishing, Max.Cell = NCELL,
                                                   Month=MONTH, Select=selectivity, Population=YearlyTotal, Year=YEAR)
      
      } else { }
      
    } # End Mortality
    
    ## Recruitment
    if(MONTH==10){
      Recruits <- recruitment.func(Population=YearlyTotal, settlement=recruitment, Max.Cell=NCELL, 
                                   BHa=alpha, BHb=beta, Mature=maturity, Weight=weight, PF=0.5)
      
      YearlyTotal[ ,1,1] <- Recruits
      
    } else { }
    # End Recruitment
  } #End bracket for months
  
  PopTotal[ , , YEAR] <- rowSums(YearlyTotal[,,Ages], dim=2) # This flattens the matrix to give you the number of fish present in the population each month, with layers representing the years
  
  
  print(YEAR)
  water$pop <- PopTotal[ , 12, YEAR] # We just want the population at the end of the year
  
  ## Plotting ##
  Total[YEAR,1] <- sum(water$pop)
  print(Total[YEAR,1])
  
  if(BurnIn==F & YEAR==59|BurnIn==T & YEAR==100){
    Total <- as.data.frame(Total)
    Total$Year <- seq(1960, 2018, by=1)
    TotalPlot <- total.plot.func(pop=Total)
    print(TotalPlot)
  } else { }
  
  if(BurnIn==F & YEAR %%5==0|BurnIn==F & YEAR==59){
    TimesPlotted <- TimesPlotted+1
    SpatialPlots[[TimesPlotted]] <- spatial.plot.func(area=water, pop=Total, pop.breaks=pop.groups, colours="PuBu")
    AgePlots[[TimesPlotted]] <- age.plot.func(pop=YearlyTotal, NTZs=NoTake)
    #LengthPlots[[TimesPlotted]] <- length.plot.func()
  } else { }
  
  filename <- paste("YearlyTotal", YEAR, sep=".")
  saveRDS(YearlyTotal, file=filename)
  
  Sys.sleep(3)
}


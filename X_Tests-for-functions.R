###################################################

# Script that has small chunks of code to test that
# the functions are doing what you want
# Creates smaller parts of the model data frames
# with easy numbers so that you can very quickly
# check your function outputs

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
library(Rcpp)
library(RcppArmadillo)


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
sourceCpp("X_Model_RccpArm.cpp")

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
start.pop <- readRDS("Starting_Pop")

#### PARAMETER VALUES ####

## Natural Mortality
# We have instantaneous mortality from Marriot et al 2011 and we need to convert that into monthly mortality
yearly_surv=exp(-0.146)
monthly_mort=1-(yearly_surv^(1/12))

M <- monthly_mort # Natural mortality rate per month

# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
alpha <- 0.3245958
beta <- 0.0001910434

#Fish movement parameters
SwimSpeed <- 1.0 # Swim 1km in a day - this is completely made up 

## Fishing mortality parameters
# A50 <- 4 # For L. miniatus from Williams et al. 2010 # L. miniatus becomes vulnerable to fishing at about age two
# A95 <- 6 # For L. miniatus from Williams et al. 2012
# q <- 0.5 # Apparently this is what lots of stock assessments set q to be...

NCELL <- nrow(water)
Ages <- seq(4,30) #These are the ages you want to plot 
Time <- seq(1,100) #This is how long you want the model to run for
PlotTotal <- T #This is whether you want a line plot of the total or the map

Pop.Groups <- seq(1,12)

#### SET UP A SMALL INITIAL POPULATION FOR THE TESTS ####

Total <- array(NA, dim=c(length(Time),1))

YearlyTotal_test <- array(0, dim = c(NCELL,12,30)) #This is our yearly population split by age category (every layer is an age group)

start.pop.year <- start.pop %>% 
  slice(which(row_number() %% 12 == 1)) # Gives you the total in each age group at the end of the year

for(d in 1:dim(YearlyTotal_test)[3]){ # This allocates fish to cells randomly with the fish in age group summing to the total we caluclated above - beware the numbers change slightly due to rounding errors
  for(N in 1:start.pop.year[d,1]){
    cellID <- ceiling((runif(n=1, min = 0, max = 1))*NCELL)
    YearlyTotal_test[cellID,1,d] <- YearlyTotal_test[cellID,1,d]+1
  }
}

PopTotal <- array(0, dim=c(NCELL, 12, length(Time))) # This is our total population, all ages are summed and each column is a month (each layer is a year)

#### RUN MODEL ####
BurnIn = T #This is to swap the model between burn in and running the model properly

#### Mortality Function Test ####


YearlyTotal_test[ , 2, 1] <- mortality.func(Age=1, mort.50=A50, mort.95=A95, Nat.Mort=M, NTZ=NoTake, Effort=fishing, Cell=CELL, Max.Cell = NCELL,
                                            Month=2, Year=3, Population=YearlyTotal_test)
YearlyTotal_test[ , , 1] 

#### Movement Function Test ####

Movers <- movement.func(Age=2, Month=1, Population=YearlyTotal_test, Max.Cell=NCELL, Adult.Move= movement,
                                            Juv.Move=juv_movement)

sum(YearlyTotal_test[,1,2])
sum(YearlyTotal_test[,2,1]) # should be equal
sum(Movers)
#### Recruitment Function Test ####


YearlyTotal_test[,2,1] <- recruitment.func(Population=YearlyTotal_test, Age=6, mat.95=M95, mat.50=M50, settlement=recruitment, #Normally month would be 1 but for this it's easier to set it at 2
                                      Max.Cell=NCELL, BHa=a, BHb=b, Mature=maturity, Weight=weight, PF=0.5)

YearlyTotal_test[ , , 1] 

water <- water%>%
  mutate(pop = round(pop, digits=0)) %>% 
  mutate(pop_level = cut(pop, pop.breaks)) 

pop.breaks <- pop.groups




##### Testing #####
adults <- YearlyTotal[ ,10, ] %>% 
  colSums(.) # Gives us just females because they are the limiting factor for reproduction
adults <- adults * 0.5
tot.recs <- data.frame(matrix(0, nrow = NCELL, ncol = dim(YearlyTotal)[3]))

tot.recs <- lapply(1:dim(YearlyTotal)[3], function(Age){
  SB <- adults[Age] * maturity[Age,1] * weight[(Age*12)+1] #Gives us spawning biomass in each age group at the end of the year, hence the x 12+1 as it starts at 1 not zero
  TotMatBio <- sum(SB) #Gives us total mature spawning biomass
  recs <- (SB/alpha+beta*SB) # This gives us males and females to go into the next generation
})
tot.recs <- do.call(rbind, tot.recs)

settle.recs <- sum(tot.recs)

settled <- settle.recs*recruitment[,2]



Juv.Pop <- matrix(YearlyTotal[ , 1, 1]) # This gives you the fish in all the sites at time step Month-1 of age A-1

Juv.Pop2 <- sapply(seq(NCELL), function(Cell){
  Juv.Pop2 <- as.matrix(juv_movement[Cell, ] * Juv.Pop[Cell,1]) #This should give you the number of fish that move from each cell to all the other sites
})

Juv.Movement2 <- rowSums(Juv.Pop2)
Juv.Moved <- sum(Juv.Movement2)


All.Movers <- NULL
All.Movers <- cbind(All.Movers, Juv.Movement2)

sum(Juv.Pop)


#### RCPP FUNCTION TESTS ####
setwd(sg_dir)
Water <- readRDS("ningaloo_water")
CellVars <- readRDS("CellVars")
Select <- readRDS("selret")
BoatDays <- readRDS("BR_Trips") %>%  
  arrange(NumYear) 
Spatial_q <- readRDS("ningaloo_Spatial_q_NTZ") # Will need to load in the right file for the scenario
age_survived <- readRDS("age_survived") 

BoatDays <- BoatDays[,4:6] %>% 
  as.data.frame() %>% 
  as.matrix()
  
BoatDays <- array(BoatDays, dim=c(12,4,59))


DistBR <- as.matrix(CellVars[,1:4])
CellArea_log <- log(CellVars[,5])

Fished <- Water[,4:6] %>% 
  st_drop_geometry() %>% 
  unnest()

Fished <- ifelse(Fished == "Y", 1, 0) %>% 
  as.matrix() 
rownames(Fished) <-NULL

FuelPrice = array(1, dim=c(12,1))

MaxAge = 30
NBR = 4
MaxCell = nrow(CellVars)
MONTH = 1
YEAR = 1


### YOU CANNOT CHANGE THE ORDER OF THESE ####
Part1 <- effortfunc_cpp(MaxAge, NBR, MaxCell, MONTH, YEAR,
                          FuelPrice, CellArea_log, Fished,
                          DistBR, Spatial_q,
                          YearlyTotal, Select, BoatDays)

CellU_Rcpp <- Part1[[4]]
RowU_Rcpp <- Part1[[5]]
CellU_NTZ_Rcpp <- Part1[[4]]
Check <- Part1[[6]]
Effort <- Part1[[7]]
TotalEffort2 <- Part1[[8]]
dim(Part1[["CellU"]])
dim(Fished)


MaxAge = 0+1
NBR = 4
MaxCell = nrow(CellVars)
MONTH = 1+1
YEAR = 40+1



Catchable_Age <- matrix(0, nrow=MaxCell, ncol=MaxAge)

for(Cell in 1:MaxCell){
  for(Age in 1:MaxAge){
    Catchable_Age[Cell,Age] <- age_survived[Cell,MONTH,Age]*Select[Age,MONTH,YEAR]
  }
}

Catchable <- rowSums(Catchable_Age)


CellCoef <- array(0, dim=c(MaxCell, NBR))
for(Col in 1:NBR){
  for(Row in 1:MaxCell){
    CellCoef[Row, Col] <- (FuelPrice[1]*DistBR[Row,Col])+Catchable[Row];
  }
}

cellU <- matrix(NA, ncol=4, nrow=MaxCell)

for(RAMP in 1:4){
  for(cell in 1:MaxCell){
    U <- exp(CellCoef[cell,RAMP]+CellArea_log[cell])
    cellU[cell, RAMP] <- U
  }
}

cellU2 <- cellU[ , ] * Fished[,2] 

rowU2 <- as.data.frame(colSums(cellU2))

for (RAMP in 1:4){
  for (cell in 1:NCELL_6086){
    BR_U_6086[cell,RAMP] <- exp(CellCoef[cell,RAMP]+CellArea_log[cell])/rowU[RAMP,1]
  }
}
colSums(BR_U_6086)


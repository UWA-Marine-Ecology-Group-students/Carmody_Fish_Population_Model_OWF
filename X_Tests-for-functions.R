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



#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
m_dir <- paste(working.dir, "Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")

#### PRE-SETS ####

## Create colours for the plot
pop.groups <- c(0,500,1000,5000,10000,20000,30000,40000,50000,60000,
                70000,80000,90000,100000,200000,300000,400000,500000,600000)
my.colours <- "RdBu"

## Read in functions
setwd(working.dir)
source("X_Functions.R")

#### LOAD FILES ####
setwd(sg_dir)
movement <- readRDS("test_movement")
juv_movement <- readRDS("test_juvmove")
recruitment <- readRDS("recruitment")
fishing <- readRDS("test_fishing")
NoTake <- readRDS("test_NoTake")
water <- readRDS("test_water")

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
Fecundity <- 7000000

#Fish movement parameters
SwimSpeed <- 1.0 # Swim 1km in a day - this is completely made up 

## Fishing mortality parameters
A50 <- 4 # For L. miniatus from Williams et al. 2010 # L. miniatus becomes vulnerable to fishing at about age two
A95 <- 6 # For L. miniatus from Williams et al. 2012
q <- 0.5 # Apparently this is what lots of stock assessments set q to be...

NCELL <- nrow(water)
Ages <- seq(4,8) #These are the ages you want to plot 
Time <- seq(1,500) #This is how long you want the model to run for
PlotTotal <- T #This is whether you want a line plot of the total or the map

#### Set up a small population to use for the tests ####
YearlyTotal_test <- array(0, dim = c(NCELL,3,6)) #This is our yearly population split by age category (every layer is an age group)
# If you change age you have to change it in the fish mortality function too
for(d in 1:dim(YearlyTotal_test)[3]){
  YearlyTotal_test[,1,d] <- matrix(floor(runif(NCELL, 1, 100))) #50 is too few
}

PopTotal_test <- array(0, dim=c(NCELL, 12, length(Time))) # This is our total population, all ages are summed and each column is a month (each layer is a year)

#### Mortality Function Test ####

BurnIn = T

YearlyTotal_test[ , 2, 1] <- mortality.func(Age=1, mort.50=A50, mort.95=A95, Nat.Mort=M, NTZ=NoTake, Effort=fishing, Cell=CELL, Max.Cell = NCELL,
                                            Month=2, Year=3, Population=YearlyTotal_test)
YearlyTotal_test[ , , 1] 

#### Movement Function Test ####

YearlyTotal_test[ , 2, 1] <- movement.func(Age=2, Month=2, Population=YearlyTotal_test, Max.Cell=NCELL, Adult.Move= movement,
                                            Juv.Move=juv_movement)

sum(YearlyTotal_test[,1,1])
sum(YearlyTotal_test[,2,1]) # should be equal

#### Recruitment Function Test ####


YearlyTotal_test[,2,1] <- recruitment.func(Population=YearlyTotal_test, Age=6, mat.95=M95, mat.50=M50, settlement=recruitment, #Normally month would be 1 but for this it's easier to set it at 2
                                      Max.Cell=NCELL, Fcd=Fecundity, RRa=a, RRb=b)

YearlyTotal_test[ , , 1] 

water <- water%>%
  mutate(pop = round(pop, digits=0)) %>% 
  mutate(pop_level = cut(pop, pop.breaks)) 

pop.breaks <- pop.groups






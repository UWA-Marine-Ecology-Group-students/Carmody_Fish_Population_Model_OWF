###################################################

# Nothing works and I don't know why

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


#* Normal - RE-DONE ####
setwd(pop_dir)

# Whole Population 
PopTotalMov <- array(0, dim=c(NCELL, 12, 59)) 

numYear <- seq(1,59,1)

for(Y in 1:59){
  
  year <- readRDS(paste0("NM_YearlyTotal_", numYear[Y]))
  year <- rowSums(year[,,1:30], dim=2)
  
  PopTotalMov[,,Y] <- year

}

netmove2 <- array(0, dim=c(59,12))

for (Y in 1:59){
  for(M in 1:12){
    
    if(M==1){
      movement <-  PopTotalMov[,1,Y] - PopTotalMov[,12,Y-1] 
      netmoveNTZ <- sum(movement[c(NoTake[[3]])])
      netmoveF <- sum(movement[-c(NoTake[[3]])])
      
      netmove2[Y,M] <- netmoveNTZ + netmoveF 
      
    } else if(M>1){
      movement <-  PopTotalMov[,M,Y] - PopTotalMov[,M-1,Y]
      netmoveNTZ <- sum(movement[c(NoTake[[3]])])
      netmoveF <- sum(movement[-c(NoTake[[3]])])
      
      netmove2[Y,M] <- netmoveNTZ + netmoveF
      } 
  }
}

## Checking Baranov catch with catch approximation
setwd(pop_dir)
Approximation <- readRDS("ningaloo_Catch_by_Age_S03_medium_movement")
Baranov <- readRDS("ningaloo_Catch_by_Age_Baranov_S03_Medium")

Approximation <-colSums(Approximation[[1]])
Baranov <- colSums(Baranov[[1]])

plot(Approximation)
plot(Baranov)




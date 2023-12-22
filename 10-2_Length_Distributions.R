###################################################

# Script for working out the length ratios
# for different levels of F 

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
library(raster)
library(sfnetworks)
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
msy_dir <- paste(working.dir, "MSY_Outputs", sep="/")

model.name <- "ningaloo"

## Read in functions
setwd(working.dir)
sourceCpp("X_Model_RccpArm.cpp")
source("X_Functions.R")

#### READ IN DATA ####
setwd(sg_dir)
AdultMove <- readRDS(paste0(model.name, sep="_", "movement"))
Settlement <- readRDS(paste0(model.name, sep="_","recruitment")) 
Effort <- readRDS(paste0(model.name, sep="_", "fishing"))
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
water <- readRDS(paste0(model.name, sep="_","water")) %>% 
  st_make_valid()
YearlyTotal <- readRDS(paste0(model.name, sep="_", "BurnInPop"))
Selectivity <- readRDS("selret")
Mature <- readRDS("maturity")
Weight <- readRDS("weight")

#setwd(pop_dir)
#MortRate <- readRDS("Mort_Rate")


NCELL <- nrow(water)

#### READ IN AND FORMAT DATA ####
setwd(msy_dir)

Age.NTZ <- list()
Age.F <- list()

#* No fishing
Age.Dist.NTZ <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S00", sep="_", "no_fishing"))
Age.Dist.F <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S00", sep="_", "no_fishing"))

Age.NTZ[[1]] <- Age.Dist.NTZ[[1]]
Age.F[[1]] <- Age.Dist.F[[1]]


#* F = Fmsy
Age.Dist.NTZ <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S00", sep="_", "Fmsy"))
Age.Dist.F <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S00", sep="_", "Fmsy"))

Age.NTZ[[2]] <- Age.Dist.NTZ[[1]]
Age.F[[2]] <- Age.Dist.F[[1]]

#* F = 2/3M
Age.Dist.NTZ <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S00", sep="_", "two_thirds_M"))
Age.Dist.F <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S00", sep="_", "two_thirds_M"))

Age.NTZ[[3]] <- Age.Dist.NTZ[[1]]
Age.F[[3]] <- Age.Dist.F[[1]]

#* F = M
Age.Dist.NTZ <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S00", sep="_", "NatMort"))
Age.Dist.F <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S00", sep="_", "NatMort"))

Age.NTZ[[4]] <- Age.Dist.NTZ[[1]]
Age.F[[4]] <- Age.Dist.F[[1]]


#* F = 1.5M
Age.Dist.NTZ <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S00", sep="_", "1-5xNatMort"))
Age.Dist.F <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S00", sep="_", "1-5xNatMort"))


Age.NTZ[[5]] <- Age.Dist.NTZ[[1]]
Age.F[[5]] <- Age.Dist.F[[1]]

#### CONVERT AGES INTO LENGTHS ####
## Use equation to work out the larger cut off which is half way between LM and Linf
#Ave.Linf <- (591 + 556)/2 # Average of Linf for males and females

BigLM <- 350 + ((664 - 350)/2)

age.length <- age.to.length(NTZ.ages = Age.NTZ, F.ages = Age.F, max.age = 30, n.scenarios = 5, Linf = 664, k = 0.241, t0 = -0.375, LM = 350, BigLM = 507)

age.length.NTZ <- age.length[[1]]
age.length.F <- age.length[[2]]

## Put everything together into one big dataframe to have a look at it
mortalities <- c("F=0", "F=MSY", "F=2/3M", "F=M", "F=1.5M")
temp.NTZ <- NULL
temp.F <- NULL

for (M in 1:5){
  
  temp1 <- age.length.NTZ[[M]]
  temp1 <- temp1 %>% 
    mutate(Mortality = mortalities[M]) %>% 
    mutate(Zone = "NTZ")
  
  temp.NTZ <- rbind(temp.NTZ, temp1)
  
  temp1 <- age.length.F[[M]]
  temp1 <- temp1 %>% 
    mutate(Mortality = mortalities[M]) %>% 
    mutate(Zone = "F")
  
  temp.F <- rbind(temp.F, temp1)
  
}

age.ratios <- rbind(temp.NTZ, temp.F) %>% 
  filter(!Length.Group %in% "Too Small") %>% 
  group_by(Mortality, Zone, Length.Group) %>% 
  summarise(Total = sum(Count)) %>% 
  ungroup() %>% 
  group_by(Mortality, Zone) %>% 
  summarise(Ratio = Total[Length.Group == "> Big LM"]/Total[Length.Group == "> LM"])






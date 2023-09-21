###################################################

# RUnning GLMs on the model outputs to determine
# which scenarios are different to one another

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
library(gmailr)
library(beepr)

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
AdultMove <- readRDS(paste0(model.name, sep="_", "movement_fast"))
Settlement <- readRDS(paste0(model.name, sep="_","recruitment")) 
Effort <- readRDS(paste0(model.name, sep="_", "fishing"))
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
Water <- readRDS(paste0(model.name, sep="_","water")) %>% 
  st_make_valid()
BurnInPop <- readRDS(paste0(model.name, sep="_", "BurnInPop_High_M"))
Selectivity <- readRDS("selret")
Mature <- readRDS("maturity")
Weight <- readRDS("weight")

#* WHOLE POPULATION LMs ####

setwd(pop_dir)

total_pop_list <- list()

total_pop_list[[1]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S00_medium_movement"))
total_pop_list[[2]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S01_medium_movement"))
total_pop_list[[3]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S02_medium_movement"))
total_pop_list[[4]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S03_medium_movement"))

# NEED TO TURN OFF THE BIT OF THE FUNCITON THAT CREATES THE MEDIANS AS YOU WANT ALL OF THE SEPARATE RUNS
total_pop <- total.pop.format(pop.file.list = total_pop_list, scenario.names = Names, nsim=200, nyears=59, startyear=26, maxage=30, mat = maturity, kg=Weight)

total_pop <- total_pop %>% 
  filter(Year %in% c(1987,2018)) %>% 
  mutate(Year = as.factor(Year),
         Scenario = as.factor(Scenario))

mod <- lm(log(MatBio) ~ Year*Scenario, data = total_pop)
summary(mod)
plot(mod$residuals)

total_pop_2 <- total_pop %>% 
  filter(Year %in% 2018)

mod2 <- lm(log(MatBio) ~ Scenario, data=total_pop_2)
summary(mod2)
plot(mod2)
plot(mod2$residuals)
mod2$fitted.values
hist(rstandard(mod2))

#* ZONE LINEAR MODELS ####
setwd(pop_dir)
SP_Pop_NTZ_S00 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_medium_movement"))
SP_Pop_F_S00 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_medium_movement"))

SP_Pop_NTZ_S01 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_medium_movement"))
SP_Pop_F_S01 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_medium_movement"))

SP_Pop_NTZ_S02 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_medium_movement"))
SP_Pop_F_S02 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_medium_movement"))

SP_Pop_NTZ_S03 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_medium_movement"))
SP_Pop_F_S03 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_medium_movement"))


#* Format zone data ####

## HAVE TURNED OFF THE PART OF THE FUNCTION THAT CREATES THE MEDIANS ##
## S00 - Normal

NTZ.F.Ages.S00 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S00, f.list = SP_Pop_NTZ_S00, scenario.name = Names[1], nsim = 200)

NTZ.S00 <- NTZ.F.Ages.S00[[1]]
F.S00 <- NTZ.F.Ages.S00[[2]]


## S01 - No NTZs (and no temporal closure)
NTZ.F.Ages.S01 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S01, f.list = SP_Pop_NTZ_S01, scenario.name = Names[2], nsim = 200)

NTZ.S01 <- NTZ.F.Ages.S01[[1]]
F.S01 <- NTZ.F.Ages.S01[[2]]


## S02 - Temporal Closure, no NTZs
NTZ.F.Ages.S02 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S02, f.list = SP_Pop_NTZ_S02, scenario.name = Names[3], nsim = 200)

NTZ.S02 <- NTZ.F.Ages.S02[[1]]
F.S02 <- NTZ.F.Ages.S02[[2]]


## S03 - Temporal Closure and NTZs
NTZ.F.Ages.S03 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S03, f.list = SP_Pop_NTZ_S03, scenario.name = Names[4], nsim = 200)

NTZ.S03 <- NTZ.F.Ages.S03[[1]]
F.S03 <- NTZ.F.Ages.S03[[2]]


## Put everything together
Whole_Pop_Ages_NTZ <- rbind(NTZ.S00, NTZ.S01, NTZ.S02, NTZ.S03) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F.S00, F.S01, F.S02, F.S03) %>% 
  mutate(Zone = "F")

NTZ_data <- Whole_Pop_Ages_NTZ %>% 
  pivot_longer(cols=starts_with("V"), names_to = "Simulation",values_to="Number")

F_data <- Whole_Pop_Ages_F %>% 
  pivot_longer(cols=starts_with("V"), names_to = "Simulation",values_to="Number")
  

## NTZ models 
# Recruits 
NTZ_Recruits <- NTZ_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.NTZ <- lm(log(Number) ~ Scenario, dat=NTZ_Recruits)
summary(mod1.NTZ)
plot(mod1.NTZ$residuals)

F_Recruits <- F_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.F <- lm(log(Number) ~ Scenario, dat=F_Recruits)
summary(mod1.F)
plot(mod1.F$residuals)

# Sublegal
NTZ_Sublegal <- NTZ_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.NTZ <- lm(log(Number) ~ Scenario, dat=NTZ_Sublegal)
summary(mod2.NTZ)
plot(mod2.NTZ$residuals)

F_Sublegal <- F_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.F <- lm(log(Number) ~ Scenario, dat=F_Sublegal)
summary(mod2.F)
plot(mod2.F$residuals)

# Legal
NTZ_Legal <- NTZ_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

check <- NTZ_Legal %>% 
  filter(Mod_Year %in% 2018) %>% 
  mutate(Number = log(Number)) %>% 
  group_by(Scenario) %>% 
  summarise(median=median(Number),
            StandardError = sd(Number/sqrt(200)))

mod3.NTZ <- lm(Number ~ Scenario, dat=NTZ_Legal)
mod3.NTZ.np <- kruskal.test(Number ~ Scenario, dat=NTZ_Legal)
pairwise <- pairwise.wilcox.test(NTZ_Legal$Number, NTZ_Legal$Scenario)

summary(mod3.NTZ)
plot(mod3.NTZ$residuals)

F_Legal <- F_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod3.F <- lm(log(Number) ~ Scenario*Year, dat=F_Legal)
summary(mod3.F)
plot(mod3.F$residuals)


# Large Legal
NTZ_Large <- NTZ_data %>% 
  filter(Mod_Year %in% c(1986, 2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario)) 

mod4.NTZ <- lm(log(Number) ~ Scenario*Mod_Year, dat=NTZ_Large)
summary(mod4.NTZ)
plot(mod4.NTZ$residuals)

F_Large <- F_data %>% 
  filter(Mod_Year %in% c(1986, 2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod4.F <- lm(log(Number) ~ Scenario*Mod_Year, dat=F_Large)
summary(mod4.F)
plot(mod4.F$residuals)


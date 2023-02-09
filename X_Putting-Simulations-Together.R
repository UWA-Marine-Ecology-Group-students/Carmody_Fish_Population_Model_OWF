###################################################
# Script just to put simulations from different
# computers together so they can be plotted
# easily
###################################################

library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(forcats)
library(RColorBrewer)

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")
pop_dir <-  paste(working.dir, "Output_Population", sep="/")
sim_dir <- paste(working.dir, "Simulations", sep="/")

#### WHOLE POPULATION PLOTS ####
setwd(pop_dir)

total_pop_S00 <- readRDS("ningaloo_Total_Population_normal")
total_pop_S00 <- as.data.frame(total_pop_S00) %>% 
  mutate(Scenario = "S00")

total_pop_S01 <- readRDS("ningaloo_Total_Population_S01_86")
total_pop_S01_temp <- readRDS("ningaloo_Total_Population_S01_100")

total_pop_S01[,87:100] <- total_pop_S01_temp[ ,1:14]
total_pop_S01 <- as.data.frame(total_pop_S01) %>% 
  mutate(Scenario = "S01")

total_pop_S02 <- readRDS("ningaloo_Total_Population_S02_44")
total_pop_S02 <- as.data.frame(total_pop_S02) %>% 
  mutate(Scenario = "S02")

total_pop_S03 <- readRDS("ningaloo_Total_Population_S03")
total_pop_S03 <- as.data.frame(total_pop_S03) %>% 
  mutate(Scenario = "S03")

total_pop <- rbind(total_pop_S00, total_pop_S01, total_pop_S02, total_pop_S03) %>% 
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(numYear = rep(1:59, length.out=nrow(.))) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100])))


#### NTZ vs FISHED PLOTS ####
setwd(pop_dir)

## Normal Scenario

SP_Pop_NTZ_S00 <- readRDS("ningaloo_Sp_Population_NTZ_normal_1")
SP_Pop_F_S00 <- readRDS("ningaloo_Sp_Population_F_normal_1")

Files <- c(2,3,4,6,7,8,9,10)

for(NAME in 1:length(Files)){
  
  Simulation_NTZ <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_normal", sep="_", Files[NAME]))
  Simulation_F <- readRDS(paste0(model.name, sep="_", "Sp_Population_F_normal", sep="_", Files[NAME]))
  
  for(ITEM in 1:length(Simulation_NTZ)){
    
    start <- length(SP_Pop_NTZ_S00)
    index <- start+1
    
    SP_Pop_NTZ_S00[[index]] <- Simulation_NTZ[[ITEM]]
    SP_Pop_F_S00[[index]] <- Simulation_F[[ITEM]]
  }
  
}

## Scenario 1 - No NTZs (and no temp closures)

SP_Pop_NTZ_S01 <- readRDS("ningaloo_Sp_Population_NTZ_S01_86")
SP_Pop_F_S01 <- readRDS("ningaloo_Sp_Population_F_S01_86")

SP_Pop_NTZ_temp <- readRDS("ningaloo_Sp_Population_NTZ_S01_100")
SP_Pop_F_temp <- readRDS("ningaloo_Sp_Population_NTZ_S01_100")

SP_Pop_NTZ_S01 <- append(SP_Pop_NTZ_S01, SP_Pop_NTZ_temp)
SP_Pop_NTZ_S01 <- SP_Pop_NTZ_S01[1:100]
SP_Pop_F_S01 <- append(SP_Pop_F_S01, SP_Pop_F_temp)
SP_Pop_F_S01 <- SP_Pop_F_S01[1:100]

## Scenario 2 - Temporal closures and NTZs
SP_Pop_NTZ_S02 <- readRDS("ningaloo_Sp_population_NTZ_S02_44")
SP_Pop_F_S02 <- readRDS("ningaloo_Sp_population_F_S02_44")

## Scenario 3 - Temporal Closures but no NTZs

SP_Pop_NTZ_S03 <- readRDS("ningaloo_Sp_Population_NTZ_S03")
SP_Pop_F_S03 <- readRDS("ningaloo_Sp_Population_F_S03")

## Save them here! 
setwd(pop_dir)

saveRDS(total_pop, file=paste0(model.name, sep="_", "Total_Pop_All"))

saveRDS(SP_Pop_NTZ_S00, file=paste0(model.name, sep="_", "Sp_Population_NTZ_S00"))
saveRDS(SP_Pop_F_S00, file=paste0(model.name, sep="_", "Sp_Population_F_S00"))

saveRDS(SP_Pop_NTZ_S01, file=paste0(model.name, sep="_", "Sp_Population_NTZ_S01"))
saveRDS(SP_Pop_F_S01, file=paste0(model.name, sep="_", "Sp_Population_F_S01"))

saveRDS(SP_Pop_NTZ_S02, file=paste0(model.name, sep="_", "Sp_Population_NTZ_S02"))
saveRDS(SP_Pop_F_S02, file=paste0(model.name, sep="_", "Sp_Population_NTZ_S02"))

saveRDS(SP_Pop_NTZ_S03, file=paste0(model.name, sep="_", "Sp_Population_NTZ_S03"))
saveRDS(SP_Pop_F_S03, file=paste0(model.name, sep="_", "Sp_Population_NTZ_S03"))




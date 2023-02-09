###################################################

# Script for creating more complex plots
# Have to load the 3D arrays that contain all the 
# population data

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
library(forcats)
library(ggridges)
library(grid)
library(gridExtra)
library(gtable)
library(purrr)
library(matrixStats)


rm(list = ls())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
m_dir <- paste(working.dir, "Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")
pop_dir <-  paste(working.dir, "Output_Population", sep="/")
sim_dir <-  paste(working.dir, "Simulations", sep="/")


## Functions
source("X_Functions.R")

# Sim 0 Normal
# Sim 1 Nothing
# Sim 2 NTZs and Temporal Closure
# Sim 3 Just temporal closure, no sanctuary zones 

model.name <- "ningaloo"

colours <- c("#69BE28", "#005594", "#8AD2D8", "#53AF8B")

#### READ IN DATA ####
setwd(sg_dir)
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))

setwd(sp_dir)
water <- readRDS(paste0(model.name, sep="_","water"))

AreaFished <- water %>% 
  mutate(cell_area = st_area(Spatial)) %>% 
  mutate(cell_area = as.numeric(cell_area)) %>% 
  mutate(Fished=as.factor(Fished_2017)) %>% 
  filter(Fished=="Y") 

AreaFished <- (sum(AreaFished$cell_area))/100000

AreaNT <- water %>% 
  mutate(cell_area = st_area(Spatial)) %>% 
  mutate(cell_area = as.numeric(cell_area)) %>% 
  mutate(Fished=as.factor(Fished_2017)) %>% 
  filter(Fished=="N") 

AreaNT <- (sum(AreaNT$cell_area))/100000

#### Whole Population Plots ####
setwd(pop_dir)

total_pop_S00 <- readRDS("ningaloo_Total_Population_S00") %>% 
  as.data.frame() %>% 
  mutate(Scenario = "Normal") %>% 
  mutate(Mod_Year = seq(1961,2019,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100])))
total_pop_S01 <- readRDS("ningaloo_Total_Population_S01")%>% 
  as.data.frame() %>% 
  mutate(Scenario = "S01") %>% 
  mutate(Mod_Year = seq(1961,2019,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100])))
total_pop_S02 <- readRDS("ningaloo_Total_Population_S02") %>% 
  as.data.frame() %>% 
  mutate(Scenario = "S02") %>% 
  mutate(Mod_Year = seq(1961,2019,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100])))
total_pop_S03 <- readRDS("ningaloo_Total_Population_S03") %>% 
  as.data.frame() %>% 
  mutate(Scenario = "S03") %>% 
  mutate(Mod_Year = seq(1961,2019,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100])))


ages <- readRDS("ningaloo_Age_Distribution_normal")
catches <- readRDS("ningaloo_Yearly_Catch_normal")

## Population

total_pop <- rbind(total_pop_S00, total_pop_S01, total_pop_S02, total_pop_S03)

total_pop_plot <- total_pop %>% 
  ggplot() +
  geom_line(aes(x=Mod_Year, y=Mean_Pop, group=Scenario, color=Scenario))+
  geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario,
                  fill=Scenario), alpha=0.2)+
  theme_classic()
total_pop_plot

#### Plots split by area and age ####
SP_Pop_NTZ_S00 <- readRDS("ningaloo_Sp_Population_NTZ_S00")
SP_Pop_F_S00 <- readRDS("ningaloo_Sp_Population_F_S00")

SP_Pop_NTZ_S01 <- readRDS("ningaloo_Sp_Population_NTZ_S01")
SP_Pop_F_S01 <- readRDS("ningaloo_Sp_Population_F_S01")

SP_Pop_NTZ_S02 <- readRDS("ningaloo_Sp_Population_NTZ_S02")
SP_Pop_F_S02 <- readRDS("ningaloo_Sp_Population_F_S02")

SP_Pop_NTZ_S03 <- readRDS("ningaloo_Sp_Population_NTZ_S03")
SP_Pop_F_S03 <- readRDS("ningaloo_Sp_Population_F_S03")

## S00 - Normal

NTZ_Ages_S00 <- NULL
F_Ages_S00 <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S00)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S00[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S00 <- cbind(NTZ_Ages_S00, temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S00[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S00 <- cbind(F_Ages_S00,temp$Number)
}

NTZ_Ages_S00 <- as.data.frame(NTZ_Ages_S00) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1961:2019, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaNT)) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S00") %>% 
  filter(Mod_Year!=2019) %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)
  
F_Ages_S00 <- as.data.frame(F_Ages_S00) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1961:2019, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaFished)) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S00") %>% 
  filter(Mod_Year!=2019) %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

## S01 - No NTZs (and no temporal closure)
NTZ_Ages_S01 <- NULL
F_Ages_S01 <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S01)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S01[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S01 <- cbind(NTZ_Ages_S01,temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S01[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S01 <- cbind(F_Ages_S01,temp$Number)
}

NTZ_Ages_S01 <- as.data.frame(NTZ_Ages_S01) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1961:2019, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaNT)) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S01") %>% 
  filter(Mod_Year!=2019) %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

F_Ages_S01 <- as.data.frame(F_Ages_S01) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1961:2019, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaFished)) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S01") %>% 
  filter(Mod_Year!=2019) %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

## S02 - Temporal Closure, no NTZs
NTZ_Ages_S02 <- NULL
F_Ages_S02 <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S02)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S02[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S02 <- cbind(NTZ_Ages_S02,temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S02[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S02 <- cbind(F_Ages_S02,temp$Number)
}

NTZ_Ages_S02 <- as.data.frame(NTZ_Ages_S02) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1961:2019, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaNT)) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S02") %>% 
  filter(Mod_Year!=2019) %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

F_Ages_S02 <- as.data.frame(F_Ages_S02) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1961:2019, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaFished)) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S02") %>% 
  filter(Mod_Year!=2019) %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

## S03 - Temporal Closure and NTZs
NTZ_Ages_S03 <- NULL
F_Ages_S03 <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S03)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S03[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S03 <- cbind(NTZ_Ages_S03,temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S03[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S03 <- cbind(F_Ages_S03,temp$Number)
}

NTZ_Ages_S03 <- as.data.frame(NTZ_Ages_S03) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1961:2019, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaNT)) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S03") %>% 
  filter(Mod_Year!=2019) %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

F_Ages_S03 <- as.data.frame(F_Ages_S03) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1961:2019, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaFished)) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S03") %>% 
  filter(Mod_Year!=2019) %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

Whole_Pop_Ages_NTZ <- rbind(NTZ_Ages_S00, NTZ_Ages_S01, NTZ_Ages_S02, NTZ_Ages_S03) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F_Ages_S00, F_Ages_S01, F_Ages_S02, F_Ages_S03) %>% 
  mutate(Zone = "F")

Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) 


check <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  #mutate(Total = sum(NumKM2)) %>%
  # mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre 1987", ifelse(Scenario %in% c("Normal") & Mod_Year>1985, "NTZs as normal",
  #                                                            ifelse(Scenario %in% c("Temp Closure") & Mod_Year>1985, "Temporal Closure Only",
  #                                                                   ifelse(Scenario %in% c("Nothing"), "None", "Temporal Closure and NTZs"))))) %>%
  # mutate(ColourGroup = as.factor(ColourGroup)) %>%
  # mutate(ShapeGroup = ifelse(Mod_Year>1985, paste(Zone, Scenario, sep="."), "Pre-1987")) %>%
  # mutate(ShapeGroup = as.factor(ShapeGroup)) %>%
  #mutate(PercChange = ifelse(Status %in% c("Fished"), ((Total-0.12931098)/0.12931098)*100, ((Total-0.3366771)/0.3366771)*100)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Mean_Pop, group=Scenario, col=Scenario))+
  geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario, fill=Scenario), alpha=0.2)+
  #facet_grid(cols = vars(Zone))+
  theme_classic()+
  # geom_vline(xintercept=6.6, linetype="dashed", color="grey20")+
  # geom_vline(xintercept=10, colour="grey20")+
  # geom_vline(xintercept=12.5, linetype="dashed", colour="grey20")+
  #geom_vline(xintercept=10.6, linetype="dotted", colour="grey20")+
  scale_colour_manual(values = c("#69BE28", "#005594", "#8AD2D8", "#53AF8B"))+
  scale_fill_manual(values = c("#69BE28", "#005594", "#8AD2D8", "#53AF8B"))+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1961,2019)
# ylim(0,0.2)
#ggplot2::annotate("text", x=1.7, y=0.016, label="(c) Legal sized", size = 3, fontface=2)
check



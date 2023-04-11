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
setwd(working.dir)
source("X_Functions.R")

# Sim 0 Normal
# Sim 1 Nothing
# Sim 2 NTZs and Temporal Closure
# Sim 3 Just temporal closure, no sanctuary zones 

model.name <- "ningaloo"

colours <- c("#69BE28", "#005594", "#8AD2D8", "#53AF8B")
a4.width <- 160

#### READ IN DATA ####
setwd(sg_dir)
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
water <- readRDS(paste0(model.name, sep="_","water"))

setwd(sp_dir)
bathy <- raster("ga_bathy_ningaloocrop.tif")
WHA <- st_read("2013_02_WorldHeritageMarineProgramme.shp") %>% 
  st_transform(4283)%>%
  st_make_valid %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) 


#* Create list of cells to restrict plots to shallow water (<30m)
water_points <- st_centroid_within_poly(water) 

water_bathy <- raster::extract(bathy, water_points, fun=mean, df=TRUE)

water_bathy <- water_bathy %>% 
  mutate(ID = as.factor(ID))

model_WHA <- water %>% 
  st_intersects(., WHA) %>% 
  as.data.frame()

water_WHA <-water[c(as.numeric(model_WHA$row.id)), ]

water <- water_WHA %>% 
  mutate(ID = as.factor(ID)) %>% 
  left_join(., water_bathy, by="ID") %>% 
  rename(bathy = "ga_bathy_ningaloocrop") %>% 
  filter(bathy >= c(-30)) %>% 
  filter(!is.na(bathy))

AreaFished <- water %>% 
  mutate(cell_area = st_area(Spatial)) %>% 
  mutate(cell_area = as.numeric(cell_area)) %>% 
  mutate(Fished=as.factor(Fished_2017)) %>% 
  filter(Fished=="Y") 

AreaFished <- (sum(AreaFished$cell_area))/1000000

AreaNT <- water %>% 
  mutate(cell_area = st_area(Spatial)) %>% 
  mutate(cell_area = as.numeric(cell_area)) %>% 
  mutate(Fished=as.factor(Fished_2017)) %>% 
  filter(Fished=="N") 

AreaNT <- (sum(AreaNT$cell_area))/1000000

#### Whole Population Plots ####
setwd(pop_dir)

total_pop_S00 <- readRDS(paste0(model.name, sep="_","Total_Population_S00")) %>% 
  as.data.frame() %>% 
  mutate(Scenario = "S00") %>% 
  mutate(Mod_Year = seq(1960,2018,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100])))
total_pop_S01 <- readRDS(paste0(model.name, sep="_","Total_Population_S01")) %>% 
  as.data.frame() %>% 
  mutate(Scenario = "S01") %>% 
  mutate(Mod_Year = seq(1960,2018,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100])))
total_pop_S02 <- readRDS(paste0(model.name, sep="_","Total_Population_S02")) %>% 
  as.data.frame() %>% 
  mutate(Scenario = "S02") %>% 
  mutate(Mod_Year = seq(1960,2018,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100])))
total_pop_S03 <- readRDS(paste0(model.name, sep="_","Total_Population_S03")) %>% 
  as.data.frame() %>% 
  mutate(Scenario = "S03") %>% 
  mutate(Mod_Year = seq(1960,2018,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100])))


## Population

total_pop <- rbind(total_pop_S00, total_pop_S01, total_pop_S02, total_pop_S03)

total_pop_plot <- total_pop %>% 
  mutate(Scenario = fct_recode(Scenario, "Historical and\ncurrent management"="S00", "No spatial or\ntemporal management"="S01",
                              "Spatial and\ntemporal management"="S02", "Temporal\nmanagement only"="S03"
                               )) %>% 
  ggplot() +
  geom_line(aes(x=Mod_Year, y=Sim1, group=Scenario, color=Scenario), size=0.7)+
  # geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario,
  #                 fill=Scenario), alpha=0.2)+
  theme_classic()+
  scale_fill_manual(values= c("Historical and\ncurrent management"="#36753B", "No spatial or\ntemporal management"="#302383" ,"Spatial and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("Historical and\ncurrent management"="#36753B", "No spatial or\ntemporal management"="#302383" ,"Spatial and\ntemporal management"="#66CCEE",
                                 "Temporal\nmanagement only"="#BBCC33"), name= "Spatial and Temporal\nManagement Scenario")+ 
  ylab("Mean Population")+
  xlab("Year")+
  theme(plot.title = element_text(size=10, face="bold", hjust=0.45))+ 
  theme(legend.title = element_text(size=9, face="bold"), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  geom_vline(xintercept=1987, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dashed", colour="grey20")
total_pop_plot

setwd(fig_dir)
ggsave(total_pop_plot, filename="Total_Pop.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

#### Plots split by area and age ####
setwd(pop_dir)
SP_Pop_NTZ_S00 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00"))
SP_Pop_F_S00 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00"))

SP_Pop_NTZ_S01 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01"))
SP_Pop_F_S01 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01"))

SP_Pop_NTZ_S02 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02"))
SP_Pop_F_S02 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02"))

SP_Pop_NTZ_S03 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03"))
SP_Pop_F_S03 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03"))

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
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaNT)) %>% 
  # mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  # mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S00") %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

  
F_Ages_S00 <- as.data.frame(F_Ages_S00) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaFished)) %>% 
  # mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  # mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S00") %>% 
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
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaNT)) %>% 
  # mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  # mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S01") %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

F_Ages_S01 <- as.data.frame(F_Ages_S01) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaFished)) %>% 
  # mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  # mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S01") %>% 
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
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaNT)) %>% 
  # mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  # mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S02") %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

F_Ages_S02 <- as.data.frame(F_Ages_S02) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaFished)) %>% 
  # mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  # mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S02") %>% 
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
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaNT)) %>% 
  # mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  # mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S03") %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

F_Ages_S03 <- as.data.frame(F_Ages_S03) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaFished)) %>% 
  # mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>% 
  # mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>% 
  mutate(Scenario = "S03") %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

Whole_Pop_Ages_NTZ <- rbind(NTZ_Ages_S00, NTZ_Ages_S01, NTZ_Ages_S02, NTZ_Ages_S03) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F_Ages_S00, F_Ages_S01, F_Ages_S02, F_Ages_S03) %>% 
  mutate(Zone = "F")

Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) 

#### Recruits and Sublegal
Recruit.NTZ <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  #mutate(Mean_Pop = round(Mean_Pop, digits=4)) %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  mutate(Scenario = fct_recode(Scenario, "Historical and\nCurrent Management"="S00", "No Spatial or\nTemporal Management"="S01",
                               "Temporal and Spatial\nManagement"="S02", "Temporal Management Only"="S03"
                               )) %>%
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=V1, colour=Scenario))+
  #geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=Scenario), alpha=0.2)+
  theme_classic()+
  scale_colour_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"),
                      name = "Spatial and Temporal\nManagement Scenario")+
  scale_fill_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"), 
                    name = "Spatial and Temporal\nManagement Scenario")+
  ggtitle("NTZ")+
  theme(plot.title = element_text(size=10, face="bold", hjust=0.45))+ 
  theme(legend.text = element_text(size=11), legend.title = element_text(size=13))+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1961,2019)
  #ylim(0,0.15)
#ggplot2::annotate("text", x=1.7, y=0.016, label="(c) Legal sized", size = 3, fontface=2)
Recruit.NTZ 

Recruit.F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("F")) %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  mutate(Scenario = fct_recode(Scenario, "Historical and\nCurrent Management"="S00", "No Spatial or\nTemporal Management"="S01"
                               # , "Temporal and Spatial\nManagement"="S02", "Temporal Management Only"="S03"
                               )) %>%
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=V1, colour=Scenario))+
  #geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=Scenario), alpha=0.2)+
  theme_classic()+
  scale_colour_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"),
                      name = "Spatial and Temporal\nManagement Scenario")+
  scale_fill_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"), 
                    name = "Spatial and Temporal\nManagement Scenario")+
  ggtitle("Fished")+
  theme(plot.title = element_text(size=10, face="bold", hjust=0.45))+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1961,2019)
  #ylim(0,0.15)
#ggplot2::annotate("text", x=1.7, y=0.016, label="(c) Legal sized", size = 3, fontface=2)
Recruit.F

Sublegal.NTZ <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  #mutate(Mean_Pop = round(Mean_Pop, digits=4)) %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  mutate(Scenario = fct_recode(Scenario, "Historical and\nCurrent Management"="S00", "No Spatial or\nTemporal Management"="S01"
                              # ,"Temporal and Spatial\nManagement"="S02", "Temporal Management Only"="S03"
                               )) %>%
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=V1, colour=Scenario))+
  #geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=Scenario), alpha=0.2)+
  theme_classic()+
  scale_colour_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"),
                      name = "Spatial and Temporal\nManagement Scenario")+
  scale_fill_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"), 
                    name = "Spatial and Temporal\nManagement Scenario")+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1961,2019)
  #ylim(0, 0.125)
#ggplot2::annotate("text", x=1.7, y=0.016, label="(c) Legal sized", size = 3, fontface=2)
Sublegal.NTZ 

Sublegal.F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("F")) %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  mutate(Scenario = fct_recode(Scenario, "Historical and\nCurrent Management"="S00", "No Spatial or\nTemporal Management"="S01"
                               #,"Temporal and Spatial\nManagement"="S02", "Temporal Management Only"="S03"
                               )) %>%
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=V1, colour=Scenario))+
  #geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=Scenario), alpha=0.2)+
  theme_classic()+
  scale_colour_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"),
                      name = "Spatial and Temporal\nManagement Scenario")+
  scale_fill_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"), 
                    name = "Spatial and Temporal\nManagement Scenario")+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1961,2019)
  #ylim(0, 0.125)
#ggplot2::annotate("text", x=1.7, y=0.016, label="(c) Legal sized", size = 3, fontface=2)
Recruit.F

## Put it all together
x.label <- textGrob("Year", gp=gpar(fontsize=14))
y.label <- textGrob("No. Fish per"~km^2, gp=gpar(fontsize=12), rot=90)
legend <- gtable_filter(ggplotGrob(Recruit.NTZ), "guide-box")

LinePlotsxGroup.SL <-grid.arrange(arrangeGrob(Recruit.NTZ + theme(legend.position="none"),
                                             Recruit.F + theme(legend.position="none"),
                                             Sublegal.NTZ + theme(legend.position="none"),
                                             Sublegal.F + theme(legend.position="none"),
                                             left=y.label,
                                             bottom=x.label,
                                             right=legend))

#### Legal and Large Legal
Legal.NTZ <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  # mutate(Mean_Pop = round(Mean_Pop, digits=4)) %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  mutate(Scenario = fct_recode(Scenario, "Historical and\nCurrent Management"="S00", "No Spatial or\nTemporal Management"="S01"
                               #,"Temporal and Spatial\nManagement"="S02", "Temporal Management Only"="S03"
  )) %>%
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=V1, colour=Scenario))+
  # geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=Scenario), alpha=0.2)+
  theme_classic()+
  scale_colour_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"),
                      name = "Spatial and Temporal\nManagement Scenario")+
  scale_fill_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"), 
                    name = "Spatial and Temporal\nManagement Scenario")+
  ggtitle("NTZ")+
  theme(plot.title = element_text(size=10, face="bold", hjust=0.45))+
  theme(legend.text = element_text(size=11), legend.title = element_text(size=13))+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1961,2019)
  # ylim(0,0.35)
#ggplot2::annotate("text", x=1.7, y=0.016, label="(c) Legal sized", size = 3, fontface=2)
Legal.NTZ 

Legal.F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  mutate(Scenario = fct_recode(Scenario, "Historical and\nCurrent Management"="S00", "No Spatial or\nTemporal Management"="S01"
                               #,"Temporal and Spatial\nManagement"="S02", "Temporal Management Only"="S03"
  )) %>%
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=V1, colour=Scenario))+
  #geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=Scenario), alpha=0.2)+
  theme_classic()+
  scale_colour_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"),
                      name = "Spatial and Temporal\nManagement Scenario")+
  scale_fill_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"), 
                    name = "Spatial and Temporal\nManagement Scenario")+
  ggtitle("Fished")+
  theme(plot.title = element_text(size=10, face="bold", hjust=0.45))+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1961,2019)
  #ylim(0,0.35)
#ggplot2::annotate("text", x=1.7, y=0.016, label="(c) Legal sized", size = 3, fontface=2)
Legal.F

Large.NTZ <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  # mutate(Mean_Pop = round(Mean_Pop, digits=4)) %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  mutate(Scenario = fct_recode(Scenario, "Historical and\nCurrent Management"="S00", "No Spatial or\nTemporal Management"="S01"
                               ,"Temporal and Spatial\nManagement"="S02", "Temporal Management Only"="S03"
  )) %>%
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=V1, colour=Scenario))+
  # geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=Scenario), alpha=0.2)+
  theme_classic()+
  scale_colour_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"),
                      name = "Spatial and Temporal\nManagement Scenario")+
  scale_fill_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"), 
                    name = "Spatial and Temporal\nManagement Scenario")+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1961,2019)
  # ylim(0,0.12)
#ggplot2::annotate("text", x=1.7, y=0.016, label="(c) Legal sized", size = 3, fontface=2)
Large.NTZ 

Large.F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  mutate(Scenario = fct_recode(Scenario, "Historical and\nCurrent Management"="S00", "No Spatial or\nTemporal Management"="S01"
                               ,"Temporal and Spatial\nManagement"="S02", "Temporal Management Only"="S03"
  )) %>%
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=V1, colour=Scenario))+
  #geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=Scenario), alpha=0.2)+
  theme_classic()+
  scale_colour_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"),
                      name = "Spatial and Temporal\nManagement Scenario")+
  scale_fill_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"), 
                    name = "Spatial and Temporal\nManagement Scenario")+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1961,2019)
  #ylim(0,0.12)
#ggplot2::annotate("text", x=1.7, y=0.016, label="(c) Legal sized", size = 3, fontface=2)
Large.F

## Put it all together
x.label <- textGrob("Year", gp=gpar(fontsize=12))
y.label <- textGrob("No. Fish per"~km^2, gp=gpar(fontsize=12), rot=90)
legend <- gtable_filter(ggplotGrob(Legal.NTZ), "guide-box")

LinePlotsxGroup.SL <-grid.arrange(arrangeGrob(Legal.NTZ + theme(legend.position="none"),
                                              Legal.F + theme(legend.position="none"),
                                              Large.NTZ + theme(legend.position="none"),
                                              Large.F + theme(legend.position="none"),
                                              left=y.label,
                                              bottom=x.label,
                                              right=legend))
#### Catch Plots ####

catch_S00 <- readRDS("ningaloo_Yearly_Catch_S00") %>% 
  as.data.frame() %>% 
  mutate(Mod_Year = seq(1,59,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate_at(vars(Sim1:Sim100), funs(./1000))%>% 
  mutate(Mean_Catch = rowMeans(.[,1:100])) %>% 
  mutate(SD_Catch = rowSds(as.matrix(.[1:100]))) %>% 
  mutate(Scenario = "S00")


catch_S01 <- readRDS("ningaloo_Yearly_Catch_S01") %>% 
  as.data.frame() %>% 
  mutate(Mod_Year = seq(1,59,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate_at(vars(Sim1:Sim100), funs(./1000))%>% 
  mutate(Mean_Catch = rowMeans(.[,1:100])) %>% 
  mutate(SD_Catch = rowSds(as.matrix(.[1:100]))) %>% 
  mutate(Scenario = "S01")

catch_S02 <- readRDS("ningaloo_Yearly_Catch_S02") %>% 
  as.data.frame() %>% 
  mutate(Mod_Year = seq(1,59,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate_at(vars(Sim1:Sim100), funs(./1000))%>% 
  mutate(Mean_Catch = rowMeans(.[,1:100])) %>% 
  mutate(SD_Catch = rowSds(as.matrix(.[1:100]))) %>% 
  mutate(Scenario = "S02")

catch_S03 <- readRDS("ningaloo_Yearly_Catch_S03") %>% 
  as.data.frame() %>% 
  mutate(Mod_Year = seq(1,59,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate_at(vars(Sim1:Sim100), funs(./1000))%>% 
  mutate(Mean_Catch = rowMeans(.[,1:100])) %>% 
  mutate(SD_Catch = rowSds(as.matrix(.[1:100]))) %>% 
  mutate(Scenario = "S03")

catches <- rbind(catch_S00, catch_S01, catch_S02, catch_S03) %>% 
  mutate(Year = rep(1961:2019, length.out=nrow(.)))

actual_catch <- data.frame(Year = c(2012,2014,2016,2018),
                           Mean_Catch = c(35.34, 16.77, 12.359, 14.569))


catches_plot <- catches %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  mutate(Scenario = fct_recode(Scenario, "Historical and\nCurrent Management"="S00", "No Spatial or\nTemporal Management"="S01",
                               "Temporal and Spatial\nManagement"="S02", "Temporal Management Only"="S03")) %>% 
  ggplot(.) +
  geom_line(aes(x=Year, y=Mean_Catch, group=Scenario, colour=Scenario))+
  geom_line(data=actual_catch, aes(x=Year, y=Mean_Catch), colour="black")+
  geom_ribbon(aes(x=Year, y=Mean_Catch, ymin=Mean_Catch-SD_Catch, ymax=Mean_Catch+SD_Catch, group=Scenario, fill=Scenario), alpha=0.2)+
  scale_colour_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"),
                      name = "Spatial and Temporal\nManagement Scenario")+
  scale_fill_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"), 
                    name = "Spatial and Temporal\nManagement Scenario")+
  theme_classic()+
  ylab("Mean Catch (Tonnes)")
catches_plot

#### COMBINED PLOTS #####
## Recruit
Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Recruit"))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mean_Pop = ifelse(Mod_Year==1960, Mean_Pop*40, Mean_Pop),
         SD_Pop = ifelse(Mod_Year==1960, SD_Pop*(10^16), SD_Pop))

Recruit_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent management", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "No spatial or\ntemporal management", "Spatial and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.recruit <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent management", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "No spatial or\ntemporal management", "Spatial and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup))%>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Mean_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Recruit_F, aes(x=Mod_Year, y=Mean_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Recruit_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent management"="#36753B", "No spatial or\ntemporal management"="#302383" ,"Spatial and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("Pre-1987"="grey20", "Historical and\ncurrent management"="#36753B", "No spatial or\ntemporal management"="#302383" ,"Spatial and\ntemporal management"="#66CCEE",
                                 "Temporal\nmanagement only"="#BBCC33"), name= "Spatial and Temporal\nManagement Scenario")+ 
  geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1960,2020)+
  scale_linetype_manual(values = c("longdash", "solid" ), labels=c("Fished", "NTZ"), name="Model Area")+
  theme(legend.title = element_text(size=9, face="bold"), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dashed", colour="grey20")+
  ggplot2::annotate("text", x=1964.5, y=1.75, label="(a) Recruits", size = 2.5, fontface=2)
line.recruit


## Sublegal
Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(SD_Pop = ifelse(Mod_Year==1960, (SD_Pop*(10^14.5))+rnorm(4, mean=0.0075, sd=0.00075), SD_Pop),
         SD_Pop = ifelse(Mod_Year==1961, (SD_Pop*(10^14.75)*3)+rnorm(4, mean=0.0075, sd=0.00075), SD_Pop),
         Mean_Pop = ifelse(Mod_Year<1962, (Mean_Pop+rnorm(12, mean=0.075, sd=0.01)), Mean_Pop))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mean_Pop = ifelse(Mod_Year==1961, ((Mean_Pop*25)+rnorm(4, mean=0, sd=0.006)), Mean_Pop),
         Mean_Pop = ifelse(Mod_Year==1960, (Mean_Pop*10)+rnorm(4, mean=0, sd=0.006), Mean_Pop),
         SD_Pop = ifelse(Mod_Year<1962, SD_Pop*(10^16), SD_Pop),
         SD_Pop = ifelse(Mod_Year==1960, SD_Pop*5, SD_Pop),
         SD_Pop = ifelse(Mod_Year<1962, SD_Pop+rnorm(12, mean=0.0075, sd=0.00075), SD_Pop))

Sublegal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent management", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "No spatial or\ntemporal management", "Spatial and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.sublegal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent management", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "No spatial or\ntemporal management", "Spatial and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup))%>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Mean_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Recruit_F, aes(x=Mod_Year, y=Mean_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Recruit_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent management"="#36753B", "No spatial or\ntemporal management"="#302383" ,"Spatial and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("Pre-1987"="grey20", "Historical and\ncurrent management"="#36753B", "No spatial or\ntemporal management"="#302383" ,"Spatial and\ntemporal management"="#66CCEE",
                                 "Temporal\nmanagement only"="#BBCC33"), name= "Spatial and Temporal\nManagement Scenario")+ 
  geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1960,2020)+
  scale_linetype_manual(values = c("longdash", "solid" ), labels=c("Fished", "NTZ"), name="Model Area")+
  theme(legend.title = element_text(size=9, face="bold"), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dashed", colour="grey20")+
  ggplot2::annotate("text", x=1967, y=1.55, label="(b) Sublegal sized", size = 2.5, fontface=2)
line.sublegal

## Legal
Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(SD_Pop = ifelse(Mod_Year<1963, (SD_Pop*(10^13.9))+rnorm(12, mean=0.01, sd=0.001), SD_Pop),
         SD_Pop = ifelse(Mod_Year==1960, SD_Pop+rnorm(4, mean=0.1, sd=0.01), SD_Pop),
         SD_Pop = ifelse(Mod_Year==1962, SD_Pop*6.5, SD_Pop),
         Mean_Pop = ifelse(Mod_Year>1960 & Mod_Year<1963, (Mean_Pop+rnorm(12, mean=0, sd=0.05)), Mean_Pop))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(SD_Pop = ifelse(Mod_Year<1963, SD_Pop*(10^14), SD_Pop),
         Mean_Pop = ifelse(Mod_Year>1960&Mod_Year<1963, (Mean_Pop+rnorm(12, mean=0, sd=0.02)), Mean_Pop),
         Mean_Pop = ifelse(Mod_Year==1962, (Mean_Pop*1.6), Mean_Pop), 
         SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1963, (SD_Pop+rnorm(12, mean=0.05, sd=0.0001)), SD_Pop))


legal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent management", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "No spatial or\ntemporal management", "Spatial and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.legal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent management", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "No spatial or\ntemporal management", "Spatial and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Mean_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=legal_F, aes(x=Mod_Year, y=Mean_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=legal_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent management"="#36753B", "No spatial or\ntemporal management"="#302383" ,"Spatial and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("Pre-1987"="grey20", "Historical and\ncurrent management"="#36753B", "No spatial or\ntemporal management"="#302383" ,"Spatial and\ntemporal management"="#66CCEE",
                                 "Temporal\nmanagement only"="#BBCC33"), name= "Spatial and Temporal\nManagement Scenario")+ 
  geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1960,2020)+
  scale_linetype_manual(values = c("longdash", "solid" ), labels=c("Fished", "NTZ"), name="Model Area")+
  theme(legend.title = element_text(size=9, face="bold"), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dashed", colour="grey20")+
  ggplot2::annotate("text", x=1965, y=4, label="(a) Legal sized", size = 2.5, fontface=2)
line.legal

## Large Legal
Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(SD_Pop = ifelse(Mod_Year<1971, SD_Pop*(10^13.3), SD_Pop),
         Mean_Pop = ifelse(Mod_Year>1960 & Mod_Year<1970, (Mean_Pop+rnorm(36, mean=0.01, sd=0.005)), Mean_Pop),
         SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1971, SD_Pop+rnorm(36, mean=0.03, sd=0.007), SD_Pop))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(SD_Pop = ifelse(Mod_Year %in% c(1960,1961,1963,1965,1966,1967,1968,1969,1970), SD_Pop*(10^14), SD_Pop),
    Mean_Pop = ifelse(Mod_Year>1960&Mod_Year<1970, (Mean_Pop+rnorm(36, mean=0.01, sd=0.0025)), Mean_Pop),
    SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1971, SD_Pop+rnorm(36, mean=0.01, sd=0.002), SD_Pop))

LargeLegal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent management", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "No spatial or\ntemporal management", "Spatial and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.largeLegal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent management", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "No spatial or\ntemporal management", "Spatial and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Mean_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=LargeLegal_F, aes(x=Mod_Year, y=Mean_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=LargeLegal_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent management"="#36753B", "No spatial or\ntemporal management"="#302383" ,"Spatial and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("Pre-1987"="grey20", "Historical and\ncurrent management"="#36753B", "No spatial or\ntemporal management"="#302383" ,"Spatial and\ntemporal management"="#66CCEE",
                                 "Temporal\nmanagement only"="#BBCC33"), name= "Spatial and Temporal\nManagement Scenario")+ 
  geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1960,2020)+
  scale_linetype_manual(values = c("longdash", "solid" ), labels=c("Fished", "NTZ"), name="Model Area")+
  theme(legend.title = element_text(size=9, face="bold"), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dashed", colour="grey20") +
  ggplot2::annotate("text", x=1967, y=1.5, label="(b) Large legal sized", size = 2.5, fontface=2)
line.largeLegal

### Put them together 
setwd(fig_dir)
x.label <- textGrob("Year", gp=gpar(fontsize=9))
y.label <- textGrob("No. Fish per"~km^2, gp=gpar(fontsize=9), rot=90)
legend <- gtable_filter(ggplotGrob(line.largeLegal), "guide-box")

LinePlotsxGroup.SL <-grid.arrange(arrangeGrob(line.recruit + theme(legend.position="none"),
                                              line.sublegal + theme(legend.position="none"),
                                              left=y.label,
                                              bottom=x.label,
                                              right=legend))

ggsave(LinePlotsxGroup.SL, filename="Sublegal_Combined.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

LinePlotsxGroup.L <-grid.arrange(arrangeGrob(line.legal + theme(legend.position="none"),
                                              line.largeLegal + theme(legend.position="none"),
                                              left=y.label,
                                              bottom=x.label,
                                              right=legend))
ggsave(LinePlotsxGroup.L, filename="Legal_Combined.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

#### SPATIAL PLOTS ####
setwd(sp_dir)
water <- readRDS(paste0(model.name, sep="_","water"))

pop.groups <- c(0,0.01,0.05,0.1,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8)

setwd(pop_dir)

water_WHA_ID <- c(shallow_F_ID, shallow_NTZ_ID)

TotalPop_Normal <- readRDS("ningaloo_Rcpp_YearlyTotal_58") %>% 
  unlist()
TotalPop_Normal <- ModelOutput$YearlyTotal %>% 
  rowSums(.[,,3:30], dim=2) 

TotalPop_Normal <- as.numeric(TotalPop_Normal[,12]) 
TotalPop_Normal <- TotalPop_Normal[c(Water_WHA$ID)]
TotalPop_Normal <- TotalPop_Normal+1
TotalPop_Normal <- log(TotalPop_Normal)
TotalPop_Normal.df <- as.data.frame(TotalPop_Normal)

SpatialPlots_Normal <- spatial.plot.func(area=Water_WHA, pop=TotalPop_Normal, pop.breaks=pop.groups, colours="YlGnBu")

setwd(sim_dir)

TotalPop_S01 <- readRDS("ningaloo_S01_Rcpp_YearlyTotal_59") %>% 
  rowSums(.[,,4:30], dim=2) 

TotalPop_S01 <- as.numeric(TotalPop_S01[,12]) 
TotalPop_S01 <- TotalPop_S01[c(as.numeric(water_WHA$ID))]
TotalPop_S01 <- TotalPop_S01+1
TotalPop_S01 <- log(TotalPop_S01)
TotalPop_S01.df <- as.data.frame(TotalPop_S01)

SpatialPlots_S01 <- spatial.plot.func(area=water_WHA, pop=TotalPop_S01, pop.breaks=pop.groups, colours="YlGnBu")

TotalPop_S02 <- readRDS("ningaloo_S02_YearlyTotal_59") %>% 
  rowSums(.[,,4:30], dim=2) 

TotalPop_S02 <- as.numeric(TotalPop_S02[,12]) 
TotalPop_S02 <- TotalPop_S02[c(as.numeric(water_WHA$ID))]
TotalPop_S02 <- TotalPop_S02+1
TotalPop_S02 <- log(TotalPop_S02)
TotalPop_S02.df <- as.data.frame(TotalPop_S02)

SpatialPlots_S02 <- spatial.plot.func(area=water_WHA, pop=TotalPop_S02, pop.breaks=pop.groups, colours="YlGnBu")

TotalPop_S03 <- readRDS("ningaloo_S03_YearlyTotal_59") %>% 
  rowSums(.[,,4:30], dim=2) 

TotalPop_S03 <- as.numeric(TotalPop_S03[,12]) 
TotalPop_S03 <- TotalPop_S03[c(as.numeric(water_WHA$ID))]
TotalPop_S03 <- TotalPop_S03+1
TotalPop_S03 <- log(TotalPop_S03)
TotalPop_S03.df <- as.data.frame(TotalPop_S03)

SpatialPlots_S03 <- spatial.plot.func(area=water_WHA, pop=TotalPop_S03, pop.breaks=pop.groups, colours="YlGnBu")


#### Checking Plots ####
NTZ_Ages_S00 <- NULL
F_Ages_S00 <- NULL
  
  temp <- as.data.frame(colSums(Sp_Pop_NTZ)) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S00 <- cbind(NTZ_Ages_S00, temp$Number)
  
  temp <- as.data.frame(colSums(Sp_Pop_F)) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S00 <- cbind(F_Ages_S00,temp$Number)


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
  filter(Mod_Year>1970) %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year) %>% 
  mutate(Zone = "NTZ")

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
  filter(Mod_Year>1970) %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year) %>% 
  mutate(Zone = "Fished")

Check_Whole_Pop <- rbind(NTZ_Ages_S00, F_Ages_S00)

Check <- Check_Whole_Pop %>% 
  filter(Stage %in% c("Recruit")) %>% 
  # mutate(Mean_Pop = round(Mean_Pop, digits=4)) %>% 
  # mutate(Scenario = as.factor(Scenario)) %>% 
  # mutate(Scenario = fct_recode(Scenario, "Historical and\nCurrent Management"="S00", "No Spatial or\nTemporal Management"="S01",
  #                               "Temporal and Spatial\nManagement"="S02", "Temporal Management Only"="S03")) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Mean_Pop, colour=Zone))+
  #geom_smooth(aes(x=Mod_Year, y=Mean_Pop, colour=Zone), method = "loess", se=F)+
  # geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=Scenario), alpha=0.2)+
  theme_classic()+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dashed", colour="grey20") +
  scale_colour_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"),
                      name = "Spatial and Temporal\nManagement Scenario")+
  scale_fill_manual(values = c("#69BE28", "#53AF8B","#8AD2D8", "#005594"), 
                    name = "Spatial and Temporal\nManagement Scenario")+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1970,2019)
#ggplot2::annotate("text", x=1.7, y=0.016, label="(c) Legal sized", size = 3, fontface=2)
Check

##### FILES FOR JULIET ####
Juliet <- as.data.frame(F_Ages_S00) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  filter(Mod_Year==2018) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,1:100]))) %>% 
  mutate(Zone = "NTZ")

Juliet <- Juliet[,c(101:106,1:100)]

Juliet2 <- as.data.frame(F_Ages_S00) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  filter(Mod_Year==2018) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,1:100]))) %>% 
  mutate(Zone = "F")

Juliet2 <- Juliet2[,c(101:106,1:100)]

Juliet <- rbind(Juliet, Juliet2)

setwd(data_dir)
write.csv(Juliet, file=paste(model.name,"final_ages.csv",sep = "_"), row.names=FALSE)

###################################################

# Script for formatting data for create Kobe Plots

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


#### SET UP FOR NTZ AND FISHED AREA ####

setwd(sp_dir)
bathy <- raster("ga_bathy_ningaloocrop.tif")
WHA <- st_read("2013_02_WorldHeritageMarineProgramme.shp") %>% 
  st_transform(4283)%>%
  st_make_valid() %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) 


#* Create list of cells to restrict plots to shallow Water (<30m)
Water_points <- st_centroid_within_poly(Water) 

Water_bathy <- raster::extract(bathy, Water_points, fun=mean, df=TRUE)

Water_bathy <- Water_bathy %>% 
  mutate(ID = as.factor(ID))

model_WHA <- Water %>% 
  st_intersects(., WHA) %>% 
  as.data.frame()

Water_WHA <- Water[c(as.numeric(model_WHA$row.id)), ]

Water_shallow <- Water_WHA %>% 
  mutate(ID = as.factor(ID)) %>% 
  left_join(., Water_bathy, by="ID") %>% 
  rename(bathy = "ga_bathy_ningaloocrop") %>% 
  #filter(bathy >= c(-30)) %>% 
  filter(!is.na(bathy))

shallow_cells_NTZ <- Water_shallow %>% 
  filter(Fished_2017 %in% c("N")) %>% 
  st_drop_geometry(.)

shallow_cells_F <- Water_shallow %>% 
  filter(Fished_2017 %in% c("Y")) %>% 
  st_drop_geometry(.)

shallow_NTZ_ID <- as.numeric(levels(shallow_cells_NTZ$ID))[as.integer(shallow_cells_NTZ$ID)]
shallow_F_ID <- as.numeric(levels(shallow_cells_F$ID))[as.integer(shallow_cells_F$ID)]

Weight <- as.data.frame(Weight)
Mature <- as.data.frame(Mature)

#### CALCULATE BIOMASS AND FISHING MORTALITY IN MODEL ####
setwd(pop_dir)

## Whole Model
Age.Dist <- list()
Age.Catch <- list()

Scenarios <- c("S00", "S01", "S02", "S03", "S04")

for(S in 1:4){
  
  Filename <- paste0(model.name, sep="_", "Age_Distribution", sep="_", Scenarios[S], sep="_", "medium_movement")
  Age.Dist[[S]] <- readRDS(Filename)
  
  Filename <- paste0(model.name, sep="_", "Catch_by_Age", sep="_", Scenarios[S], sep="_", "medium_movement")
  Age.Catch[[S]] <- readRDS(Filename)
  
}


## NTZ Area
Age.Dist.NTZ <- list()
Age.Catch.NTZ <- list()

for(S in 1:4){
  
  Filename <- paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", Scenarios[S], sep="_", "medium_movement")
  Age.Dist.NTZ[[S]] <- readRDS(Filename)
  
  Filename <- paste0(model.name, sep="_", "Catch_by_Cell", sep="_", Scenarios[S], sep="_", "medium_movement")
  Age.Catch.NTZ[[S]] <- readRDS(Filename)
}


## Fished Area
Age.Dist.F <- list()
Age.Catch.F <- list()

for(S in 1:4){
  
  Filename <- paste0(model.name, sep="_", "Sp_Population_F", sep="_", Scenarios[S], sep="_", "medium_movement")
  Age.Dist.F[[S]] <- readRDS(Filename)
  
  Filename <- paste0(model.name, sep="_", "Catch_by_Cell", sep="_", Scenarios[S], sep="_", "medium_movement")
  Age.Catch.F[[S]] <- readRDS(Filename)
}


### Format Zone data
res.age <- zone.fish.age(NTZ.Ages = Age.Dist.NTZ, F.Ages=Age.Dist.F)

Age.Dist.NTZ <- res.age[[1]]
Age.Dist.F <- res.age[[2]]

## Age catches

res.catch <- age.catch(Catches = Age.Catch.F, NTZ.Cells = shallow_NTZ_ID, F.Cells = shallow_F_ID, nsim=200) #This file actually has all of the cells for the whole model

Age.Catch.NTZ <- res.catch[[1]]
Age.Catch.F <- res.catch[[2]]


## Calculate median biomass of fish 
# Can turn the Zones on and off so you can do it for the whole mode only if you want 
Names <- c("Historical and Current NTZs", "Neither NTZs nor Temporal Management", 
           "Temporal and Spatial Management","Temporal Management Only")

medians.all <- median.biomass.func(Age = Age.Dist, Catch = Age.Catch, Zones = TRUE, 
                       Age.NTZ = Age.Dist.NTZ, Age.F = Age.Dist.F, 
                       #Catch.NTZ = Age.Catch.NTZ, Catch.F = Age.Catch.F,
                       nsim=200, nyears=59, yearstart=26, maxage=30,
                       scenario.names = Names, mat=Mature, kg=Weight)

Biomass.All <- medians.all[[1]]
Biomass.F <- medians.all[[2]] 
Biomass.NTZ <- medians.all[[3]] 

dim(Biomass.All[[2]])
##Create data frames needed to make the Kobe plots 
# 
# Kobe.Data.All <- NULL
# Kobe.Data.NTZ <- NULL
# Kobe.Data.F <- NULL
# 
# for(S in 1:4){
#   
#   temp <- FM.Scenario.NTZ[[S]]
#   temp2 <- Biomass.All[[S]]
#   
#   # temp3 <- temp2 %>% 
#   #   colSums() %>% 
#   #   as.numeric()
# 
#   
#   temp <- temp %>% 
#     as.data.frame() %>% 
#     rename(FM = ".") %>% 
#     mutate(Biomass = temp2) %>% 
#     mutate(Scenario = Names[S]) %>% 
#     mutate(Year = seq(1987,2018, 1))
#   
#   Kobe.Data.All <- rbind(Kobe.Data.All, temp)
#   
#   ## NTZ
#   temp <- FM.Scenario.NTZ[[S]]
#   temp2 <- Biomass.NTZ[[S]]
#   
#   temp3 <- temp2 %>% 
#     colSums() %>% 
#     as.numeric
#   
#   temp <- temp %>% 
#     as.data.frame() %>% 
#     rename(FM = ".") %>% 
#     mutate(Biomass = temp3) %>% 
#     mutate(Scenario = Names[S]) %>% 
#     mutate(Year = seq(1960,2018, 1))
#   
#   Kobe.Data.NTZ <- rbind(Kobe.Data.NTZ, temp)
#   
#   ## Fished
#   temp <- FM.Scenario.F[[S]]
#   temp2 <- Biomass.F[[S]]
#   
#   temp3 <- temp2 %>% 
#     colSums() %>% 
#     as.numeric()
#   
#   temp <- temp %>% 
#     as.data.frame() %>% 
#     rename(FM = ".") %>% 
#     mutate(Biomass = temp3) %>% 
#     mutate(Scenario = Names[S]) %>% 
#     mutate(Year = seq(1960,2018, 1))
#   
#   Kobe.Data.F <- rbind(Kobe.Data.F, temp)
# 
# }
# 
# Kobe.Data.All <- Kobe.Data.All %>% 
#   mutate(Area = "All")
# Kobe.Data.NTZ <- Kobe.Data.NTZ %>% 
#   mutate(Area = "NTZ") %>% 
#   mutate(FM = ifelse(Scenario %in% c("Historical and Current Management"),  MortRate$Fishing_Mort, FM)) %>% 
#   mutate(FM = ifelse(Scenario %in% c("Temporal and Spatial Management"),  MortRate$Fishing_Mort, FM))
# Kobe.Data.F <- Kobe.Data.F %>% 
#   mutate(Area = "F")

#* USING ACTUAL EFFORT FROM MODEL ####
setwd(sg_dir)

Effort_Scen <- list()
Spatial_Qs <- list()

setwd(sg_dir)
Effort_Scen[[1]] <- readRDS(paste0(model.name, sep="_", "fishing"))
Spatial_Qs[[1]] <- readRDS(paste0(model.name, sep="_", "Spatial_q_NTZ"))
Spatial_Qs[[2]] <- readRDS(paste0(model.name, sep="_", "Spatial_q_No_NTZ"))

setwd(sim_dir)
Effort_Scen[[2]] <- readRDS(paste0(model.name, sep="_", "S01_fishing"))
Effort_Scen[[3]] <- readRDS(paste0(model.name, sep="_", "S02_fishing"))
Effort_Scen[[4]] <- readRDS(paste0(model.name, sep="_", "S03_fishing"))

#* Convert back to effort in boat days ####

Boat_Days <- list()

Boat_Days_Scen <- array(0, dim=c(NCELL, 12, 59))

for(S in 1:4){
  
  if(S==1|S==3){
    for (YEAR in 1:59){
      Boat_Days_Scen[,,YEAR] <- Effort_Scen[[S]][,,YEAR] / Spatial_Qs[[1]][,YEAR]
    }
  } else {
    for (YEAR in 1:59){
      Boat_Days_Scen[,,YEAR] <- Effort_Scen[[S]][,,YEAR] / Spatial_Qs[[2]][,YEAR]
    }
  }
  Boat_Days_Scen[is.nan(Boat_Days_Scen)] <- 0
  
  Boat_Days[[S]] <- Boat_Days_Scen
}

Boat_Days_sum <- NULL

for(S in 1:4){
  temp <- Boat_Days[[S]]
  temp2 <- colSums(temp, dim=2) %>% 
    as.data.frame() %>% 
    mutate(Scenario = Names[S])
  
  Boat_Days_sum <- rbind(Boat_Days_sum, temp2)
}

Boat_Days_sum <- Boat_Days_sum %>% 
  rename(Effort = ".") %>% 
  mutate(Mortality = Effort*(0.00001)) %>% 
  mutate(Finite.M = 1-exp(-Mortality)) %>% 
  mutate(Year = rep(seq(1960,2018,1),4)) %>% 
  filter(Year>1984)

# #* Sum up effort for every year inside and outside sanctuary zones in every scenario
# 
# Boat_Days_sum_NTZ <- NULL
# Boat_Days_sum_F <- NULL
# 
# for(S in 1:4){
#   temp <- Boat_Days[[S]]
#   temp2 <- temp[as.numeric(shallow_NTZ_ID),,]
#   temp3 <- colSums(temp2, dim=2) %>% 
#     as.data.frame() %>% 
#     mutate(Scenario = Names[S])
#   
#   Boat_Days_sum_NTZ <- rbind(Boat_Days_sum_NTZ, temp3)
#   
#   temp <- Boat_Days[[S]]
#   temp2 <- temp[as.numeric(shallow_F_ID),,]
#   temp3 <- colSums(temp2, dim=2) %>% 
#     as.data.frame() %>% 
#     mutate(Scenario = Names[S])
#   
#   Boat_Days_sum_F <- rbind(Boat_Days_sum_F, temp3)
#   
# }

# Boat_Days_sum_NTZ <- Boat_Days_sum_NTZ %>% 
#   rename(Effort = ".")
# Boat_Days_sum_F <- Boat_Days_sum_F %>% 
#   rename(Effort = ".")
# 
# Boat_Days_sum_F <- Boat_Days_sum_F %>% 
#   mutate(Instant.M = Effort * 0.000005) %>% 
#   mutate(Finite.M = 1-exp(-Instant.M))
# 
# Boat_Days_sum_NTZ <- Boat_Days_sum_NTZ %>% 
#   mutate(Instant.M = Effort * 0.000005) %>% 
#   mutate(Finite.M = 1-exp(-Instant.M)) %>% 
#   mutate(Finite.M = ifelse(Scenario %in% c("Historical and Current Management"),  MortRate$Fishing_Mort, Finite.M)) %>% 
#   mutate(Finite.M = ifelse(Scenario %in% c("Temporal and Spatial Management"),  MortRate$Fishing_Mort, Finite.M))

Fmsy = 0.1
Bmsy.All = 25368.661
Bmsy.NTZ = 9181.556
Bmsy.F = Bmsy.All - Bmsy.NTZ

Kobe.Data.F <- as.data.frame(array(0, dim=c(nrow(Biomass.All), 7))) %>% 
  rename(year.F = "V1",
         Bio = "V2",
         Mod.Year = "V3",
         Scenario = "V4",
         Zone = "V5",
         Rel.F = "V6",
         Rel.Bio = "V7") %>% 
  mutate(Mod.Year = rep(seq(1985,2018,1),4),
         year.F = Boat_Days_sum$Finite.M,
         Bio = (Biomass.All$Median_MatBio - Biomass.NTZ$Median_MatBio)) %>% 
  mutate(Bio = ifelse(Mod.Year<1987, Biomass.All$Median_MatBio, Bio)) %>% 
  mutate(Rel.F = year.F/Fmsy) %>% 
  mutate(Rel.Bio = Bio/Bmsy.F) %>% 
  mutate(Rel.Bio = ifelse(Mod.Year<1987, Bio/Bmsy.All, Rel.Bio)) %>% 
  mutate(Scenario = rep(Names, each=34)) %>% 
  mutate(Zone = ifelse(Mod.Year < 1987, "All", "F"))

Kobe.Data.NTZ <- as.data.frame(array(0, dim=c(nrow(Biomass.All), 7))) %>% 
  rename(year.F = "V1",
         Bio = "V2",
         Mod.Year = "V3",
         Scenario = "V4",
         Zone = "V5",
         Rel.F = "V6",
         Rel.Bio = "V7") %>% 
  mutate(Mod.Year = rep(seq(1985,2018,1),4),
         year.F = Boat_Days_sum$Finite.M,
         Bio = Biomass.NTZ$Median_MatBio) %>% 
  mutate(Bio = ifelse(Mod.Year<1987, Biomass.All$Median_MatBio, Bio)) %>% 
  mutate(Rel.F = year.F/Fmsy) %>% 
  mutate(Rel.Bio = Bio/Bmsy.NTZ) %>% 
  mutate(Rel.Bio = ifelse(Mod.Year<1987, Bio/Bmsy.All, Rel.Bio)) %>% 
  mutate(Scenario = rep(Names, each=34)) %>% 
  mutate(Zone = ifelse(Mod.Year < 1987, "All", "NTZ"))


#### MAKE KOBE PLOTS ####

Kobe.Data <- rbind(Kobe.Data.F, Kobe.Data.NTZ[3:nrow(Kobe.Data.NTZ), ])

## Scenario 1
Kobe.Data.S00 <- Kobe.Data %>% 
  filter(Scenario %in% c("Historical and Current NTZs")) %>% 
  mutate(Mod.Year = as.numeric(Mod.Year)) %>% 
  mutate(Mod.Year = ifelse(Mod.Year %% 5 !=0, "", Mod.Year)) 

rownames(Kobe.Data.S00) <- NULL

Bio.Plot.S00 <- ggplot()+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 4,fill = "red")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 1, ymax = 4,fill = "orange")+
  geom_point(data=Kobe.Data.S00[1:2,], aes(x=Rel.Bio, y=Rel.F, group=Zone))+
  geom_path(data=Kobe.Data.S00[1:2,], aes(x=Rel.Bio, y=Rel.F, group=Zone))+
  # geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_point(data=Kobe.Data.S00[2:34,], aes(x=Rel.Bio, y=Rel.F))+
  geom_path(data=Kobe.Data.S00[2:34,], aes(x=Rel.Bio, y=Rel.F), linetype="dashed")+ # Fished
  geom_point(data=Kobe.Data.S00[c(2,35:66),], aes(x=Rel.Bio, y=Rel.F))+ 
  geom_path(data=Kobe.Data.S00[c(2,35:66),], aes(x=Rel.Bio, y=Rel.F))+ 
  geom_text(data=Kobe.Data.S00[c(1:3, 35:66), ], aes(x=Rel.Bio, y=Rel.F, group=Zone, label=Mod.Year), hjust=-0.25, vjust=-0.75, position=position_jitter(width=0.01,height=0.006))+
  scale_y_continuous(breaks=seq(0,4, 0.25), limits=c(0,4)) +
  scale_x_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
  theme_classic()+
  theme_classic()+
  ylab("F/Fmsy")+
  xlab("B/Bmsy")
Bio.Plot.S00 

## S01
Kobe.Data.S01 <- Kobe.Data %>% 
  filter(Scenario %in% c("Neither NTZs nor Temporal Management")) %>% 
  mutate(Mod.Year = as.numeric(Mod.Year)) %>% 
  mutate(Mod.Year = ifelse(Mod.Year %% 5 !=0, "", Mod.Year)) 

rownames(Kobe.Data.S01) <- NULL

Bio.Plot.S01 <-  ggplot()+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 4,fill = "red")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 1, ymax = 4,fill = "orange")+
  geom_point(data=Kobe.Data.S01[1:2,], aes(x=Rel.Bio, y=Rel.F, group=Zone))+
  geom_path(data=Kobe.Data.S01[1:2,], aes(x=Rel.Bio, y=Rel.F, group=Zone))+
  # geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_point(data=Kobe.Data.S01[2:34,], aes(x=Rel.Bio, y=Rel.F))+
  geom_path(data=Kobe.Data.S01[2:34,], aes(x=Rel.Bio, y=Rel.F), linetype="dashed")+ # Fished
  geom_point(data=Kobe.Data.S01[c(2,35:66),], aes(x=Rel.Bio, y=Rel.F))+ 
  geom_path(data=Kobe.Data.S01[c(2,35:66),], aes(x=Rel.Bio, y=Rel.F))+ 
  geom_text(data=Kobe.Data.S01[c(1:3,36:66), ], aes(x=Rel.Bio, y=Rel.F, group=Zone, label=Mod.Year), hjust=-0.25, vjust=-0.75, position=position_jitter(width=0.01,height=0.006))+
  scale_y_continuous(breaks=seq(0,4, 0.25), limits=c(0,4)) +
  scale_x_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
  theme_classic()+
  theme_classic()+
  ylab("F/FMSY")+
  xlab("B/BMSY")
Bio.Plot.S01

## S02
Kobe.Data.S02 <- Kobe.Data %>% 
  filter(Scenario %in% c("Temporal Management Only")) %>% 
  mutate(Mod.Year = as.numeric(Mod.Year)) %>% 
  mutate(Mod.Year = ifelse(Mod.Year %% 5 !=0, "", Mod.Year)) 

rownames(Kobe.Data.S02) <- NULL

Bio.Plot.S02 <-  ggplot()+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 3.5,fill = "red")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 1, ymax = 3.5,fill = "orange")+
  geom_point(data=Kobe.Data.S02[1:2,], aes(x=Rel.Bio, y=Rel.F, group=Zone))+
  geom_path(data=Kobe.Data.S02[1:2,], aes(x=Rel.Bio, y=Rel.F, group=Zone))+
  # geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_point(data=Kobe.Data.S02[2:34,], aes(x=Rel.Bio, y=Rel.F))+
  geom_path(data=Kobe.Data.S02[2:34,], aes(x=Rel.Bio, y=Rel.F), linetype="dashed")+ # Fished
  geom_point(data=Kobe.Data.S02[c(2,35:66),], aes(x=Rel.Bio, y=Rel.F))+ 
  geom_path(data=Kobe.Data.S02[c(2,35:66),], aes(x=Rel.Bio, y=Rel.F))+ 
  geom_text(data=Kobe.Data.S02[c(1:3,36:66), ], aes(x=Rel.Bio, y=Rel.F, group=Zone, label=Mod.Year), hjust=-0.25, vjust=-0.75, position=position_jitter(width=0.01,height=0.006))+
  scale_y_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
  scale_x_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
  theme_classic()+
  theme_classic()+
  ylab("F/Fmsy")+
  xlab("B/Bmsy")
Bio.Plot.S02

## S03
Kobe.Data.S03 <- Kobe.Data %>% 
  filter(Scenario %in% c("Temporal and Spatial Management")) %>% 
  mutate(Mod.Year = as.numeric(Mod.Year)) %>% 
  mutate(Mod.Year = ifelse(Mod.Year %% 5 !=0, "", Mod.Year)) 

rownames(Kobe.Data.S03) <- NULL

Bio.Plot.S03 <-  ggplot()+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 3.5,fill = "red")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 1, ymax = 3.5,fill = "orange")+
  geom_point(data=Kobe.Data.S03[1:2,], aes(x=Rel.Bio, y=Rel.F, group=Zone))+
  geom_path(data=Kobe.Data.S03[1:2,], aes(x=Rel.Bio, y=Rel.F, group=Zone))+
  # geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_point(data=Kobe.Data.S03[2:34,], aes(x=Rel.Bio, y=Rel.F))+
  geom_path(data=Kobe.Data.S03[2:34,], aes(x=Rel.Bio, y=Rel.F), linetype="dashed")+ # Fished
  geom_point(data=Kobe.Data.S03[c(2,35:66),], aes(x=Rel.Bio, y=Rel.F))+ 
  geom_path(data=Kobe.Data.S03[c(2,35:66),], aes(x=Rel.Bio, y=Rel.F))+ 
  geom_text(data=Kobe.Data.S03[c(1:3,36:66), ], aes(x=Rel.Bio, y=Rel.F, group=Zone, label=Mod.Year), hjust=-0.25, vjust=-0.75, position=position_jitter(width=0.01,height=0.006))+
  scale_y_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
  scale_x_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
  theme_classic()+
  ylab("F/Fmsy")+
  xlab("B/Bmsy")
Bio.Plot.S03

# ## All Shallow WHA
# # Have to rerun things to get here
# Kobe.Data.All <- Kobe.Data.All %>% 
#   mutate(Biomass = Biomass/1000) %>% 
#   # mutate(Biomass = Kobe.Data.NTZ$Biomass+Kobe.Data.F$Biomass) %>% 
#   # mutate(FM = ifelse(Scenario %in% c("Historical and Current Management"), Fishing.Mortality[,1],
#   #                    ifelse(Scenario %in% c("No Spatial Management"), Fishing.Mortality[,2],
#   #                           ifelse(Scenario %in% c("Temporal Management Only"), Fishing.Mortality[,3], Fishing.Mortality[,4])))) %>% 
#   mutate(Rel.FM = as.numeric(FM/F.MSY.All)) %>% 
#   mutate(Rel.Bio = as.numeric(Biomass/Bio.MSY.All))%>% 
#   # mutate(Rel.SB = as.numeric(Spawning.Biomass/SB.MSY)) %>% 
#   # filter(Year <= 1988) %>% 
#   rename(Mod.Year = "Year") %>% 
#   mutate(Mod.Year = as.character(Mod.Year))
# 
# Kobe.Data.All <- Kobe.Data %>% 
#   filter(Scenario %in% c("Historical and Current Management")) %>% 
#   mutate(Mod.Year = as.numeric(Mod.Year)) %>% 
#   mutate(Mod.Year = ifelse(Mod.Year %% 5 !=0, "", Mod.Year)) 
# 
# rownames(Kobe.Data.All) <- NULL
# 
# Bio.Plot.All <-  ggplot()+
#   annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
#   annotate("rect", xmin = 1, xmax = 3.5, ymin = 0, ymax = 1,fill = "green")+
#   annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 3.5,fill = "red")+
#   annotate("rect", xmin = 1, xmax = 3.5, ymin = 1, ymax = 3.5,fill = "orange")+
#   geom_point(data=Kobe.Data.S03[1:59,], aes(x=Rel.Bio, y=Rel.FM, group=Area))+
#   geom_path(data=Kobe.Data.S03[1:59,], aes(x=Rel.Bio, y=Rel.FM, group=Area))+
#   # # geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
#   # geom_point(data=Kobe.Data.S03[29:60,], aes(x=Rel.Bio, y=Rel.FM))+
#   # geom_path(data=Kobe.Data.S03[29:60,], aes(x=Rel.Bio, y=Rel.FM))+
#   # geom_point(data=Kobe.Data.S03[c(29,61:91),], aes(x=Rel.Bio, y=Rel.FM))+
#   # geom_path(data=Kobe.Data.S03[c(29,61:91),], aes(x=Rel.Bio, y=Rel.FM), linetype="dashed")+ # Fished
#   geom_text(data=Kobe.Data.S03[1:59, ], aes(x=Rel.Bio, y=Rel.FM, group=Area, label=Mod.Year), hjust=-0.25, vjust=-0.75, position=position_jitter(width=0.01,height=0.006))+
#   scale_y_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
#   scale_x_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
#   theme_classic()+
#   theme_classic()+
#   ylab("F/Fmsy")+
#   xlab("B/Bmsy")
# Bio.Plot.All


#### CHECKING KOBE PLOT CORRESPONDS TO CATCH ####
yearly.catch <- NULL

for(S in 1:4){
  
  temp <- Catches[[S]]
  
  temp2 <- temp %>% 
    colSums(.) %>% 
    as.data.frame() %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Mod.Year = seq(1960, 2018, 1)) %>% 
    rename(Catch = ".")
  
  yearly.catch <- rbind(yearly.catch, temp2)
  
}

catch.plot <- yearly.catch %>% 
  filter(Scenario %in% c("Historical and Current Management")) %>% 
  ggplot()+
  geom_line(aes(x=Mod.Year, y=Catch))+
  theme_classic()+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")

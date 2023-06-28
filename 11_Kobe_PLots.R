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
water <- readRDS(paste0(model.name, sep="_","water"))
YearlyTotal <- readRDS(paste0(model.name, sep="_", "BurnInPop"))
Selectivity <- readRDS("selret")
Mature <- readRDS("maturity")
Weight <- readRDS("weight")
MortRate <- readRDS("Mort_Rate")


NCELL <- nrow(water)


#### CALCULATE F AND BIOMASS FROM MY MODEL DATA ####
setwd(pop_dir)

#### SET UP FOR NTZ AND FISHED AREA ####
setwd(sp_dir)
bathy <- raster("ga_bathy_ningaloocrop.tif")
WHA <- st_read("2013_02_WorldHeritageMarineProgramme.shp") %>% 
  st_transform(4283)%>%
  st_make_valid %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) 


#* Create list of cells to restrict plots to shallow Water (<30m)
Water_points <- st_centroid_within_poly(water) 

Water_bathy <- raster::extract(bathy, Water_points, fun=mean, df=TRUE)

Water_bathy <- Water_bathy %>% 
  mutate(ID = as.factor(ID))

model_WHA <- water %>% 
  st_intersects(., WHA) %>% 
  as.data.frame()

Water_WHA <-water[c(as.numeric(model_WHA$row.id)), ]

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

# Each layer is a simulation, rows are cells and columns are years
Age.Dist[[1]] <- readRDS(paste0(model.name, sep="_", "Age_Distribution", sep="_", "S00"))
Age.Dist[[2]] <- readRDS(paste0(model.name, sep="_", "Age_Distribution", sep="_", "S01"))
Age.Dist[[3]] <- readRDS(paste0(model.name, sep="_", "Age_Distribution", sep="_", "S02"))
Age.Dist[[4]] <- readRDS(paste0(model.name, sep="_", "Age_Distribution", sep="_", "S03"))

Age.Catch[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Age", sep="_", "S00"))
Age.Catch[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Age", sep="_", "S01"))
Age.Catch[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Age", sep="_", "S02"))
Age.Catch[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Age", sep="_", "S03"))

## NTZ Area
Age.Dist.NTZ <- list()
Age.Catch.NTZ <- list()

# List of lists - 100 items in list, each item in list is rows as cells, columns as ages and layers are years
Age.Dist.NTZ[[1]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S00"))
Age.Dist.NTZ[[2]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S01"))
Age.Dist.NTZ[[3]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S02"))
Age.Dist.NTZ[[4]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S03"))


# Each layer is a simulation, rows are cells and columns are years
Age.Catch.NTZ[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S00"))
Age.Catch.NTZ[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S01"))
Age.Catch.NTZ[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S02"))
Age.Catch.NTZ[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S03"))

## Fished Area
Age.Dist.F <- list()
Age.Catch.F <- list()

# List of lists - 100 items in list, each item in list is rows as cells, columns as ages and layers are years
Age.Dist.F[[1]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S00"))
Age.Dist.F[[2]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S01"))
Age.Dist.F[[3]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S02"))
Age.Dist.F[[4]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S03"))


# Each layer is a simulation, rows are cells and columns are years
Age.Catch.F[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S00"))
Age.Catch.F[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S01"))
Age.Catch.F[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S02"))
Age.Catch.F[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S03"))


### Format Zone data

NTZ_Ages <- array(0, dim=c(30,59,100))
F_Ages <- array(0, dim=c(30,59,100))

## Numbers of fish by age
for(S in 1:4){
  
  SP_Pop_NTZ <- Age.Dist.NTZ[[S]]
  SP_Pop_F <- Age.Dist.F[[S]]
  
  for(SIM in 1:length(SP_Pop_NTZ)){
    
    temp <- as.data.frame(colSums(SP_Pop_NTZ[[SIM]])) %>% 
      mutate(Age = seq(1:30)) %>% 
      pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
      mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", ""))) %>% 
      mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
      group_by(Age, Mod_Year) %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      ungroup() %>%
      dplyr::select(!Num_Year) %>%
      pivot_wider(names_from="Mod_Year",id_col="Age", values_from="Number",values_fn = list(count=list)) %>% 
      dplyr::select(!Age) %>% 
      unlist()
    
    temp2 <- array(temp, dim=c(30,59))

    NTZ_Ages[,,SIM] <- temp2
    
    temp <- as.data.frame(colSums(SP_Pop_F[[SIM]])) %>% 
      mutate(Age = seq(1:30)) %>% 
      pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
      mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", ""))) %>% 
      mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
      group_by(Age, Mod_Year) %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      ungroup() %>% 
      pivot_wider(names_from="Mod_Year",id_col="Age", values_from="Number") %>% 
      dplyr::select(!Age) %>% 
      unlist()
   
     temp2 <- array(temp, dim=c(30,59))
    
    F_Ages[,,SIM] <- temp2
  }
  
  Age.Dist.NTZ[[S]] <- NTZ_Ages
  Age.Dist.F[[S]] <- F_Ages
  
}

## Age catches
catch <- array(0, dim=c(NCELL,59,100))

for (S in 1:4){
  
  temp <- Age.Catch.NTZ[[S]] # NTZ and F are the same files 
  
  for(SIM in 1:100){
    catch[,,SIM] <- temp[[SIM]]
  }
  
  temp2 <- catch[as.numeric(shallow_NTZ_ID),,]
  temp3 <- catch[as.numeric(shallow_F_ID),,]
  
  Age.Catch.NTZ[[S]] <- temp2
  
  Age.Catch.F[[S]] <- temp3
  
}

Names <- c("Historical and Current Management", "No Spatial Management", 
           "Temporal and Spatial Management","Temporal Management Only" )


## Calculate mean number of fish and biomass in each year of the model
Ages.All <- list()
Ages.F <- list()
Ages.NTZ <- list()

Catches.All <- list()
Catches.NTZ <- list()
Catches.F <- list()

catch <- array(0, dim=c(30,59,100))
catch.NTZ <- array(0, dim=c(30,59,100))
catch.F <- array(0, dim=c(30,59,100))

Biomass.All <- list()
Biomass.NTZ <- list()
Biomass.F <- list()

for(S in 1:4){
  
  #*************AGES***************
  temp.All <- Age.Dist[[S]]
  
  temp2 <- temp.All %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) 
  
  Ages.All[[S]] <- temp2
  
  # Ages NTZ
  temp.NTZ <- Age.Dist.NTZ[[S]]
  
  temp2 <- temp.NTZ %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) 
  
  Ages.NTZ[[S]] <- temp2
  
  # Ages F
  temp.F <- Age.Dist.F[[S]]
  
  temp2 <- temp.F %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) 
  
  Ages.F[[S]] <- temp2
  #*********************************
  #*************CATCHES*************
  ## All
  temp.c <- Age.Catch[[S]]
  for(SIM in 1:100){
    
    catch[,,SIM] <- temp.c[[SIM]]
    
  }
  
  temp2.c <- catch %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) 
  
  Catches.All[[S]] <- temp2.c
  
  ## NTZ
  temp.c.NTZ <- Age.Catch.NTZ[[S]]
  for(SIM in 1:100){
    
    catch.NTZ[,,SIM] <- temp.c.NTZ[[SIM]]
    
  }
  
  temp2.c.NTZ <- catch.NTZ %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) 
  
  if(S==1|S==3){
    temp2.c.NTZ[,] <- 0
  } else { }
  
  Catches.NTZ[[S]] <- temp2.c.NTZ
  
  ## F
  temp.c.F <- Age.Catch.F[[S]]
  for(SIM in 1:100){
    
    catch.F[,,SIM] <- temp.c.F[[SIM]]
    
  }
  
  temp2.c.F <- catch.F %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) 
  
  Catches.F[[S]] <- temp2.c.F
  
  #*********************************  
  #***********BIOMASS***************
  ## All
  temp3 <- temp.All %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) %>% 
    mutate_all(.,function(col){Weight$V12*col}) 
  
  Biomass.All[[S]] <- temp3
  
  ## NTZ
  temp3.NTZ <- temp.NTZ %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) %>% 
    mutate_all(.,function(col){Weight$V12*col}) 
  
  Biomass.NTZ[[S]] <- temp3.NTZ
  
  ## F
  temp3.F<- temp.F %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) %>% 
    mutate_all(.,function(col){Weight$V12*col}) 
  
  Biomass.F[[S]] <- temp3.F
  #*********************************  
}

## Biomass and Spawning Biomass
Kobe.Data.All <- NULL
Kobe.Data.NTZ <- NULL
Kobe.Data.F <- NULL

for(S in 1:4){
  
  temp <- FM.Scenario.All[[S]]
  temp2 <- Biomass.All[[S]]
  
  temp3 <- temp2 %>% 
    colSums() %>% 
    as.numeric()

  
  temp <- temp %>% 
    as.data.frame() %>% 
    rename(FM = ".") %>% 
    mutate(Biomass = temp3) %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Year = seq(1960,2018, 1))
  
  Kobe.Data.All <- rbind(Kobe.Data.All, temp)
  
  ## NTZ
  temp <- FM.Scenario.NTZ[[S]]
  temp2 <- Biomass.NTZ[[S]]
  
  temp3 <- temp2 %>% 
    colSums() %>% 
    as.numeric
  
  temp <- temp %>% 
    as.data.frame() %>% 
    rename(FM = ".") %>% 
    mutate(Biomass = temp3) %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Year = seq(1960,2018, 1))
  
  Kobe.Data.NTZ <- rbind(Kobe.Data.NTZ, temp)
  
  ## Fished
  temp <- FM.Scenario.F[[S]]
  temp2 <- Biomass.F[[S]]
  
  temp3 <- temp2 %>% 
    colSums() %>% 
    as.numeric()
  
  temp <- temp %>% 
    as.data.frame() %>% 
    rename(FM = ".") %>% 
    mutate(Biomass = temp3) %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Year = seq(1960,2018, 1))
  
  Kobe.Data.F <- rbind(Kobe.Data.F, temp)

}

Kobe.Data.All <- Kobe.Data.All %>% 
  mutate(Area = "All")
Kobe.Data.NTZ <- Kobe.Data.NTZ %>% 
  mutate(Area = "NTZ") %>% 
  mutate(FM = ifelse(Scenario %in% c("Historical and Current Management"),  MortRate$Fishing_Mort, FM)) %>% 
  mutate(FM = ifelse(Scenario %in% c("Temporal and Spatial Management"),  MortRate$Fishing_Mort, FM))
Kobe.Data.F <- Kobe.Data.F %>% 
  mutate(Area = "F")

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
  mutate(Mortality = Effort*(0.000005)) %>% 
  mutate(Finite.M = 1-exp(-Mortality))

#* Sum up effort for every year inside and outside sanctuary zones in every scenario

Boat_Days_sum_NTZ <- NULL
Boat_Days_sum_F <- NULL

for(S in 1:4){
  temp <- Boat_Days[[S]]
  temp2 <- temp[as.numeric(shallow_NTZ_ID),,]
  temp3 <- colSums(temp2, dim=2) %>% 
    as.data.frame() %>% 
    mutate(Scenario = Names[S])
  
  Boat_Days_sum_NTZ <- rbind(Boat_Days_sum_NTZ, temp3)
  
  temp <- Boat_Days[[S]]
  temp2 <- temp[as.numeric(shallow_F_ID),,]
  temp3 <- colSums(temp2, dim=2) %>% 
    as.data.frame() %>% 
    mutate(Scenario = Names[S])
  
  Boat_Days_sum_F <- rbind(Boat_Days_sum_F, temp3)
  
}

Boat_Days_sum_NTZ <- Boat_Days_sum_NTZ %>% 
  rename(Effort = ".")
Boat_Days_sum_F <- Boat_Days_sum_F %>% 
  rename(Effort = ".")

Boat_Days_sum_F <- Boat_Days_sum_F %>% 
  mutate(Instant.M = Effort * 0.000005) %>% 
  mutate(Finite.M = 1-exp(-Instant.M))

Boat_Days_sum_NTZ <- Boat_Days_sum_NTZ %>% 
  mutate(Instant.M = Effort * 0.000005) %>% 
  mutate(Finite.M = 1-exp(-Instant.M)) %>% 
  mutate(Finite.M = ifelse(Scenario %in% c("Historical and Current Management"),  MortRate$Fishing_Mort, Finite.M)) %>% 
  mutate(Finite.M = ifelse(Scenario %in% c("Temporal and Spatial Management"),  MortRate$Fishing_Mort, Finite.M))

Kobe.Data.All <- Kobe.Data.All %>% 
  mutate(Area = "All") %>% 
  mutate(FM = Boat_Days_sum$Finite.M)
Kobe.Data.NTZ <- Kobe.Data.NTZ %>% 
  mutate(Area = "NTZ") %>% 
  mutate(FM = ifelse(Scenario %in% c("Historical and Current Management"),  MortRate$Fishing_Mort, FM)) %>% 
  mutate(FM = ifelse(Scenario %in% c("Temporal and Spatial Management"),  MortRate$Fishing_Mort, FM))
Kobe.Data.F <- Kobe.Data.F %>% 
  mutate(Area = "F") %>% 
  mutate(FM = Boat_Days_sum_F$Finite.M)


#### MAKE KOBE PLOTS ####
# Currently set using the 33% of the biomass inside NTZ etc.
F.MSY.All = 0.1 
F.Cons.All = 0.09
Bio.MSY.All =  43653.160 # If we find out biomass in real model at this levels of F
Bio.Cons.All = 61.17121818

F.MSY.NTZ = 0.03
F.Cons.NTZ = 0.015
Bio.MSY.NTZ = 14405.54 # 7858.7788 If we find out biomass in real model at this levels of F
Bio.Cons.NTZ = 65.37503736

F.MSY.F = 0.03
F.Cons.F = 0.01
Bio.MSY.F = 28811.09 # 4127.5699 If we find out biomass in real model at this levels of F
Bio.Cons.F = 107.66974383


Kobe.Data.All.87 <- Kobe.Data.All %>% 
  # mutate(Biomass = Kobe.Data.NTZ$Biomass+Kobe.Data.F$Biomass) %>%
  #mutate(Biomass = Biomass/1000) %>%
  # mutate(FM = ifelse(Scenario %in% c("Historical and Current Management"), Fishing.Mortality[,1],
  #                    ifelse(Scenario %in% c("No Spatial Management"), Fishing.Mortality[,2],
  #                           ifelse(Scenario %in% c("Temporal Management Only"), Fishing.Mortality[,3], Fishing.Mortality[,4])))) %>%
  mutate(Rel.FM = as.numeric(FM/F.MSY.All)) %>% 
  mutate(Rel.Bio = as.numeric(Biomass/Bio.MSY.All))%>% 
  filter(Year <= 1988) %>% 
  rename(Mod.Year = "Year") %>% 
  mutate(Mod.Year = as.character(Mod.Year))

Kobe.Data.NTZ <- Kobe.Data.NTZ %>% 
  # mutate(Biomass = Biomass/1000) %>% 
  mutate(Rel.FM = as.numeric(FM/F.MSY.NTZ)) %>% 
  mutate(Rel.Bio = as.numeric(Biomass/Bio.MSY.NTZ)) %>% 
  filter(Year >= 1988) %>% 
  rename(Mod.Year = "Year") %>% 
  mutate(Mod.Year = as.character(Mod.Year))

Kobe.Data.F <- Kobe.Data.F %>% 
  # mutate(Biomass = Biomass/1000) %>% 
  mutate(FM = Boat_Days_sum_F$Finite.M) %>% 
  mutate(Rel.FM = as.numeric(FM/F.MSY.F)) %>% 
  mutate(Rel.Bio = as.numeric(Biomass/Bio.MSY.F))%>% 
  filter(Year >= 1988) %>% 
  rename(Mod.Year = "Year") %>% 
  mutate(Mod.Year = as.character(Mod.Year))

Kobe.Data <- rbind(Kobe.Data.All.87, Kobe.Data.NTZ, Kobe.Data.F)

## Scenario 1
Kobe.Data.S00 <- Kobe.Data %>% 
  filter(Scenario %in% c("Historical and Current Management")) %>% 
  mutate(Mod.Year = as.numeric(Mod.Year)) %>% 
  mutate(Mod.Year = ifelse(Mod.Year %% 5 !=0, "", Mod.Year)) 

rownames(Kobe.Data.S00) <- NULL

Bio.Plot.S00 <- ggplot()+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 3.5,fill = "red")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 1, ymax = 3.5,fill = "orange")+
  geom_point(data=Kobe.Data.S00[1:29,], aes(x=Rel.Bio, y=Rel.FM, group=Area))+
  geom_path(data=Kobe.Data.S00[1:29,], aes(x=Rel.Bio, y=Rel.FM, group=Area))+
  # geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_point(data=Kobe.Data.S00[29:60,], aes(x=Rel.Bio, y=Rel.FM))+
  geom_path(data=Kobe.Data.S00[29:60,], aes(x=Rel.Bio, y=Rel.FM))+
  geom_point(data=Kobe.Data.S00[c(29,61:91),], aes(x=Rel.Bio, y=Rel.FM))+
  geom_path(data=Kobe.Data.S00[c(29,61:91),], aes(x=Rel.Bio, y=Rel.FM), linetype="dashed")+ # Fished
  geom_text(data=Kobe.Data.S00[c(1:29,30:60,61:91), ], aes(x=Rel.Bio, y=Rel.FM, group=Area, label=Mod.Year), hjust=-0.25, vjust=-0.75, position=position_jitter(width=0.01,height=0.006))+
  scale_y_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
  scale_x_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
  theme_classic()+
  theme_classic()+
  ylab("F/Fmsy")+
  xlab("B/Bmsy")
Bio.Plot.S00 

## S01
Kobe.Data.S01 <- Kobe.Data %>% 
  filter(Scenario %in% c("No Spatial Management")) %>% 
  mutate(Mod.Year = as.numeric(Mod.Year)) %>% 
  mutate(Mod.Year = ifelse(Mod.Year %% 5 !=0, "", Mod.Year)) 

rownames(Kobe.Data.S01) <- NULL

Bio.Plot.S01 <-  ggplot()+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 3.5,fill = "red")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 1, ymax = 3.5,fill = "orange")+
  geom_point(data=Kobe.Data.S01[1:29,], aes(x=Rel.Bio, y=Rel.FM, group=Area))+
  geom_path(data=Kobe.Data.S01[1:29,], aes(x=Rel.Bio, y=Rel.FM, group=Area))+
  # geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_point(data=Kobe.Data.S01[29:60,], aes(x=Rel.Bio, y=Rel.FM))+
  geom_path(data=Kobe.Data.S01[29:60,], aes(x=Rel.Bio, y=Rel.FM))+
  geom_point(data=Kobe.Data.S01[c(29,61:91),], aes(x=Rel.Bio, y=Rel.FM))+
  geom_path(data=Kobe.Data.S01[c(29,61:91),], aes(x=Rel.Bio, y=Rel.FM), linetype="dashed")+ # Fished
  geom_text(data=Kobe.Data.S01[c(1:29,30:60,61:91), ], aes(x=Rel.Bio, y=Rel.FM, group=Area, label=Mod.Year), hjust=-0.25, vjust=-0.75, position=position_jitter(width=0.01,height=0.006))+
  scale_y_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
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
  geom_point(data=Kobe.Data.S02[1:29,], aes(x=Rel.Bio, y=Rel.FM, group=Area))+
  geom_path(data=Kobe.Data.S02[1:29,], aes(x=Rel.Bio, y=Rel.FM, group=Area))+
  # geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_point(data=Kobe.Data.S02[29:60,], aes(x=Rel.Bio, y=Rel.FM))+
  geom_path(data=Kobe.Data.S02[29:60,], aes(x=Rel.Bio, y=Rel.FM))+
  geom_point(data=Kobe.Data.S02[c(29,61:91),], aes(x=Rel.Bio, y=Rel.FM))+
  geom_path(data=Kobe.Data.S02[c(29,61:91),], aes(x=Rel.Bio, y=Rel.FM), linetype="dashed")+ # Fished
  geom_text(data=Kobe.Data.S02[c(1:29,30:60,61:91), ], aes(x=Rel.Bio, y=Rel.FM, group=Area, label=Mod.Year), hjust=-0.25, vjust=-0.75, position=position_jitter(width=0.01,height=0.006))+
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
  geom_point(data=Kobe.Data.S03[1:29,], aes(x=Rel.Bio, y=Rel.FM, group=Area))+
  geom_path(data=Kobe.Data.S03[1:29,], aes(x=Rel.Bio, y=Rel.FM, group=Area))+
  # geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_point(data=Kobe.Data.S03[29:60,], aes(x=Rel.Bio, y=Rel.FM))+
  geom_path(data=Kobe.Data.S03[29:60,], aes(x=Rel.Bio, y=Rel.FM))+
  geom_point(data=Kobe.Data.S03[c(29,61:91),], aes(x=Rel.Bio, y=Rel.FM))+
  geom_path(data=Kobe.Data.S03[c(29,61:91),], aes(x=Rel.Bio, y=Rel.FM), linetype="dashed")+ # Fished
  geom_text(data=Kobe.Data.S03[c(1:29,30:60,61:91), ], aes(x=Rel.Bio, y=Rel.FM, group=Area, label=Mod.Year), hjust=-0.25, vjust=-0.75, position=position_jitter(width=0.01,height=0.006))+
  scale_y_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
  scale_x_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
  theme_classic()+
  theme_classic()+
  ylab("F/Fmsy")+
  xlab("B/Bmsy")
Bio.Plot.S03

## All Shallow WHA
# Have to rerun things to get here
Kobe.Data.All <- Kobe.Data.All %>% 
  mutate(Biomass = Biomass/1000) %>% 
  # mutate(Biomass = Kobe.Data.NTZ$Biomass+Kobe.Data.F$Biomass) %>% 
  # mutate(FM = ifelse(Scenario %in% c("Historical and Current Management"), Fishing.Mortality[,1],
  #                    ifelse(Scenario %in% c("No Spatial Management"), Fishing.Mortality[,2],
  #                           ifelse(Scenario %in% c("Temporal Management Only"), Fishing.Mortality[,3], Fishing.Mortality[,4])))) %>% 
  mutate(Rel.FM = as.numeric(FM/F.MSY.All)) %>% 
  mutate(Rel.Bio = as.numeric(Biomass/Bio.MSY.All))%>% 
  # mutate(Rel.SB = as.numeric(Spawning.Biomass/SB.MSY)) %>% 
  # filter(Year <= 1988) %>% 
  rename(Mod.Year = "Year") %>% 
  mutate(Mod.Year = as.character(Mod.Year))

Kobe.Data.All <- Kobe.Data %>% 
  filter(Scenario %in% c("Historical and Current Management")) %>% 
  mutate(Mod.Year = as.numeric(Mod.Year)) %>% 
  mutate(Mod.Year = ifelse(Mod.Year %% 5 !=0, "", Mod.Year)) 

rownames(Kobe.Data.All) <- NULL

Bio.Plot.All <-  ggplot()+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 3.5,fill = "red")+
  annotate("rect", xmin = 1, xmax = 3.5, ymin = 1, ymax = 3.5,fill = "orange")+
  geom_point(data=Kobe.Data.S03[1:59,], aes(x=Rel.Bio, y=Rel.FM, group=Area))+
  geom_path(data=Kobe.Data.S03[1:59,], aes(x=Rel.Bio, y=Rel.FM, group=Area))+
  # # geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  # geom_point(data=Kobe.Data.S03[29:60,], aes(x=Rel.Bio, y=Rel.FM))+
  # geom_path(data=Kobe.Data.S03[29:60,], aes(x=Rel.Bio, y=Rel.FM))+
  # geom_point(data=Kobe.Data.S03[c(29,61:91),], aes(x=Rel.Bio, y=Rel.FM))+
  # geom_path(data=Kobe.Data.S03[c(29,61:91),], aes(x=Rel.Bio, y=Rel.FM), linetype="dashed")+ # Fished
  geom_text(data=Kobe.Data.S03[1:59, ], aes(x=Rel.Bio, y=Rel.FM, group=Area, label=Mod.Year), hjust=-0.25, vjust=-0.75, position=position_jitter(width=0.01,height=0.006))+
  scale_y_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
  scale_x_continuous(breaks=seq(0,3.5, 0.25), limits=c(0,3.5)) +
  theme_classic()+
  theme_classic()+
  ylab("F/Fmsy")+
  xlab("B/Bmsy")
Bio.Plot.All


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

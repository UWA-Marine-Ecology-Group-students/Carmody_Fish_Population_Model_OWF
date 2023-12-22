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
library(pwr)

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

Names <- c("Historical and Current NTZs", "Neither NTZs nor Temporal Management", 
           "Temporal and Spatial Management","Temporal Management Only")

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
total_pop <- total.pop.format.full(pop.file.list = total_pop_list, scenario.names = Names, nsim=200, nyears=59, startyear=26, maxage=30, mat = Mature, kg=Weight)

total_pop <- total_pop %>% 
  filter(Year %in% c(2018)) %>% 
  #filter(Scenario %in% c("Historical and Current NTZs", "Temporal Management Only")) %>% 
  mutate(Year = as.factor(Year),
         Scenario = as.factor(Scenario)) %>% 
  mutate(Movement = "Medium")

mod <- lm(log(MatBio) ~ Scenario, data = total_pop)
summary(mod)

pwr.f2.test(u = 2, v = 397, f2 = (0.02917/(1-0.02917)), sig.level = 0.05)
# If I've done this right my power is 0.88

total_pop %>% group_by(Scenario) %>% 
  summarise(mean = mean(MatBio),
            sd = sd(MatBio))

ggplot()+
  geom_boxplot(data=total_pop, aes(x=Scenario, y=log(MatBio)), notch=T)


summary(mod)
plot(mod$residuals)

mod.aov <- aov(mod)
mod.tukey <- TukeyHSD(mod.aov)


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

## S00 - Normal

NTZ.F.Ages.S00 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S00, f.list = SP_Pop_NTZ_S00, scenario.name = Names[1], nsim = 200)

NTZ.S00 <- NTZ.F.Ages.S00[[1]] %>% 
  mutate(Movement = "Medium")
F.S00 <- NTZ.F.Ages.S00[[2]]


## S01 - No NTZs (and no temporal closure)
NTZ.F.Ages.S01 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S01, f.list = SP_Pop_NTZ_S01, scenario.name = Names[2], nsim = 200)

NTZ.S01 <- NTZ.F.Ages.S01[[1]] %>% 
  mutate(Movement = "Medium")
F.S01 <- NTZ.F.Ages.S01[[2]]


## S02 - Temporal Closure, no NTZs
NTZ.F.Ages.S02 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S02, f.list = SP_Pop_NTZ_S02, scenario.name = Names[3], nsim = 200)

NTZ.S02 <- NTZ.F.Ages.S02[[1]] %>% 
  mutate(Movement = "Medium")
F.S02 <- NTZ.F.Ages.S02[[2]]


## S03 - Temporal Closure and NTZs
NTZ.F.Ages.S03 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S03, f.list = SP_Pop_NTZ_S03, scenario.name = Names[4], nsim = 200)

NTZ.S03 <- NTZ.F.Ages.S03[[1]] %>% 
  mutate(Movement = "Medium")
F.S03 <- NTZ.F.Ages.S03[[2]]


## Put everything together
Whole_Pop_Ages_NTZ <- rbind(NTZ.S00, NTZ.S01, NTZ.S02, NTZ.S03) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F.S00, F.S01, F.S02, F.S03) %>% 
  mutate(Zone = "F")

NTZ_data_Medium <- Whole_Pop_Ages_NTZ %>% 
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

mod1.NTZ.aov <- aov(mod1.NTZ)
mod1.NTZ.tukey <- TukeyHSD(mod1.NTZ.aov)
mod1.NTZ.tukey

F_Recruits <- F_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.F <- lm(log(Number) ~ Scenario, dat=F_Recruits)
summary(mod1.F)
plot(mod1.F$residuals)

mod1.F.aov <- aov(mod1.F)
mod1.F.tukey <- TukeyHSD(mod1.F.aov)
mod1.F.tukey

# Sublegal
NTZ_Sublegal <- NTZ_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.NTZ <- lm(log(Number) ~ Scenario, dat=NTZ_Sublegal)
summary(mod2.NTZ)
plot(mod2.NTZ$residuals)

mod2.NTZ.aov <- aov(mod2.NTZ)
mod2.NTZ.tukey <- TukeyHSD(mod2.NTZ.aov)
mod2.NTZ.tukey


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

mod3.NTZ <- lm(log(Number) ~ Scenario, dat=NTZ_Legal)
summary(mod3.NTZ)
plot(mod3.NTZ$residuals)

mod3.NTZ.aov <- aov(mod3.NTZ)
mod3.NTZ.tukey <- TukeyHSD(mod3.NTZ.aov)
mod3.NTZ.tukey


F_Legal <- F_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod3.F <- lm(log(Number) ~ Scenario, dat=F_Legal)
summary(mod3.F)
plot(mod3.F$residuals)

mod3.F.aov <- aov(mod3.F)
mod3.F.tukey <- TukeyHSD(mod3.F.aov)
mod3.F.tukey



# Large Legal
NTZ_Large <- NTZ_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario)) 

mod4.NTZ <- lm(log(Number) ~ Scenario, dat=NTZ_Large)
summary(mod4.NTZ)
plot(mod4.NTZ$residuals)

mod4.NTZ.aov <- aov(mod4.NTZ)
mod4.NTZ.tukey <- TukeyHSD(mod4.NTZ.aov)
mod4.NTZ.tukey


F_Large <- F_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod4.F <- lm(log(Number) ~ Scenario, dat=F_Large)
summary(mod4.F)
plot(mod4.F$residuals)

mod4.F.aov <- aov(mod4.F)
mod4.F.tukey <- TukeyHSD(mod4.F.aov)
mod4.F.tukey


#### DISTANCE ABUNDANCE AND CATCH GLMS ####
#* Work out which cells are within 100km of each of the boat ramps ####
setwd(sp_dir)
BR <- st_read("Boat_Ramps.shp") %>% 
  st_transform(4283)%>%
  st_make_valid() 
BR <- BR[1:4,]

network <- st_read(paste0(model.name, sep="_","network.shapefile.shp"))

setwd(sg_dir)

water <- readRDS(paste0(model.name, sep="_","water"))
BR <- st_as_sf(BR)
st_crs(BR) <- NA 

BR <- BR %>% 
  mutate(name = c("Bundegi","Exmouth","Tantabiddi","CoralBay"))

centroids <- st_centroid_within_poly(water)
points <- as.data.frame(st_coordinates(centroids))%>% #The points start at the bottom left and then work their way their way right
  mutate(ID=row_number()) 
points_sf <- st_as_sf(points, coords = c("X", "Y")) 
st_crs(points_sf) <- NA


network <- as_sfnetwork(network, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

net <- activate(network, "nodes")
network_matrix <- st_network_cost(net, from=BR, to=points_sf)
network_matrix <- network_matrix*111
dim(network_matrix)

DistBR <- as.data.frame(t(network_matrix)) %>% 
  rename("Bd_BR"=V1) %>% 
  rename("ExM_BR" = V2) %>% 
  rename("Tb_BR" = V3) %>% 
  rename("CrB_BR"=V4) %>% 
  mutate(CellID = row_number()) 

DistBR.WHA <- DistBR[c(as.numeric(model_WHA$row.id)), ]

NCELL.WHA <- nrow(DistBR.WHA)
BR.10km <- NULL
temp <- array(0, dim=c(NCELL.WHA,2)) %>% 
  as.data.frame() %>% 
  rename(.,"Distance" = V1,
         "CellID" = V2)

for(RAMP in 1:4){
  temp[,1:2] <- DistBR.WHA[,c(RAMP, 5)]
  
  Dist10 <- temp %>% 
    filter(Distance <=10)
  
  BR.10km <- rbind(BR.10km, Dist10)
}

BR.10km <- unique(BR.10km$CellID)

BR.50km <- NULL
temp <- array(0, dim=c(NCELL.WHA,2)) %>% 
  as.data.frame() %>% 
  rename(.,"Distance" = V1,
         "CellID" = V2)

for(RAMP in 1:4){
  temp[,1:2] <- DistBR.WHA[,c(RAMP, 5)]
  
  Dist50 <- temp %>% 
    filter(Distance>10 & Distance <=50)
  
  BR.50km <- rbind(BR.50km, Dist50)
}

BR.50km <- unique(BR.50km$CellID)

BR.100km <- NULL
temp <- array(0, dim=c(NCELL.WHA,2)) %>% 
  as.data.frame() %>% 
  rename(.,"Distance" = V1,
         "CellID" = V2)

for(RAMP in 1:4){
  temp[,1:2] <- DistBR.WHA[,c(RAMP, 5)]
  
  Dist100 <- temp %>% 
    filter(Distance >50 & Distance <=100)
  
  BR.100km <- rbind(BR.100km, Dist100)
}

BR.100km <- unique(BR.100km$CellID)

Distances <- list()

Distances[[1]] <- BR.10km
Distances[[2]] <- BR.50km
Distances[[3]] <- BR.100km

Dist.Names <- as.character(c("0-10 km", "10-50 km", "50-100 km"))

#* Abundance ####
setwd(pop_dir)

# Each layer is a simulation, rows are cells and columns are years
Pop.Dist <- list()

Pop.Dist[[1]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S00", sep="_", "medium_movement"))
Pop.Dist[[2]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S01", sep="_", "medium_movement"))
Pop.Dist[[3]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S02", sep="_", "medium_movement"))
Pop.Dist[[4]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S03", sep="_", "medium_movement"))

Names <- c("Historical and Current NTZs", "Neither NTZs nor Temporal Management", 
           "Temporal and Spatial Management","Temporal Management Only")

Dists.res <- dist.abundance.full(Pops = Pop.Dist, n.year=59, n.sim=200, n.cell=NCELL,
                                 n.scenario=4, fished.cells=water, distances=Distances,
                                 dist.names=Dist.Names, scen.names=Names, mod.years = seq(1960,2018,1), target.year = 59)


# S00.Abundance <- Dists.res[[1]]
# S01.Abundance <- Dists.res[[2]]
# S02.Abundance <- Dists.res[[3]]
# S03.Abundance <- Dists.res[[4]]

# Function is doing something really weird so I've done it manually by running the inside of the function

S00.Abundance <- Pop.Abundance[[1]]
S01.Abundance <- Pop.Abundance[[2]]
S02.Abundance <- Pop.Abundance[[3]]
S03.Abundance <- Pop.Abundance[[4]]

Dist.Abundance <- rbind(S00.Abundance, S01.Abundance, S02.Abundance, S03.Abundance)

## Abundance model 10km
Dist.Abundance.10km <- Dist.Abundance %>% 
  filter(Distance %in% "0-10 km") %>% 
  mutate(Scenario = as.factor(Scenario)) 

mod1.10km <- lm(log(Total) ~ Scenario, dat=Dist.Abundance.10km)
summary(mod1.10km)
plot(mod1.10km$residuals)

## Abundance model 50km
Dist.Abundance.50km <- Dist.Abundance %>% 
  filter(Distance %in% "10-50 km") %>% 
  mutate(Scenario = as.factor(Scenario)) 

mod1.50km <- lm(log(Total) ~ Scenario, dat=Dist.Abundance.50km)
summary(mod1.50km)
plot(mod1.10km$residuals)

## Abundance model 100km
Dist.Abundance.100km <- Dist.Abundance %>% 
  filter(Distance %in% "50-100 km") %>% 
  mutate(Scenario = as.factor(Scenario)) 

mod1.100km <- lm(log(Total) ~ Scenario, dat=Dist.Abundance.100km)
summary(mod1.100km)
plot(mod1.10km$residuals)
#* Catch ####
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

## Convert back to effort in boat days 

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

#* Sum up effort for every year in each of the scenarios
Effort_Dist <-list()

for(D in 1:3){
  
  Boat_Days_sum <- NULL
  
  for(S in 1:4){
    temp <- Boat_Days[[S]]
    temp2 <- temp[as.numeric(Distances[[D]]), ,] %>% 
      colSums(., dim=2) %>% 
      as.data.frame() %>% 
      mutate(Scenario = Names[S]) %>% 
      mutate(Distance = Dist.Names[D]) %>% 
      mutate(Year = seq(1960,2018,1))
    
    Boat_Days_sum <- rbind(Boat_Days_sum, temp2)
  }
  
  Boat_Days_sum <- Boat_Days_sum %>% 
    rename(Effort = ".")
  Effort_Dist[[D]] <- Boat_Days_sum
} 



setwd(pop_dir)

Pop.Catch <- list()

# Each layer is a simulation, rows are cells and columns are years
Pop.Catch[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S00", sep="_", "medium_movement"))
Pop.Catch[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S01", sep="_", "medium_movement"))
Pop.Catch[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S02", sep="_", "medium_movement"))
Pop.Catch[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S03", sep="_", "medium_movement"))

## NEED TO TURN INTO CPUE
Catch.res <- dist.catch.full(Pops = Pop.Catch, n.year=59, n.sim=200, n.cell=NCELL,
                                      n.scenario=4, fished.cells=water, distances=Distances,
                                       dist.names=Dist.Names, scen.names=Names, mod.years = seq(1960,2018,1), target.year=59)


# Pop.Catch.S00 <- Catch.res[[1]]
# Pop.Catch.S01 <- Catch.res[[2]]
# Pop.Catch.S02 <- Catch.res[[3]]
# Pop.Catch.S03 <- Catch.res[[4]]

## Function is running weirdly so I've run the inside of the function manually as I can't be bothered to fix it right now

Pop.Catch.S00 <- Pop.Catch[[1]]
Pop.Catch.S01 <- Pop.Catch[[2]]
Pop.Catch.S02 <- Pop.Catch[[3]]
Pop.Catch.S03 <- Pop.Catch[[4]]

Dist.Catch <- rbind(Pop.Catch.S00, Pop.Catch.S01, Pop.Catch.S02, Pop.Catch.S03)

## Catch model 10km
Effort.10km <- Effort_Dist[[1]]

Dist.Catch.10km <- Dist.Catch %>% 
  filter(Year %in% 2018) %>% 
  filter(Distance %in% "0-10 km") %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  left_join(., Effort.10km, by=c("Scenario", "Distance", "Year")) %>% 
  mutate(CPUE = Total/Effort)

mod1.10km <- lm(log(CPUE) ~ Scenario, dat=Dist.Catch.10km)
summary(mod1.10km)
plot(mod1.10km$residuals)

## Catch model 50km
Effort.50km <- Effort_Dist[[2]]

Dist.Catch.50km <- Dist.Catch %>% 
  filter(Year %in% 2018) %>% 
  filter(Distance %in% "10-50 km") %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  left_join(., Effort.50km, by=c("Scenario", "Distance", "Year")) %>% 
  mutate(CPUE = Total/Effort)

mod2.50km <- lm(log(CPUE) ~ Scenario, dat=Dist.Catch.50km)
summary(mod2.50km)
plot(mod2.50km$residuals)

## Catch model 100km
Effort.100km <- Effort_Dist[[3]]

Dist.Catch.100km <- Dist.Catch %>% 
  filter(Year %in% 2018) %>% 
  filter(Distance %in% "50-100 km") %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  left_join(., Effort.100km, by=c("Scenario", "Distance", "Year")) %>% 
  mutate(CPUE = Total/Effort)

mod3.100km <- lm(log(CPUE) ~ Scenario, dat=Dist.Catch.100km)
summary(mod3.100km)
plot(mod3.100km$residuals)

#### MOVEMENT SCENARIO GLMS ####
setwd(pop_dir)

#* Slow ####
total_pop_list_slow <- list()

total_pop_list_slow[[1]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S00_slow_movement"))
total_pop_list_slow[[2]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S01_slow_movement"))
total_pop_list_slow[[3]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S02_slow_movement"))
total_pop_list_slow[[4]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S03_slow_movement"))

# NEED TO TURN OFF THE BIT OF THE FUNCITON THAT CREATES THE MEDIANS AS YOU WANT ALL OF THE SEPARATE RUNS
total_pop_slow <- total.pop.format.full(pop.file.list = total_pop_list, scenario.names = Names, nsim=100, nyears=59, startyear=26, maxage=30, mat = Mature, kg=Weight)

total_pop_slow <- total_pop_slow %>% 
  filter(Year %in% c(2018)) %>% 
  mutate(Year = as.factor(Year),
         Scenario = as.factor(Scenario)) %>% 
  mutate(Movement = "Slow")

mod.slow <- lm(log(MatBio) ~ Scenario, data=total_pop_slow)
summary(mod.slow)
plot(mod.slow$residuals)

mod.slow.aov <- aov(mod.slow)
tukey.slow <- TukeyHSD(mod.slow.aov)
tukey.slow 

#* Fast ####
total_pop_list_fast <- list()

total_pop_list_fast[[1]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S00_fast_movement"))
total_pop_list_fast[[2]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S01_fast_movement"))
total_pop_list_fast[[3]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S02_fast_movement"))
total_pop_list_fast[[4]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S03_fast_movement"))

total_pop_fast <- total.pop.format.full(pop.file.list = total_pop_list_fast, scenario.names = Names, nsim=100, nyears=59, startyear=26, maxage=30, mat = Mature, kg=Weight)

total_pop_fast <- total_pop_fast %>% 
  filter(Year %in% c(2018)) %>% 
  mutate(Year = as.factor(Year),
         Scenario = as.factor(Scenario)) %>% 
  mutate(Movement = "Fast")

mod.fast <- lm(log(MatBio) ~ Scenario, data=total_pop_fast)
summary(mod.fast)
plot(mod.fast$residuals)

## Comparing same scenario but different movement
total_pop_movement <- rbind(total_pop, total_pop_slow, total_pop_fast)

S00_movement <- total_pop_movement %>% 
  filter(Scenario %in% "Historical and Current NTZs")

mod.S00 <- lm(log(MatBio) ~ Movement, data=S00_movement)

mod.S00.aov <- aov(mod.S00)
tukey.S00 <- TukeyHSD(mod.S00.aov)
tukey.S00 

S01_movement <- total_pop_movement %>% 
  filter(Scenario %in% "Neither NTZs nor Temporal Management")

mod.S01 <- lm(log(MatBio) ~ Movement, data=S01_movement)

mod.S01.aov <- aov(mod.S01)
tukey.S01 <- TukeyHSD(mod.S01.aov)
tukey.S01 

S02_movement <- total_pop_movement %>% 
  filter(Scenario %in% "Temporal Management Only")

mod.S02 <- lm(log(MatBio) ~ Movement, data=S02_movement)

mod.S02.aov <- aov(mod.S02)
tukey.S02 <- TukeyHSD(mod.S02.aov)
tukey.S02 

S03_movement <- total_pop_movement %>% 
  filter(Scenario %in% "Temporal and Spatial Management")

mod.S03 <- lm(log(MatBio) ~ Movement, data=S03_movement)

mod.S03.aov <- aov(mod.S03)
tukey.S03 <- TukeyHSD(mod.S03.aov)
tukey.S03 


# #* ZONE LINEAR MODELS SLOW 
# setwd(pop_dir)
# SP_Pop_NTZ_S00 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_slow_movement"))
# SP_Pop_F_S00 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_fast_movement"))
# 
# SP_Pop_NTZ_S01 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_slow_movement"))
# SP_Pop_F_S01 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_fast_movement"))
# 
# SP_Pop_NTZ_S02 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_slow_movement"))
# SP_Pop_F_S02 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_slow_movement"))
# 
# SP_Pop_NTZ_S03 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_slow_movement"))
# SP_Pop_F_S03 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_slow_movement"))
# 


#* ZONE LINEAR MODELS - SLOW ####
setwd(pop_dir)
SP_Pop_NTZ_S00 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_slow_movement"))
SP_Pop_F_S00 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_slow_movement"))

SP_Pop_NTZ_S01 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_slow_movement"))
SP_Pop_F_S01 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_slow_movement"))

SP_Pop_NTZ_S02 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_slow_movement"))
SP_Pop_F_S02 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_slow_movement"))

SP_Pop_NTZ_S03 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_slow_movement"))
SP_Pop_F_S03 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_slow_movement"))


#* Format zone data ####

## HAVE TURNED OFF THE PART OF THE FUNCTION THAT CREATES THE MEDIANS ##
## S00 - Normal

NTZ.F.Ages.S00 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S00, f.list = SP_Pop_NTZ_S00, scenario.name = Names[1], nsim = 100)

NTZ.S00 <- NTZ.F.Ages.S00[[1]] %>% 
  mutate(Movement = "Slow")
F.S00 <- NTZ.F.Ages.S00[[2]]


## S01 - No NTZs (and no temporal closure)
NTZ.F.Ages.S01 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S01, f.list = SP_Pop_NTZ_S01, scenario.name = Names[2], nsim = 100)

NTZ.S01 <- NTZ.F.Ages.S01[[1]] %>% 
  mutate(Movement = "Slow")
F.S01 <- NTZ.F.Ages.S01[[2]]


## S02 - Temporal Closure, no NTZs
NTZ.F.Ages.S02 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S02, f.list = SP_Pop_NTZ_S02, scenario.name = Names[3], nsim = 100)

NTZ.S02 <- NTZ.F.Ages.S02[[1]] %>% 
  mutate(Movement = "Slow")
F.S02 <- NTZ.F.Ages.S02[[2]]


## S03 - Temporal Closure and NTZs
NTZ.F.Ages.S03 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S03, f.list = SP_Pop_NTZ_S03, scenario.name = Names[4], nsim = 100)

NTZ.S03 <- NTZ.F.Ages.S03[[1]] %>% 
  mutate(Movement = "Slow")
F.S03 <- NTZ.F.Ages.S03[[2]]


## Put everything together
Whole_Pop_Ages_NTZ <- rbind(NTZ.S00, NTZ.S01, NTZ.S02, NTZ.S03) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F.S00, F.S01, F.S02, F.S03) %>% 
  mutate(Zone = "F")

NTZ_data_slow <- Whole_Pop_Ages_NTZ %>% 
  pivot_longer(cols=starts_with("V"), names_to = "Simulation",values_to="Number")

F_data_slow <- Whole_Pop_Ages_F %>% 
  pivot_longer(cols=starts_with("V"), names_to = "Simulation",values_to="Number")


## NTZ models 
# Recruits 
NTZ_Recruits_slow <- NTZ_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.NTZ.slow <- lm(log(Number) ~ Scenario, dat=NTZ_Recruits_slow)
summary(mod1.NTZ.slow)
plot(mod1.NTZ.slow$residuals)

mod1.NTZ.slow.aov <- aov(mod1.NTZ.slow)
mod1.NTZ.slow.tukey <- TukeyHSD(mod1.NTZ.slow.aov)
mod1.NTZ.slow.tukey

F_Recruits_slow <- F_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.F.slow <- lm(log(Number) ~ Scenario, dat=F_Recruits_slow)
summary(mod1.F.slow)
plot(mod1.F.slow$residuals)

mod1.F.slow.aov <- aov(mod1.F.slow)
mod1.F.slow.tukey <- TukeyHSD(mod1.F.slow.aov)
mod1.F.slow.tukey

# Sublegal
NTZ_Sublegal_slow <- NTZ_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.NTZ.slow <- lm(log(Number) ~ Scenario, dat=NTZ_Sublegal_slow)
summary(mod2.NTZ.slow)
plot(mod2.NTZ.slow$residuals)

mod2.NTZ.slow.aov <- aov(mod2.NTZ.slow)
mod2.NTZ.slow.tukey <- TukeyHSD(mod2.NTZ.slow.aov)
mod2.NTZ.slow.tukey


F_Sublegal_slow <- F_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.F.slow <- lm(log(Number) ~ Scenario, dat=F_Sublegal_slow)
summary(mod2.F.slow)
plot(mod2.F.slow$residuals)

mod2.F.slow.aov <- aov(mod2.F.slow)
mod2.F.slow.tukey <- TukeyHSD(mod2.F.slow.aov)
mod2.F.slow.tukey

# Legal
NTZ_Legal_slow <- NTZ_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod3.NTZ.slow <- lm(log(Number) ~ Scenario, dat=NTZ_Legal_slow)
summary(mod3.NTZ.slow)
plot(mod3.NTZ.slow$residuals)

mod3.NTZ.slow.aov <- aov(mod3.NTZ.slow)
mod3.NTZ.slow.tukey <- TukeyHSD(mod3.NTZ.slow.aov)
mod3.NTZ.slow.tukey


F_Legal_slow <- F_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod3.F.slow <- lm(log(Number) ~ Scenario, dat=F_Legal_slow)
summary(mod3.F.slow)
plot(mod3.F.slow$residuals)

mod3.F.slow.aov <- aov(mod3.F.slow)
mod3.F.slow.tukey <- TukeyHSD(mod3.F.slow.aov)
mod3.F.slow.tukey



# Large Legal
NTZ_Large_slow <- NTZ_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario)) 

mod4.NTZ.slow <- lm(log(Number) ~ Scenario, dat=NTZ_Large_slow)
summary(mod4.NTZ.slow)
plot(mod4.NTZ.slow$residuals)

mod4.NTZ.slow.aov <- aov(mod4.NTZ.slow)
mod4.NTZ.slow.tukey <- TukeyHSD(mod4.NTZ.slow.aov)
mod4.NTZ.slow.tukey


F_Large_slow <- F_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod4.F.slow <- lm(log(Number) ~ Scenario, dat=F_Large_slow)
summary(mod4.F.slow)
plot(mod4.F.slow$residuals)

mod4.F.slow.aov <- aov(mod4.F.slow)
mod4.F.slow.tukey <- TukeyHSD(mod4.F.slow.aov)
mod4.F.slow.tukey

#* ZONE LINEAR MODELS - FAST ####
setwd(pop_dir)
SP_Pop_NTZ_S00 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_fast_movement"))
SP_Pop_F_S00 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_fast_movement"))

SP_Pop_NTZ_S01 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_fast_movement"))
SP_Pop_F_S01 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_fast_movement"))

SP_Pop_NTZ_S02 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_fast_movement"))
SP_Pop_F_S02 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_fast_movement"))

SP_Pop_NTZ_S03 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_fast_movement"))
SP_Pop_F_S03 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_fast_movement"))


#* Format zone data ####

## S00 - Normal

NTZ.F.Ages.S00 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S00, f.list = SP_Pop_NTZ_S00, scenario.name = Names[1], nsim = 100)

NTZ.S00 <- NTZ.F.Ages.S00[[1]] %>% 
  mutate(Movement = "Fast")
F.S00 <- NTZ.F.Ages.S00[[2]]


## S01 - No NTZs (and no temporal closure)
NTZ.F.Ages.S01 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S01, f.list = SP_Pop_NTZ_S01, scenario.name = Names[2], nsim = 100)

NTZ.S01 <- NTZ.F.Ages.S01[[1]] %>% 
  mutate(Movement = "Fast")
F.S01 <- NTZ.F.Ages.S01[[2]]


## S02 - Temporal Closure, no NTZs
NTZ.F.Ages.S02 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S02, f.list = SP_Pop_NTZ_S02, scenario.name = Names[3], nsim = 100)

NTZ.S02 <- NTZ.F.Ages.S02[[1]] %>% 
  mutate(Movement = "Fast")
F.S02 <- NTZ.F.Ages.S02[[2]]


## S03 - Temporal Closure and NTZs
NTZ.F.Ages.S03 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S03, f.list = SP_Pop_NTZ_S03, scenario.name = Names[4], nsim = 100)

NTZ.S03 <- NTZ.F.Ages.S03[[1]] %>% 
  mutate(Movement = "Fast")
F.S03 <- NTZ.F.Ages.S03[[2]]


## Put everything together
Whole_Pop_Ages_NTZ <- rbind(NTZ.S00, NTZ.S01, NTZ.S02, NTZ.S03) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F.S00, F.S01, F.S02, F.S03) %>% 
  mutate(Zone = "F")

NTZ_data_fast <- Whole_Pop_Ages_NTZ %>% 
  pivot_longer(cols=starts_with("V"), names_to = "Simulation",values_to="Number")

F_data_fast <- Whole_Pop_Ages_F %>% 
  pivot_longer(cols=starts_with("V"), names_to = "Simulation",values_to="Number")


## NTZ models 
# Recruits 
NTZ_Recruits_fast <- NTZ_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.NTZ.fast <- lm(log(Number) ~ Scenario, dat=NTZ_Recruits_fast)
summary(mod1.NTZ.fast)
plot(mod1.NTZ.fast$residuals)

mod1.NTZ.fast.aov <- aov(mod1.NTZ.fast)
mod1.NTZ.fast.tukey <- TukeyHSD(mod1.NTZ.fast.aov)
mod1.NTZ.fast.tukey

F_Recruits_fast <- F_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.F.fast <- lm(log(Number) ~ Scenario, dat=F_Recruits_fast)
summary(mod1.F.fast)
plot(mod1.F.fast$residuals)

mod1.F.fast.aov <- aov(mod1.F.fast)
mod1.F.fast.tukey <- TukeyHSD(mod1.F.fast.aov)
mod1.F.fast.tukey

# Sublegal
NTZ_Sublegal_fast <- NTZ_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.NTZ.fast <- lm(log(Number) ~ Scenario, dat=NTZ_Sublegal_fast)
summary(mod2.NTZ.fast)
plot(mod2.NTZ.fast$residuals)

mod2.NTZ.fast.aov <- aov(mod2.NTZ.fast)
mod2.NTZ.fast.tukey <- TukeyHSD(mod2.NTZ.fast.aov)
mod2.NTZ.fast.tukey


F_Sublegal_fast <- F_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.F.fast <- lm(log(Number) ~ Scenario, dat=F_Sublegal_fast)
summary(mod2.F.fast)
plot(mod2.F.fast$residuals)

mod2.F.fast.aov <- aov(mod2.F.fast)
mod2.F.fast.tukey <- TukeyHSD(mod2.F.fast.aov)
mod2.F.fast.tukey

# Legal
NTZ_Legal_fast <- NTZ_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod3.NTZ.fast <- lm(log(Number) ~ Scenario, dat=NTZ_Legal_fast)
summary(mod3.NTZ.fast)
plot(mod3.NTZ.fast$residuals)

mod3.NTZ.fast.aov <- aov(mod3.NTZ.fast)
mod3.NTZ.fast.tukey <- TukeyHSD(mod3.NTZ.fast.aov)
mod3.NTZ.fast.tukey


F_Legal_fast <- F_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod3.F.fast <- lm(log(Number) ~ Scenario, dat=F_Legal_fast)
summary(mod3.F.fast)
plot(mod3.F.fast$residuals)

mod3.F.fast.aov <- aov(mod3.F.fast)
mod3.F.fast.tukey <- TukeyHSD(mod3.F.fast.aov)
mod3.F.fast.tukey



# Large Legal
NTZ_Large_fast <- NTZ_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario)) 

mod4.NTZ.fast <- lm(log(Number) ~ Scenario, dat=NTZ_Large_fast)
summary(mod4.NTZ.fast)
plot(mod4.NTZ.fast$residuals)

mod4.NTZ.fast.aov <- aov(mod4.NTZ.fast)
mod4.NTZ.fast.tukey <- TukeyHSD(mod4.NTZ.fast.aov)
mod4.NTZ.fast.tukey


F_Large_fast <- F_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod4.F.fast <- lm(log(Number) ~ Scenario, dat=F_Large_fast)
summary(mod4.F.fast)
plot(mod4.F.fast$residuals)

mod4.F.fast.aov <- aov(mod4.F.fast)
mod4.F.fast.tukey <- TukeyHSD(mod4.F.fast.aov)
mod4.F.fast.tukey

#* Comparing movement across scenarios ####
All_movement_NTZ <- rbind(NTZ_data_Medium, NTZ_data_slow, NTZ_data_fast) %>% 
  filter(Mod_Year %in% 2018)

S00_movement <- All_movement_NTZ %>% 
  filter(Scenario %in% "Historical and Current NTZs") %>% 
  filter(Stage %in% "Large Legal")

mod.S00 <- lm(log(Number) ~ Movement, data=S00_movement)

mod.S00.aov <- aov(mod.S00)
tukey.S00 <- TukeyHSD(mod.S00.aov)
tukey.S00 

S01_movement <- All_movement_NTZ %>% 
  filter(Scenario %in% "Neither NTZs nor Temporal Management") %>% 
  filter(Stage %in% "Large Legal")

mod.S01 <- lm(log(Number) ~ Movement, data=S01_movement)

mod.S01.aov <- aov(mod.S01)
tukey.S01 <- TukeyHSD(mod.S01.aov)
tukey.S01 

S02_movement <- All_movement_NTZ %>% 
  filter(Scenario %in% "Temporal Management Only") %>% 
  filter(Stage %in% "Large Legal")

mod.S02 <- lm(log(Number) ~ Movement, data=S02_movement)

mod.S02.aov <- aov(mod.S02)
tukey.S02 <- TukeyHSD(mod.S02.aov)
tukey.S02 

S03_movement <- All_movement_NTZ %>% 
  filter(Scenario %in% "Temporal and Spatial Management")%>% 
  filter(Stage %in% "Large Legal")

mod.S03 <- lm(log(Number) ~ Movement, data=S03_movement)

mod.S03.aov <- aov(mod.S03)
tukey.S03 <- TukeyHSD(mod.S03.aov)
tukey.S03 













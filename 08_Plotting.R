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


# Normal
# Sim 1 Nothing
# Sim 2 NTZs and Temporal Closure
# Sim 3 Just temporal closure, no sanctuary zones 

model.name <- "ningaloo"
#### READ IN DATA ####
setwd(sg_dir)
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))


setwd(sp_dir)
water <- readRDS(paste0(model.name, sep="_","water"))
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

shallow_cells_NTZ <- water %>% 
  filter(Fished_2017 %in% c("N")) %>% 
  st_drop_geometry(.)

shallow_cells_F <- water %>% 
  filter(Fished_2017 %in% c("Y")) %>% 
  st_drop_geometry(.)

shallow_NTZ_ID <- as.numeric(shallow_cells_NTZ$ID)
shallow_F_ID <- as.numeric(shallow_cells_F$ID)

#### CALCULATE TOTAL AREA OF FISHED AND NO-TAKE ####

AreaFished <- Water_shallow %>% 
  mutate(cell_area = st_area(Spatial)) %>% 
  mutate(cell_area = as.numeric(cell_area)) %>% 
  mutate(Fished=as.factor(Fished_2017)) %>% 
  filter(Fished=="Y") 

AreaFished <- (sum(AreaFished$cell_area))/100000

AreaNT <- Water_shallow %>% 
  mutate(cell_area = st_area(Spatial)) %>% 
  mutate(cell_area = as.numeric(cell_area)) %>% 
  mutate(Fished=as.factor(Fished_2017)) %>% 
  filter(Fished=="N") 

AreaNT <- (sum(AreaNT$cell_area))/100000

#* Normal - RE-DONE ####
setwd(pop_dir)

# Whole Population 
TotalPop <- array(0, dim=c(59,2))
TotalPop <- as.data.frame(TotalPop)
TotalPop <- TotalPop %>% 
  rename(Year="V1") %>% 
  mutate(Year=seq(1960,2018,1)) %>% 
  rename(Tot.Pop="V2")

numYear <- seq(0,58,1)

for(Y in 1:59){
  
  year <- readRDS(paste0(model.name, sep="_", "Rcpp_YearlyTotal_", numYear[Y])) %>% 
    unlist()
  year <- array(year, dim=c(1834,12,30))
  year <- rowSums(year[,,1:30], dim=2)
  year <- sum(year[,12])
  
  TotalPop[Y,2] <- year
  
}

# Separated by Fished and No-Take
NoTakeAges <- array(0, dim=c(30,13))
NoTakeAges <- as.data.frame(NoTakeAges)
FishedAges <- NoTakeAges


numYears <- seq(0,59,5)
numYears[1] <- 0
numYears[13] <- 58

for(YEAR in 1:13){
  
  Population <-  readRDS(paste0(model.name, sep="_", "Rcpp_YearlyTotal_", numYears[YEAR])) %>% 
    unlist()
  Population <- array(Population, dim=c(1834, 12,30))
  
  Population.NT <- Population[c(as.numeric(shallow_NTZ_ID)),12, ] %>% 
    colSums(.)
  Population.F <- Population[c(as.numeric(shallow_F_ID)),12, ] %>% 
    colSums(.)
  
  NoTakeAges[,YEAR] <- Population.NT
  FishedAges[,YEAR] <- Population.F
  
  if (YEAR==13){
    NoTakeAges <- NoTakeAges %>% 
      mutate(Status = "NTZ") %>% 
      mutate(Age = seq(1:30)) %>% 
      mutate(Scenario = "Normal") %>% 
      mutate(Scenario = as.factor(Scenario))
    colnames(NoTakeAges) <- c("1960", "1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018", "Status", "Age", "Scenario")
    
    FishedAges <- FishedAges %>% 
      mutate(Status = "Fished") %>% 
      mutate(Age = seq(1:30)) %>% 
      mutate(Scenario = "Normal") %>% 
      mutate(Scenario = as.factor(Scenario))
    colnames(FishedAges) <- c("1960", "1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018", "Status", "Age", "Scenario")
    
    WholePop <- rbind(NoTakeAges, FishedAges) %>%  
      pivot_longer(cols=-c(Age, Status, Scenario), names_to="Year", values_to="Number") %>% 
      mutate(Year = as.factor(Year)) %>% 
      mutate(Year = fct_relevel(Year, c("1960", "1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018"))) %>% 
      mutate(Status = as.factor(Status)) 
    
    NoRecruits <- WholePop %>% 
      filter(Age!=1) 
  }
  
}

#* Scenario 1 - RE-DONE ####
setwd(sim_dir)

# Whole Population 
s1_TotalPop <- array(0, dim=c(59,2))
s1_TotalPop <- as.data.frame(s1_TotalPop)
s1_TotalPop <- s1_TotalPop %>% 
  rename(Year="V1") %>% 
  mutate(Year=seq(1960,2018,1)) %>% 
  rename(s1_Tot.Pop="V2")

numYear <- seq(0,58,1)

for(Y in 1:59){
  
  year <- readRDS(paste0(model.name, sep="_", "S01_Rcpp_YearlyTotal_", numYear[Y])) %>% 
  unlist()
  year <- array(year, dim=c(1834,12,30))
  year <- rowSums(year[,,1:30], dim=2)
  year <- sum(year[,12])
  
  s1_TotalPop[Y,2] <- year
  
}


# Separated by Fished and No-Take
s1_NoTakeAges <- array(0, dim=c(30,13))
s1_NoTakeAges <- as.data.frame(s1_NoTakeAges)
s1_FishedAges <- s1_NoTakeAges


numYears <- seq(0,59,5)
numYears[1] <- 0
numYears[13] <- 58

for(YEAR in 1:13){
  
  Population <-  readRDS(paste0(model.name, sep="_", "S01_Rcpp_YearlyTotal_", numYears[YEAR]))%>% 
    unlist()
  Population <- array(Population, dim=c(1834, 12,30))
  
  Population.NT <- Population[c(as.numeric(shallow_NTZ_ID)),12, ] %>% 
    colSums(.)
  Population.F <- Population[-c(as.numeric(shallow_NTZ_ID)),12, ] %>% 
    colSums(.)
  
  s1_NoTakeAges[,YEAR] <- Population.NT
  s1_FishedAges[,YEAR] <- Population.F
  
  if (YEAR==13){
    s1_NoTakeAges <- s1_NoTakeAges %>% 
      mutate(Status = "NTZ") %>% 
      mutate(Age = seq(1:30)) %>% 
      mutate(Scenario = "Nothing") %>% 
      mutate(Scenario = as.factor(Scenario))
    colnames(s1_NoTakeAges) <- c("1960", "1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018", "Status", "Age", "Scenario")
    
    s1_FishedAges <- s1_FishedAges %>% 
      mutate(Status = "Fished") %>% 
      mutate(Age = seq(1:30)) %>% 
      mutate(Scenario = "Nothing") %>% 
      mutate(Scenario = as.factor(Scenario))
    colnames(s1_FishedAges) <- c("1960", "1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018", "Status", "Age", "Scenario")
    
    s1_WholePop <- rbind(s1_NoTakeAges, s1_FishedAges) %>%  
      pivot_longer(cols=-c(Age, Status, Scenario), names_to="Year", values_to="Number") %>% 
      mutate(Year = as.factor(Year)) %>% 
      mutate(Year = fct_relevel(Year, c("1960", "1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018"))) %>% 
      mutate(Status = as.factor(Status)) 
    
    s1_NoRecruits <- s1_WholePop %>% 
      filter(Age!=1) 
  }
  
}


#* Scenario 2 - Temp Closure with NTZ RE-DONE ####
# Whole Population 
s2_TotalPop <- array(0, dim=c(59,2))
s2_TotalPop <- as.data.frame(s2_TotalPop)
s2_TotalPop <- s2_TotalPop %>% 
  rename(Year="V1") %>% 
  mutate(Year=seq(1960,2018,1)) %>% 
  rename(s2_Tot.Pop="V2")

numYear <- seq(0,58,1)

for(Y in 1:59){
  
  year <- readRDS(paste0(model.name, sep="_", "S02_Rcpp_YearlyTotal_", numYear[Y])) %>% 
  unlist()
  year <- array(year, dim=c(1834,12,30))
  year <- rowSums(year[,,1:30], dim=2)
  year <- sum(year[,12])
  s2_TotalPop[Y,2] <- year
  
}


# Separated by Fished and No-Take
s2_NoTakeAges <- array(0, dim=c(30,13))
s2_NoTakeAges <- as.data.frame(s2_NoTakeAges)
s2_FishedAges <- s2_NoTakeAges


numYears <- seq(0,59,5)
numYears[1] <- 0
numYears[13] <- 58

for(YEAR in 1:13){
  
  Population <-  readRDS(paste0(model.name, sep="_", "S02_Rcpp_YearlyTotal_", numYears[YEAR]))%>% 
    unlist()
  Population <- array(Population, dim=c(1834, 12,30))
  
  Population.NT <- Population[c(as.numeric(shallow_NTZ_ID)),12, ] %>% 
    colSums(.)
  Population.F <- Population[-c(as.numeric(shallow_NTZ_ID)),12, ] %>% 
    colSums(.)
  
  s2_NoTakeAges[,YEAR] <- Population.NT
  s2_FishedAges[,YEAR] <- Population.F
  
  if (YEAR==13){
    s2_NoTakeAges <- s2_NoTakeAges %>% 
      mutate(Status = "NTZ") %>% 
      mutate(Age = seq(1:30)) %>% 
      mutate(Scenario = "Temp Closure and NTZ") %>% 
      mutate(Scenario = as.factor(Scenario))
    colnames(s2_NoTakeAges) <- c("1960", "1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018", "Status", "Age", "Scenario")
    
    s2_FishedAges <- s2_FishedAges %>% 
      mutate(Status = "Fished") %>% 
      mutate(Age = seq(1:30)) %>% 
      mutate(Scenario = "Temp Closure and NTZ") %>% 
      mutate(Scenario = as.factor(Scenario))
    colnames(s2_FishedAges) <- c("1960", "1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018", "Status", "Age", "Scenario")
    
    s2_WholePop <- rbind(s2_NoTakeAges, s2_FishedAges) %>%  
      pivot_longer(cols=-c(Age, Status, Scenario), names_to="Year", values_to="Number") %>% 
      mutate(Year = as.factor(Year)) %>% 
      mutate(Year = fct_relevel(Year, c("1960", "1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018"))) %>% 
      mutate(Status = as.factor(Status)) 
    
    s2_NoRecruits <- s2_WholePop %>% 
      filter(Age!=1) 
  }
  
}

#* Scenario 3 - Just Temp Closure COMPLETED ####
# Whole Population 
s3_TotalPop <- array(0, dim=c(59,2))
s3_TotalPop <- as.data.frame(s3_TotalPop)
s3_TotalPop <- s3_TotalPop %>% 
  rename(Year="V1") %>% 
  mutate(Year=seq(1960,2018,1)) %>% 
  rename(s3_Tot.Pop="V2")

numYear <- seq(0,58,1)

for(Y in 1:59){
  
  year <- readRDS(paste0(model.name, sep="_", "S03_Rcpp_YearlyTotal_", numYear[Y]))%>% 
    unlist()
  year <- array(year, dim=c(1834,12,30))
  year <- rowSums(year[,,1:30], dim=2)
  year <- sum(year[,12])
  
  s3_TotalPop[Y,2] <- year
  
}


# Separated by Fished and No-Take
s3_NoTakeAges <- array(0, dim=c(30,13))
s3_NoTakeAges <- as.data.frame(s3_NoTakeAges)
s3_FishedAges <- s3_NoTakeAges


numYears <- seq(0,59,5)
numYears[1] <- 0
numYears[13] <- 58

for(YEAR in 1:13){
  
  Population <-  readRDS(paste0(model.name, sep="_", "S03_Rcpp_YearlyTotal_", numYears[YEAR]))%>% 
    unlist()
  Population <- array(Population, dim=c(1834, 12,30))
  
  Population.NT <- Population[c(as.numeric(shallow_NTZ_ID)),12, ] %>% 
    colSums(.)
  Population.F <- Population[-c(as.numeric(shallow_NTZ_ID)),12, ] %>% 
    colSums(.)
  
  s3_NoTakeAges[,YEAR] <- Population.NT
  s3_FishedAges[,YEAR] <- Population.F
  
  if (YEAR==13){
    s3_NoTakeAges <- s3_NoTakeAges %>% 
      mutate(Status = "NTZ") %>% 
      mutate(Age = seq(1:30)) %>% 
      mutate(Scenario = "Temp Closure") %>% 
      mutate(Scenario = as.factor(Scenario))
    colnames(s3_NoTakeAges) <- c("1960", "1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018", "Status", "Age", "Scenario")
    
    s3_FishedAges <- s3_FishedAges %>% 
      mutate(Status = "Fished") %>% 
      mutate(Age = seq(1:30)) %>% 
      mutate(Scenario = "Temp Closure") %>% 
      mutate(Scenario = as.factor(Scenario))
    colnames(s3_FishedAges) <- c("1960", "1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018", "Status", "Age", "Scenario")
    
    s3_WholePop <- rbind(s3_NoTakeAges, s3_FishedAges) %>%  
      pivot_longer(cols=-c(Age, Status, Scenario), names_to="Year", values_to="Number") %>% 
      mutate(Year = as.factor(Year)) %>% 
      mutate(Year = fct_relevel(Year, c("1960", "1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018"))) %>% 
      mutate(Status = as.factor(Status)) 
    
    s3_NoRecruits <- s3_WholePop %>% 
      filter(Age!=1) 
  }
  
}


#### TIME SERIES PLOT WITH ALL SCENARIOS ####
options(repr.plot.width =9, repr.plot.height =9)

TotalPop <- reduce(list(TotalPop, s1_TotalPop, s2_TotalPop, s3_TotalPop), dplyr::left_join ,by='Year') %>% 
  dplyr::rename(Normal = "Tot.Pop",
         Nothing = "s1_Tot.Pop",
         Temp.Closure.NTZ = "s2_Tot.Pop",
         Temp.Closure = "s3_Tot.Pop") %>% 
  pivot_longer(cols=-c(Year), names_to="Scenario", values_to="Total.Population")

TimeSeries <- TotalPop %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre 1987", ifelse(Scenario %in% c("Normal") & Year>1985, "NTZs as normal", 
                                                             ifelse(Scenario %in% c("Temp.Closure") & Year>1985, "Temporal Closure Only", 
                                                                    ifelse(Scenario %in% c("Nothing"), "None", "Temporal Closure and NTZs"))))) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ColourGroup = fct_relevel(ColourGroup, c("Pre 1987", "NTZs as normal", "None" ,"Temporal Closure Only", "Temporal Closure and NTZs"))) %>% 
  ggplot()+
  geom_line(aes(x=Year, y=Total.Population, group=Scenario, colour=ColourGroup)) +
  scale_x_continuous("Year", breaks = c(1960, 1970, 1980, 1990, 2000, 2010, 2020))+
  scale_colour_manual("Spatial and Temporal\nManagement Scenario",values=c( "gray20",  "#69BE28", "#005594", "#8AD2D8", "#53AF8B"), labels=c("Pre-1987", "Historical and current\nmanagement", 
                                                                                                            "No spatial or temporal\nmanagement" ,"Temporal management\nonly", 
                                                                                                            "Spatial and temporal\nmanagement"))+
  xlab("Year")+
  ylab("Total Population")+
  geom_vline(xintercept=1987, linetype="dotted", color="grey20")+
  geom_vline(xintercept=2005, linetype="dashed", colour="grey20")+
  geom_vline(xintercept=2008, colour="grey20")+
  theme_classic()+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12)) + #change legend text font size
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
TimeSeries
  
#### LINE PLOTS BY AGE GROUP OF ALL SCENARIOS ####

AllNoTake <- rbind(NoTakeAges, s1_NoTakeAges, s2_NoTakeAges, s3_NoTakeAges) 
AllFished <- rbind(FishedAges, s1_FishedAges, s2_FishedAges, s3_FishedAges)

ScenarioWholePop <- rbind(AllNoTake, AllFished) %>%  
  pivot_longer(cols=-c(Age, Status, Scenario), names_to="Year", values_to="Number") %>% 
  mutate(NumKM2 = ifelse(Status %in% c("Fished"), Number/AreaFished, Number/AreaNT)) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<=4, "Sublegal",
                               ifelse(Age>4 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) 

mylabels1 <- expression("Pre-1987", NULL, NULL,NULL, NULL,NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL)
mylabels2 <- expression("Historical and\ncurrent management", "No spatial or temporal\nmanagement", "Temporal\nmanagement only", "Spatial and temporal\nmanagement",
                       NULL, NULL, NULL, NULL,NULL)
mylabels3 <- expression(NULL, NULL,NULL, NULL, NULL, NULL,NULL,NULL, NULL, NULL, NULL, NULL)


line.recruits <- ScenarioWholePop %>% 
  filter(Age==1) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Normal") & Year>1985, "Historical and current management", 
                                                             ifelse(Scenario %in% c("Temp Closure") & Year>1985, "Temporal management only", 
                                                                    ifelse(Scenario %in% c("Nothing"), "No spatial or temporal management", "Spatial and temporal management"))))) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  #filter(Status %in% c("NTZ")) %>% 
  #mutate(NumKM2 = ifelse(Year==1960 & Scenario %in% ("Nothing") & Stage %in% c("Recruit") & Status %in% c("NTZ"), 1.332263391, NumKM2)) %>% 
  mutate(PercChange = ifelse(Status %in% c("Fished"), ((NumKM2-1.436692)/1.436692)*100, ((NumKM2-1.905991)/1.905991)*100)) %>% 
  mutate(PercChange = ifelse(Year==1960, 0, PercChange)) %>% 
  ggplot(., aes(x=Year, y=NumKM2, group=interaction(Status,Scenario), colour=ColourGroup))+
  geom_point(size=0.1, aes(fill=ShapeGroup, shape=ShapeGroup, group=Status))+
  geom_point(size=0.1, aes(fill=ShapeGroup, shape=ShapeGroup, group=Stage))+
  geom_point(size=2.5, aes(fill=ShapeGroup,  shape=ShapeGroup, group=interaction(Status,Scenario)))+
  geom_line(aes(colour=ColourGroup, group=interaction(Status,Scenario)))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dotted", color="grey20")+
  geom_vline(xintercept=10, linetype="dashed", colour="grey20")+
  geom_vline(xintercept=12.5,colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ`='circle', `NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled", `NTZ.Temp Closure`="diamond filled", `NTZ.Temp Closure and NTZ`="triangle down filled" , `Pre-1987`="circle",
                               `Staged`="square",`Fished`="circle", `Fished.Normal`="square filled", `Fished.Nothing`="triangle filled", `Fished.Temp Closure`="diamond filled", `Fished.Temp Closure and NTZ`="triangle down filled"), name="Area and scenario",
                     labels=mylabels1, guide="none")+
  scale_fill_manual(values= c(`Pre-1987`="grey20", `NTZ`="white", `NTZ.Normal`="#69BE28", `NTZ.Nothing`="#005594", `NTZ.Temp Closure`="#8AD2D8", `NTZ.Temp Closure and NTZ`="#53AF8B",
                              `Staged`="white",`Fished`="white", `Fished.Normal`="white", `Fished.Nothing`="white",`Fished.Temp Closure`="white", `Fished.Temp Closure and NTZ`="white"),
                    labels=mylabels1, name="Spatial and Temporal\nManagement Scenario")+
  scale_colour_manual(values = c(`Pre-1987`="grey20", `NTZ`="white", `Historical and current management`="#69BE28", `Temporal management only`="#8AD2D8", `No spatial or temporal management`="#005594", 
                                 `Spatial and temporal management`="#53AF8B", `Fished`="white", `Stage`="white"), guide="none")+ 
  guides(fill = guide_legend(nrow = 6, label.position = "right",
                             override.aes = list(shape = c("circle", "bullet","square filled", "triangle filled", "diamond filled", "triangle down filled", "bullet", "bullet","square filled", "triangle filled", "diamond filled", "triangle down filled"),
                                                 colour=c("grey20", "white","white",  "white", "white", "white", "white", "white","white",  "white", "white", "white"),
                                                 fill = c("grey20", "white","white",  "white", "white", "white", "white", "white","white",  "white", "white", "white"),
                                                 labels = mylabels1)))+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0, 0.025)+
  theme(#legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=14, face="bold"), #change legend title font size
        legend.text = element_text(size=12)) + #change legend text font size
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=14,face="bold"))
  #ggplot2::annotate("text", x=1.4, y=1.99, label="(a) Recruits", size = 3, fontface=2)
line.recruits


line.sublegal <- ScenarioWholePop %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Normal") & Year>1985, "Historical and current management", 
                                                             ifelse(Scenario %in% c("Temp Closure") & Year>1985, "Temporal management only", 
                                                                    ifelse(Scenario %in% c("Nothing"), "No spatial or temporal management", "Spatial and temporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  #filter(Status %in% c("NTZ")) %>% 
  mutate(PercChange = ifelse(Status %in% c("Fished"), ((Total-0.1732865)/0.1732865)*100, ((Total-0.4809484)/0.4809484)*100)) %>% 
  ggplot(.)+
  geom_point(aes(x=Year, y=Total, group=interaction(Scenario, Status), colour=ColourGroup, fill=ShapeGroup,  shape=ShapeGroup), size=2.5,)+
  geom_line(aes(x=Year, y=Total, group=interaction(Scenario, Status), colour=ColourGroup))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dotted", color="grey20")+
  geom_vline(xintercept=10, linetype="dashed", colour="grey20")+
  geom_vline(xintercept=12.5,colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled", `NTZ.Temp Closure`="diamond filled", `NTZ.Temp Closure and NTZ`="triangle down filled" , `Pre-1987`="circle",
                               `Fished.Normal`="square filled", `Fished.Nothing`="triangle filled", `Fished.Temp Closure`="diamond filled", `Fished.Temp Closure and NTZ`="triangle down filled"), name="Area and scenario",
                     labels=mylabels2, guide="none")+
  scale_fill_manual(values= c(`Pre-1987`="grey20",  `NTZ.Normal`="#69BE28", `NTZ.Nothing`="#005594", `NTZ.Temp Closure`="#8AD2D8", `NTZ.Temp Closure and NTZ`="#53AF8B",
                              `Fished.Normal`="white", `Fished.Nothing`="white",`Fished.Temp Closure`="white", `Fished.Temp Closure and NTZ`="white"),
                    labels=mylabels3, name="No take\nzone")+
  scale_colour_manual(values = c(`Pre-1987`="grey20", `Historical and current management`="#69BE28", `Temporal management only`="#8AD2D8", `No spatial or temporal management`="#005594", 
                                 `Spatial and temporal management`="#53AF8B"), guide="none")+ 
  guides(fill = guide_legend(nrow = 6, label.position = "right",
                             override.aes = list(shape = c("square filled", "triangle filled", "diamond filled", "triangle down filled","circle","square filled", "triangle filled", "diamond filled", "triangle down filled"),
                                                 colour=c("#69BE28",  "#005594", "#8AD2D8", "#53AF8B", "white",  "white", "white", "white","white"),
                                                 fill = c("#69BE28",  "#005594", "#8AD2D8", "#53AF8B", "white",  "white", "white", "white","white"))))+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0, 0.1)+
  theme(#legend.key.size = unit(1, 'cm'), #change legend key size
    legend.key.height = unit(1, 'cm'), #change legend key height
    legend.key.width = unit(0.1, 'cm'), #change legend key width
    legend.title = element_text(size=12, face="italic"), #change legend title font size
    legend.text = element_text(size=12), #change legend text font size
    legend.text.align = 0,
    legend.title.align = 0) + 
  theme(axis.text=element_text(size=11))
  # ggplot2::annotate("text", x=2, y=0.499, label="(b) Sublegal sized", size = 3, fontface=2)
line.sublegal

line.legal <- ScenarioWholePop %>% 
  filter(Stage %in% c("Legal")) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Normal") & Year>1985, "Historical and current management", 
                                                             ifelse(Scenario %in% c("Temp Closure") & Year>1985, "Temporal management only", 
                                                                    ifelse(Scenario %in% c("Nothing"), "No spatial or temporal management", "Spatial and temporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  #filter(Status %in% c("NTZ")) %>% 
  mutate(PercChange = ifelse(Status %in% c("Fished"), ((Total-0.1839906)/0.1839906)*100, ((Total-0.5106572)/0.5106572)*100)) %>% 
  ggplot(.)+
  geom_point(aes(x=Year, y=Total, group=interaction(Status,Scenario), colour=ColourGroup, fill=ShapeGroup,  shape=ShapeGroup), size=2.5)+
  geom_line(aes(x=Year, y=Total, group=interaction(Status,Scenario), colour=ColourGroup))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dotted", color="grey20")+
  geom_vline(xintercept=10, linetype="dashed", colour="grey20")+
  geom_vline(xintercept=12.5,colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled", `NTZ.Temp Closure`="diamond filled", `NTZ.Temp Closure and NTZ`="triangle down filled" , `Pre-1987`="circle",
                               `Fished.Normal`="square filled", `Fished.Nothing`="triangle filled", `Fished.Temp Closure`="diamond filled", `Fished.Temp Closure and NTZ`="triangle down filled"), name="Area and scenario",
                     labels=mylabels3, guide="none")+
  scale_fill_manual(values= c(`Pre-1987`="grey20",  `NTZ.Normal`="#69BE28", `NTZ.Nothing`="#005594", `NTZ.Temp Closure`="#8AD2D8", `NTZ.Temp Closure and NTZ`="#53AF8B",
                              `Fished.Normal`="white", `Fished.Nothing`="white",`Fished.Temp Closure`="white", `Fished.Temp Closure and NTZ`="white"),
                    labels=mylabels2, name="General use\n(fished)")+
  scale_colour_manual(values = c(`Pre-1987`="grey20", `Historical and current management`="#69BE28", `Temporal management only`="#8AD2D8", `No spatial or temporal management`="#005594", 
                                 `Spatial and temporal management`="#53AF8B"), guide="none")+ 
  guides(fill = guide_legend(nrow = 6, label.position = "right",
                             override.aes = list(shape = c("square filled", "triangle filled", "diamond filled", "triangle down filled","circle","square filled", "triangle filled", "diamond filled", "triangle down filled"),
                                                 colour=c("#69BE28",  "#005594", "#8AD2D8", "#53AF8B", "white",  "white", "white", "white","white"),
                                                 fill = c("white",  "white", "white", "white","white","white","white", "white", "white"))))+
  xlab(NULL)+
  ylab(NULL)+
  ylim(0, 0.4)+
  theme(legend.key.height = unit(1, 'cm'), #change legend key height
    legend.key.width = unit(1, 'cm'), #change legend key width
    legend.title = element_text(size=12, face="italic"), #change legend title font size
    legend.text = element_text(size=12), #change legend text font size
    legend.text.align = 0,
    legend.title.align = 0) +
  theme(axis.text=element_text(size=11))
  #ggplot2::annotate("text", x=1.6, y=0.599, label="(a) Legal sized", size = 3, fontface=2)
line.legal


line.biglegal <- ScenarioWholePop %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Normal") & Year>1985, "Historical and current management", 
                                                             ifelse(Scenario %in% c("Temp Closure") & Year>1985, "Temporal management only", 
                                                                    ifelse(Scenario %in% c("Nothing"), "No spatial or temporal management", "Spatial and temporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  mutate(PercChange = ifelse(Status %in% c("Fished"), ((Total-0.12422065)/0.12422065)*100, ((Total-0.3447684)/0.3447684)*100)) %>% 
  #filter(Status %in% c("NTZ")) %>% 
  ggplot(.)+
  geom_point(aes(x=Year, y=Total, group=interaction(Status,Scenario), colour=ColourGroup, fill=ShapeGroup,  shape=ShapeGroup), size=2.5)+
  geom_line(aes(x=Year, y=Total, group=interaction(Status,Scenario), colour=ColourGroup))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dotted", color="grey20")+
  geom_vline(xintercept=10, linetype="dashed", colour="grey20")+
  geom_vline(xintercept=12.5,colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled", `NTZ.Temp Closure`="diamond filled", `NTZ.Temp Closure and NTZ`="triangle down filled" , `Pre-1987`="circle",
                               `Fished.Normal`="square filled", `Fished.Nothing`="triangle filled", `Fished.Temp Closure`="diamond filled", `Fished.Temp Closure and NTZ`="triangle down filled"), name="Area and scenario",
                     labels=mylabels3, guide="none")+
  scale_fill_manual(values= c(`Pre-1987`="grey20",  `NTZ.Normal`="#69BE28", `NTZ.Nothing`="#005594", `NTZ.Temp Closure`="#8AD2D8", `NTZ.Temp Closure and NTZ`="#53AF8B",
                              `Fished.Normal`="white", `Fished.Nothing`="white",`Fished.Temp Closure`="white", `Fished.Temp Closure and NTZ`="white"),
                    labels=mylabels2, name="General use\n(fished)")+
  scale_colour_manual(values = c(`Pre-1987`="grey20", `Historical and current management`="#69BE28", `Temporal management only`="#8AD2D8", `No spatial or temporal management`="#005594", 
                                 `Spatial and temporal management`="#53AF8B"), guide="none")+ 
  guides(fill = guide_legend(nrow = 6, label.position = "right",
                             override.aes = list(shape = c("square filled", "triangle filled", "diamond filled", "triangle down filled","circle","square filled", "triangle filled", "diamond filled", "triangle down filled"),
                                                 colour=c("#69BE28",  "#005594", "#8AD2D8", "#53AF8B", "white",  "white", "white", "white","white"),
                                                 fill = c("white",  "white", "white", "white","white","white","white", "white", "white"))))+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0, 0.25)+
  theme(axis.text=element_text(size=11))
  # ggplot2::annotate("text", x=1.95, y=0.399, label="(b) Large legal sized", size = 3, fontface=2)
line.biglegal

#### PUT PLOTS TOGETHER FOR PUBLISHING ####
x.label <- textGrob("Year", gp=gpar(fontsize=14))
y.label <- textGrob("No. Fish per"~km^2, gp=gpar(fontsize=14), rot=90)
top.legend <- gtable_filter(ggplotGrob(line.recruits), "guide-box")
NTZ.legend <- gtable_filter(ggplotGrob(line.sublegal), "guide-box")
Fished.legend <- gtable_filter(ggplotGrob(line.legal), "guide-box")

mat_layout <- rbind(c(1,1,1,1,NA,NA),
                    c(1,1,1,1,NA,NA),
                    c(1,1,1,1,NA,NA),
                    c(1,1,1,1,2,2),
                    c(1,1,1,1,3,4),
                    c(1,1,1,1,NA,NA),
                    c(1,1,1,1,NA,NA),
                    c(1,1,1,1,NA,NA))

# Size for copy plot is 1425 x 900

LinePlotsxGroup.L <-grid.arrange(arrangeGrob(line.legal + theme(legend.position="none"),
                                             line.biglegal + theme(legend.position="none"),
                                             ncol=2,
                                             left=y.label,
                                             bottom=x.label),
                                 arrangeGrob(top.legend),
                                 arrangeGrob(NTZ.legend),
                                 arrangeGrob(Fished.legend),
                                 widths=c(1,1,1,1,0.61,0.553),           
                                 layout_matrix = mat_layout)

LinePlotsxGroup.SL <-grid.arrange(arrangeGrob(line.recruits + theme(legend.position="none"),
                                              line.sublegal + theme(legend.position="none"),
                                              ncol=2,
                                              left=y.label,
                                              bottom=x.label),
                                  arrangeGrob(top.legend),
                                  arrangeGrob(NTZ.legend),
                                  arrangeGrob(Fished.legend),
                                  widths=c(1,1,1,1,0.61,0.553),           
                                  layout_matrix = mat_layout)

#### PLOT OF MODEL AREA ####
water <- water %>% 
  mutate(WHA = ifelse(ID %in% water_WHA$ID, "Y", "N")) %>% 
  mutate(WHA = ifelse(Fished_2017 %in% "N", "NTZ",WHA)) %>% 
  filter(!ID==387) # getting rid of weird ploygon that sticks out for some reason

water_points <- as.data.frame(st_coordinates(water_points)) %>% 
  mutate(left_join(., water$ID))

area_plot <- water %>% 
  ggplot(.)+
  geom_sf(aes(fill=WHA, colour=Fished_2017), lwd=0.25)+
  scale_fill_manual(values=c(`NTZ`="#33A02C", `Y`="#D6CF7D", `N`="#002D89"), name="Zone Type", labels=c("No-take zone", "World Heritage Area", "Outside World Heritage\nand marine park area"))+
  scale_colour_manual(values=c("grey20", "grey20"), guide="none")+
  theme_void() +
  theme(legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=12, face="bold"), #change legend title font size
        legend.text = element_text(size=12)) #change legend text font size
area_plot


#### SPATIAL PLOTS ####
setwd(sp_dir)
water <- readRDS(paste0(model.name, sep="_","water"))

pop.groups <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8)

setwd(pop_dir)

TotalPop_Normal <- readRDS("ningaloo_Rcpp_YearlyTotal_58") %>% 
  unlist()
TotalPop_Normal <- array(TotalPop_Normal, dim=c(1834,12,30)) %>% 
  rowSums(.[,,4:30], dim=2) 

TotalPop_Normal <- as.numeric(TotalPop_Normal[,12]) 
TotalPop_Normal <- TotalPop_Normal[c(as.numeric(water_WHA$ID))]
TotalPop_Normal <- TotalPop_Normal+1
TotalPop_Normal <- log(TotalPop_Normal)
TotalPop_Normal.df <- as.data.frame(TotalPop_Normal)

SpatialPlots_Normal <- spatial.plot.func(area=water_WHA, pop=TotalPop_Normal, pop.breaks=pop.groups, colours="YlGnBu")

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

#### CHECKING PLOTS ####
WholePop <- rbind(NoTakeAges, FishedAges) %>%  
  pivot_longer(cols=-c(Age, Status, Scenario), names_to="Year", values_to="Number") %>% 
  mutate(NumKM2 = ifelse(Status %in% c("Fished"), Number/AreaFished, Number/AreaNT)) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<=4, "Sublegal",
                               ifelse(Age>4 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA)))))


check <- WholePop %>% 
  filter(Stage %in% c("Recruit")) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre 1987", ifelse(Scenario %in% c("Normal") & Year>1985, "NTZs as normal", 
                                                             ifelse(Scenario %in% c("Temp Closure") & Year>1985, "Temporal Closure Only", 
                                                                    ifelse(Scenario %in% c("Nothing"), "None", "Temporal Closure and NTZs"))))) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  #mutate(PercChange = ifelse(Status %in% c("Fished"), ((Total-0.12931098)/0.12931098)*100, ((Total-0.3366771)/0.3366771)*100)) %>% 
  ggplot(.)+
  geom_point(aes(x=Year, y=Total, group=Status, colour=ShapeGroup)
                 #, colour=ColourGroup, fill=ColourGroup,  shape=ShapeGroup), size=2.5
                 )+
  geom_line(aes(x=Year, y=Total, group=Status, colour=ShapeGroup
                #, colour=ColourGroup
                ))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=10, colour="grey20")+
  geom_vline(xintercept=12.5, linetype="dashed", colour="grey20")+
  #geom_vline(xintercept=10.6, linetype="dotted", colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled" ,`Pre-1987`="circle",
                               `Fished.Normal`="square open", `Fished.Nothing`="triangle open"), name="Area and scenario",
                     labels=c("Pre-1987", "Normal Scenario (NTZ)","Nothing (NTZ)",
                              "Normal Scenario (Fished)", "Nothing (Fished)")
                     , guide="none")+
  ylab(NULL)+
  xlab(NULL)
# ylim(0,0.2)
#ggplot2::annotate("text", x=1.7, y=0.016, label="(c) Legal sized", size = 3, fontface=2)
check


for(Y in 13:59){
  
  year <- readRDS(paste0(model.name, sep="_", "S02_YearlyTotal_", numYear[Y]))
  saveRDS(year, file=paste0(model.name, sep="_", "S03_YearlyTotal_", numYear[Y]))
  
}

#### RECRUIT PLOTS FOR SHAUN ####
setwd(pop_dir)

# Whole Population 
Recruits <- array(0, dim=c(59,3))
Recruits <- as.data.frame(Recruits)
Recruits <- Recruits %>% 
  rename(Year = "V1") %>% 
  mutate(Year = seq(1960,2018,1)) %>% 
  rename(BH_Recruits = "V2") %>% 
  rename(Var_Recruits = "V3")

numYear <- seq(0,58,1)

for(Y in 1:59){
  
  year <- readRDS(paste0(model.name, sep="_", "Rcpp_BHRecruits_", numYear[Y])) %>% 
    unlist() %>% 
    as.numeric()

  Recruits[Y,2] <- sum(year)
  
}

for(Y in 1:59){
  
  year <- readRDS(paste0(model.name, sep="_", "Rcpp_YearlyTotal_", numYear[Y])) %>% 
    unlist()
  year <- array(year, dim=c(1834,12,30))
  year <- sum(year[,12,1])
  
  Recruits[Y,3] <- year
  
}

Recruit_plot <- Recruits %>% 
  ggplot()+
  geom_line(aes(x=Year, y=BH_Recruits))+
  geom_point(aes(x=Year, y=Var_Recruits))+
  ylab("No. Recruits")+
  theme_classic()+
  ylim(0, 2100)
Recruit_plot

Recruits <- Recruits %>% 
  mutate(Percent_change = (BH_Recruits-Var_Recruits)/((BH_Recruits+Var_Recruits)/2)*100)

Recruit_perc_change <- Recruits %>% 
  ggplot()+
  geom_line(aes(x=Year, y=Percent_change))+
  theme_classic()+
  ylab("% Change in Recruits")
Recruit_perc_change

Recruit_box <- Recruits %>% # Need to change the columns around so that you can make two box plots
  mutate(Age=" ") %>% 
  pivot_longer(cols=c("BH_Recruits", "Var_Recruits"), names_to="Model", values_to="Recruits") %>% 
  mutate(Model= ifelse(Model %in% c("BH_Recruits"), "Beverton-Holt", "Log-Normal Variation")) %>% 
  ggplot()+
  geom_boxplot(aes(y=Recruits, x=Model))+
  theme_classic()+
  ylab("No. Recruits")+
  ylim(0, 2100)
Recruit_box

Recruit_hist <- Recruits %>% 
  ggplot()+
  geom_histogram(aes(BH_Recruits), binwidth=25) + 
  theme_classic()+
  #scale_x_discrete(drop=FALSE)+
  ylab(NULL)+
  xlab("No. Recruits")
Recruit_hist


Recruits <- Recruits %>% 
  mutate(Percent_change_year = 0)

for (Y in 2:59){
  Recruits[Y,5] <- ((Recruits[Y,3]-Recruits[Y-1,3])/Recruits[Y-1,3])*100
}

Recruit_perc_change_year <- Recruits %>% 
  ggplot()+
  geom_line(aes(x=Year, y=Percent_change_year))+
  theme_classic()+
  ylab("% change in recruits relative to previous year")
Recruit_perc_change_year


#### PLOTTING SIMULATIONS ####
setwd(pop_dir)

total_pop <- readRDS("ningaloo_Total_Population_normal")
ages <- readRDS("ningaloo_Age_Distribution_normal")
catches <- readRDS("ningaloo_Yearly_Catch_normal")

## Population
total_pop <- as.data.frame(total_pop) %>% 
  mutate(Mod_Year = seq(1,59,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100])))

total_pop_plot <- total_pop %>% 
  ggplot() +
  geom_line(aes(x=Mod_Year, y=Mean_Pop))+
  geom_pointrange(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop))
total_pop_plot

## Catch

catches <- as.data.frame(catches) %>% 
  mutate(Mod_Year = seq(1,59,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate_at(vars(Sim1:Sim100), funs(./1000))%>% 
  mutate(Mean_Catch = rowMeans(.[,1:100])) %>% 
  mutate(SD_Catch = rowSds(as.matrix(.[1:100])))

catches_plot <- catches %>% 
  ggplot() +
  geom_line(aes(x=Mod_Year, y=Mean_Catch))+
  geom_pointrange(aes(x=Mod_Year, y=Mean_Catch, ymin=Mean_Catch-SD_Catch, ymax=Mean_Catch+SD_Catch))
catches_plot

## Age Distribution
age_dist <- NULL
for(L in 1:100){
  
  temp <- ages[,,L]
  
  temp2 <- as.data.frame(temp) %>% 
    mutate(Mod_Year = seq(1,59,1)) %>% 
    pivot_longer(cols=V1:V30, names_to="Age",values_to="Total") %>% 
    mutate(Age = as.numeric(str_replace(Age, "V", ""))) %>% 
    filter(Age>6) %>% 
    mutate(Number = Age * Total) %>% 
    group_by(Mod_Year) %>% 
    summarise(Mean_Age = sum(Number)/sum(Total)) 
    

  age_dist <- cbind(age_dist, temp2$Mean_Age)
  
}

age_dist <- as.data.frame(age_dist) %>% 
  mutate(Mean_Age = rowMeans(.[,1:100])) %>% 
  mutate(SD_Age = rowSds(as.matrix(.[1:100]))) %>% 
  mutate(Mod_Year = seq(1,59,1)) 

age_plot <- age_dist %>% 
  ggplot() +
  geom_line(aes(x=Mod_Year, y=Mean_Age))+
  geom_pointrange(aes(x=Mod_Year, y=Mean_Age, ymin=Mean_Age-SD_Age, ymax=Mean_Age+SD_Age))
age_plot


## Inside outside by age class
SP_Pop_NTZ <- readRDS("ningaloo_Sp_Population_NTZ")
SP_Pop_F <- readRDS("ningaloo_Sp_Population_F")

NoTakeAges <- array(0, dim=c(30,59))
NoTakeAges <- as.data.frame(NoTakeAges)
FishedAges <- NoTakeAges

numYears <- seq(0,59,1)

for(SIM in 1:1){
  
  for(YEAR in 1:59){
    
    Population.NT <- Sp_Pop_NTZ[,,YEAR] %>% 
      colSums(.)

    Population.F <- Sp_Pop_F[,,YEAR] %>% 
      colSums(.)
    
    NoTakeAges[,YEAR] <- Population.NT
    FishedAges[,YEAR] <- Population.F
    
    
      NoTakeAges <- NoTakeAges %>% 
        mutate(Status = "NTZ") %>% 
        mutate(Age = seq(1:30)) %>% 
        mutate(Scenario = "Normal") %>% 
        mutate(Scenario = as.factor(Scenario))
      
      FishedAges <- FishedAges %>% 
        mutate(Status = "Fished") %>% 
        mutate(Age = seq(1:30)) %>% 
        mutate(Scenario = "Normal") %>% 
        mutate(Scenario = as.factor(Scenario))
      
      WholePop <- rbind(NoTakeAges, FishedAges) %>%
        pivot_longer(cols=-c(Age, Status, Scenario), names_to="Year", values_to="Number") %>%
        mutate(Year = rep(1960:2018, length.out=nrow(.))) %>% 
        mutate(Year = as.factor(Year)) %>%
        mutate(Year = fct_reorder(Year, as.numeric(Year))) %>%
        mutate(Status = as.factor(Status))
      
      NoRecruits <- WholePop %>%
        filter(Age!=1)
    
  }
}

WholePop <- WholePop %>%
  mutate(NumKM2 = ifelse(Status %in% c("Fished"), Number/AreaFished, Number/AreaNT)) %>%
  mutate(numYear = as.numeric(Year)) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<=3, "Sublegal",
                               ifelse(Age>3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA)))))

check <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(ColourGroup = ifelse(numYear<=24, "Pre 1987", "NTZs as normal")) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(numYear>24, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  #filter(Status %in% c("Fished")) %>% 
  #mutate(PercChange = ifelse(Status %in% c("Fished"), ((Total-0.12931098)/0.12931098)*100, ((Total-0.3366771)/0.3366771)*100)) %>% 
  ggplot(.)+
  geom_point(aes(x=Year, y=Total, group=Status, colour=ShapeGroup)
             #, colour=ColourGroup, fill=ColourGroup,  shape=ShapeGroup), size=2.5
  )+
  geom_line(aes(x=Year, y=Total, group=Status, colour=ShapeGroup
                #, colour=ColourGroup
  ))+
  theme_classic()+
  # geom_vline(xintercept=6.6, linetype="dashed", color="grey20")+
  # geom_vline(xintercept=10, colour="grey20")+
  # geom_vline(xintercept=12.5, linetype="dashed", colour="grey20")+
  #geom_vline(xintercept=10.6, linetype="dotted", colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled" ,`Pre-1987`="circle",
                               `Fished.Normal`="square open", `Fished.Nothing`="triangle open"), name="Area and scenario",
                     labels=c("Pre-1987", "Normal Scenario (NTZ)","Nothing (NTZ)",
                              "Normal Scenario (Fished)", "Nothing (Fished)")
                     , guide="none")+
  ylab(NULL)+
  xlab(NULL)
check



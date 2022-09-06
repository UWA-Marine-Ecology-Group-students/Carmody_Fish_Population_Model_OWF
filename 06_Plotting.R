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


# Normal
# Sim 1 Nothing
# Sim 2 NTZs and Temporal Closure
# Sim 3 Just temporal closure, no sanctuary zones 


#### READ IN DATA ####
setwd(sg_dir)
NoTake <- readRDS("NoTakeList")


setwd(sp_dir)
water <- readRDS("water")
bathy <- raster("ga_bathy_ningaloocrop.tif")

#* Create list of cells to restrict plots to shallow water (<30m)

water_bathy <- extract(bathy, water, fun=mean, df=TRUE)

water_bathy <- water_bathy %>% 
  mutate(ID = as.factor(ID))

water <- water %>% 
  left_join(., water_bathy, by="ID") %>% 
  rename(bathy = "ga_bathy_ningaloocrop") %>% 
  filter(bathy >= c(-30)) %>% 
  filter(!is.na(bathy))

shallow_cells_NTZ <- water %>% 
  filter(Fished %in% c("N")) %>% 
  st_drop_geometry(.)

shallow_cells_F <- water %>% 
  filter(Fished %in% c("Y")) %>% 
  st_drop_geometry(.)

shallow_NTZ_ID <- shallow_cells_NTZ$ID
shallow_F_ID <- shallow_cells_F$ID

#* Normal - RE-DONE ####
setwd(pop_dir)

# Whole Population 
TotalPop <- array(0, dim=c(59,2))
TotalPop <- as.data.frame(TotalPop)
TotalPop <- TotalPop %>% 
  rename(Year="V1") %>% 
  mutate(Year=seq(1960,2018,1)) %>% 
  rename(Tot.Pop="V2")

numYear <- seq(1,59,1)

for(Y in 1:59){
  
  year <- readRDS(paste0("YearlyTotal_", numYear[Y]))
  year <- rowSums(year[,,1:30], dim=2)
  year <- sum(year[,12])
  
  TotalPop[Y,2] <- year
  
}

# Separated by Fished and No-Take
NoTakeAges <- array(0, dim=c(30,13))
NoTakeAges <- as.data.frame(NoTakeAges)
FishedAges <- NoTakeAges


numYears <- seq(0,59,5)
numYears[1] <- 1
numYears[13] <- 59

for(YEAR in 1:13){
  
  Population <-  readRDS(paste0("YearlyTotal_", numYears[YEAR]))
  
  Population.NT <- Population[c(NoTake[[3]]),12, ] %>% 
    colSums(.)
  Population.F <- Population[-c(NoTake[[3]]),12, ] %>% 
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
  rename(Tot.Pop="V2")

numYear <- seq(1,59,1)

for(Y in 1:59){
  
  year <- readRDS(paste0("S01_YearlyTotal_", numYear[Y]))
  year <- rowSums(year[,,1:30], dim=2)
  year <- sum(year[,12])
  
  s1_TotalPop[Y,2] <- year
  
}


# Separated by Fished and No-Take
s1_NoTakeAges <- array(0, dim=c(30,13))
s1_NoTakeAges <- as.data.frame(s1_NoTakeAges)
s1_FishedAges <- s1_NoTakeAges


numYears <- seq(0,59,5)
numYears[1] <- 1
numYears[13] <- 59

for(YEAR in 1:13){
  
  Population <-  readRDS(paste0("S01_YearlyTotal_", numYears[YEAR]))
  
  Population.NT <- Population[c(NoTake[[3]]),12, ] %>% 
    colSums(.)
  Population.F <- Population[-c(NoTake[[3]]),12, ] %>% 
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
  rename(Tot.Pop="V2")

numYear <- seq(1,59,1)

for(Y in 1:59){
  
  year <- readRDS(paste0("S02_YearlyTotal_", numYear[Y]))
  year <- rowSums(year[,,1:30], dim=2)
  year <- sum(year[,12])
  
  s2_TotalPop[Y,2] <- year
  
}


# Separated by Fished and No-Take
s2_NoTakeAges <- array(0, dim=c(30,13))
s2_NoTakeAges <- as.data.frame(s2_NoTakeAges)
s2_FishedAges <- s2_NoTakeAges


numYears <- seq(0,59,5)
numYears[1] <- 1
numYears[13] <- 59

for(YEAR in 1:13){
  
  Population <-  readRDS(paste0("S02_YearlyTotal_", numYears[YEAR]))
  
  Population.NT <- Population[c(NoTake[[3]]),12, ] %>% 
    colSums(.)
  Population.F <- Population[-c(NoTake[[3]]),12, ] %>% 
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
  rename(Tot.Pop="V2")

numYear <- seq(1,59,1)

for(Y in 1:59){
  
  year <- readRDS(paste0("S03_YearlyTotal_", numYear[Y]))
  year <- rowSums(year[,,1:30], dim=2)
  year <- sum(year[,12])
  
  s3_TotalPop[Y,2] <- year
  
}


# Separated by Fished and No-Take
s3_NoTakeAges <- array(0, dim=c(30,13))
s3_NoTakeAges <- as.data.frame(s3_NoTakeAges)
s3_FishedAges <- s3_NoTakeAges


numYears <- seq(0,59,5)
numYears[1] <- 1
numYears[13] <- 59

for(YEAR in 1:13){
  
  Population <-  readRDS(paste0("S03_YearlyTotal_", numYears[YEAR]))
  
  Population.NT <- Population[c(NoTake[[3]]),12, ] %>% 
    colSums(.)
  Population.F <- Population[-c(NoTake[[3]]),12, ] %>% 
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


#### CALCULATE TOTAL AREA OF FISHED AND NO-TAKE ####

AreaFished <- water %>% 
  mutate(cell_area = st_area(Spatial)) %>% 
  mutate(cell_area = as.numeric(cell_area)) %>% 
  mutate(Fished=as.factor(Fished_17)) %>% 
  filter(Fished=="Y") 

AreaFished <- (sum(AreaFished$cell_area))/100000

AreaNT <- water %>% 
  mutate(cell_area = st_area(Spatial)) %>% 
  mutate(cell_area = as.numeric(cell_area)) %>% 
  mutate(Fished=as.factor(Fished_17)) %>% 
  filter(Fished=="N") 

AreaNT <- (sum(AreaNT$cell_area))/100000
#### TIME SERIES PLOT WITH ALL SCENARIOS ####

TotalPop <- TotalPop %>% 
  rename(Normal = "Tot.Pop",
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
  scale_colour_manual("Scenario",values=c( "gray20",  "#69BE28", "#005594", "#8AD2D8", "#53AF8B"), labels=c("Pre 1987", "NTZs in place as normal", 
                                                                                                            "No NTZs or Temporal Closure" ,"Temporal Closure Only", 
                                                                                                            "Temporal Closure and NTZs in place"))+
  xlab("Year")+
  ylab("Total Population")+
  geom_vline(xintercept=1987, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2008, linetype="dotted", colour="grey20")+
  theme_classic()

  
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


line.recruits.NTZ <- ScenarioWholePop %>% 
  filter(Age==1) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre 1987", ifelse(Scenario %in% c("Normal") & Year>1985, "NTZs as normal", 
                                                             ifelse(Scenario %in% c("Temp Closure") & Year>1985, "Temporal Closure Only", 
                                                                    ifelse(Scenario %in% c("Nothing"), "None", "Temporal Closure and NTZs"))))) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  filter(Status %in% c("NTZ")) %>% 
  ggplot(., aes(x=Year, y=NumKM2, group=Scenario, colour=ColourGroup))+
  geom_point(size=2.5, aes(fill=ShapeGroup,  shape=ShapeGroup))+
  geom_line(aes(colour=ColourGroup))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=10, colour="grey20")+
  geom_vline(xintercept=12.5, linetype="dashed", colour="grey20")+
  #geom_vline(xintercept=10.6, linetype="dotted", colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled", `NTZ.Temp Closure`="diamond filled", `NTZ.Temp Closure and NTZ`="triangle down filled" ,
                               `Fished.Normal`="square filled", `Fished.Nothing`="triangle filled", `Fished.Temp Closure`="diamond filled", `Fished.Temp Closure and NTZ`="triangle down filled"), name="Area and scenario",
                     labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal CLosure and NTZs (NTZ)",
                              "Normal Scenario (Fished)", "No Management (Fished)",
                              "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)"),
                     guide="none")+
  scale_fill_manual(values= c(`Pre-1987`="grey20", `NTZ.Normal`="#69BE28", `NTZ.Nothing`="#005594", `NTZ.Temp Closure`="#8AD2D8", `NTZ.Temp Closure and NTZ`="#53AF8B",
                              `Fished.Normal`="white", `Fished.Nothing`="white",`Fished.Temp Closure`="white", `Fished.Temp Closure and NTZ`="white"),
                    labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                             "Normal Scenario (Fished)", "No Management (Fished)",
                             "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)"), name="Area and Scenario")+
  scale_colour_manual(values = c("#005594",  "#69BE28", "grey20", "#53AF8B", "#8AD2D8"), guide="none")+ 
  guides(fill = guide_legend(override.aes = list(shape = c("circle", "square filled", "triangle filled", "diamond filled", "triangle down filled", "square filled", "triangle filled", "diamond filled", "triangle down filled"),
                                                 colour=c("grey20", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B"))))+
  ylab(NULL)+
  xlab(NULL)+
  ggplot2::annotate("text", x=1.5, y=1.7, label="(a) Recruits", size = 3, fontface=2)

line.recruits.Fished <- ScenarioWholePop %>% 
  filter(Age==1) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre 1987", ifelse(Scenario %in% c("Normal") & Year>1985, "NTZs as normal", 
                                                             ifelse(Scenario %in% c("Temp Closure") & Year>1985, "Temporal Closure Only", 
                                                                    ifelse(Scenario %in% c("Nothing"), "None", "Temporal Closure and NTZs"))))) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  filter(Status %in% c("Fished")) %>% 
  ggplot(., aes(x=Year, y=NumKM2, group=Scenario, colour=ColourGroup))+
  geom_point(size=2.5, aes(fill=ShapeGroup,  shape=ShapeGroup))+
  geom_line(aes(colour=ColourGroup))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=10, colour="grey20")+
  geom_vline(xintercept=12.5, linetype="dashed", colour="grey20")+
  #geom_vline(xintercept=10.6, linetype="dotted", colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled", `NTZ.Temp Closure`="diamond filled", `NTZ.Temp Closure and NTZ`="triangle down filled" ,
                               `Fished.Normal`="square filled", `Fished.Nothing`="triangle filled", `Fished.Temp Closure`="diamond filled", `Fished.Temp Closure and NTZ`="triangle down filled"), name="Area and scenario",
                     labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal CLosure and NTZs (NTZ)",
                              "Normal Scenario (Fished)", "No Management (Fished)",
                              "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)"),
                     guide="none")+
  scale_fill_manual(values= c(`Pre-1987`="grey20", `NTZ.Normal`="#69BE28", `NTZ.Nothing`="#005594", `NTZ.Temp Closure`="#8AD2D8", `NTZ.Temp Closure and NTZ`="#53AF8B",
                              `Fished.Normal`="white", `Fished.Nothing`="white",`Fished.Temp Closure`="white", `Fished.Temp Closure and NTZ`="white"),
                    labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                             "Normal Scenario (Fished)", "No Management (Fished)",
                             "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)"), name="Area and Scenario")+
  scale_colour_manual(values = c("#005594",  "#69BE28", "grey20", "#53AF8B", "#8AD2D8"), guide="none")+ 
  guides(fill = guide_legend(override.aes = list(shape = c("circle", "square filled", "triangle filled", "diamond filled", "triangle down filled", "square filled", "triangle filled", "diamond filled", "triangle down filled"),
                                                 colour=c("grey20", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B"))))+
  ylab(NULL)+
  xlab(NULL)+
  ggplot2::annotate("text", x=1.5, y=0.13, label="(a) Recruits", size = 3, fontface=2)
line.recruits.Fished

line.sublegal.NTZ <- ScenarioWholePop %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre 1987", ifelse(Scenario %in% c("Normal") & Year>1985, "NTZs as normal", 
                                                                ifelse(Scenario %in% c("Temp Closure") & Year>1985, "Temporal Closure Only", 
                                                                       ifelse(Scenario %in% c("Nothing"), "None", "Temporal Closure and NTZs"))))) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  filter(Status %in% c("NTZ")) %>% 
  ggplot(.)+
  geom_point(aes(x=Year, y=Total, group=Scenario, colour=ShapeGroup, fill=ShapeGroup,  shape=ShapeGroup), size=2.5,)+
  geom_line(aes(x=Year, y=Total, group=Scenario, colour=ColourGroup))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=10, colour="grey20")+
  geom_vline(xintercept=12.5, linetype="dashed", colour="grey20")+
  #geom_vline(xintercept=10.6, linetype="dotted", colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled", `NTZ.Temp Closure`="diamond filled", `NTZ.Temp Closure and NTZ`="triangle down filled" , `Pre-1987`="circle",
                               `Fished.Normal`="square filled", `Fished.Nothing`="triangle filled", `Fished.Temp Closure`="diamond filled", `Fished.Temp Closure and NTZ`="triangle down filled"), name="Area and scenario",
                     labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                              "Normal Scenario (Fished)", "No Management (Fished)",
                              "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)")
                     , guide="none")+
  scale_fill_manual(values= c(`Pre 1987`="grey20", `NTZ.Normal`="#69BE28", `NTZ.Nothing`="#005594", `NTZ.Temp Closure`="#8AD2D8", `NTZ.Temp Closure and NTZ`="#53AF8B",
                              `Fished.Normal`="white", `Fished.Nothing`="white",`Fished.Temp Closure`="white", `Fished.Temp Closure and NTZ`="white"),
                    labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                             "Normal Scenario (Fished)", "No Management (Fished)",
                             "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)"), name="Area and Scenario")+
  scale_colour_manual(values = c("#005594", "#69BE28", "#005594", "#8AD2D8", "#53AF8B", "#69BE28", "grey 20", "grey 20", "#53AF8B", "#8AD2D8"), guide="none")+ 
  guides(fill = guide_legend(override.aes = list(shape = c("circle", "square filled", "triangle filled", "diamond filled", "triangle down filled", "square filled", "triangle filled", "diamond filled", "triangle down filled"),
                                                 colour=c("grey20", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B"))))+
  ylab(NULL)+
  xlab(NULL)+
  ggplot2::annotate("text", x=1.9, y=0.24, label="(b) Sub-legal sized", size = 3, fontface=2)
line.sublegal.NTZ

line.sublegal.Fished <- ScenarioWholePop %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre 1987", ifelse(Scenario %in% c("Normal") & Year>1985, "NTZs as normal", 
                                                             ifelse(Scenario %in% c("Temp Closure") & Year>1985, "Temporal Closure Only", 
                                                                    ifelse(Scenario %in% c("Nothing"), "None", "Temporal Closure and NTZs"))))) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  filter(Status %in% c("Fished")) %>% 
  ggplot(.)+
  geom_point(aes(x=Year, y=Total, group=Scenario, colour=ShapeGroup, fill=ShapeGroup,  shape=ShapeGroup), size=2.5,)+
  geom_line(aes(x=Year, y=Total, group=Scenario, colour=ColourGroup))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=10, colour="grey20")+
  geom_vline(xintercept=12.5, linetype="dashed", colour="grey20")+
  #geom_vline(xintercept=10.6, linetype="dotted", colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled", `NTZ.Temp Closure`="diamond filled", `NTZ.Temp Closure and NTZ`="triangle down filled" , `Pre-1987`="circle",
                               `Fished.Normal`="square filled", `Fished.Nothing`="triangle filled", `Fished.Temp Closure`="diamond filled", `Fished.Temp Closure and NTZ`="triangle down filled"), name="Area and scenario",
                     labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                              "Normal Scenario (Fished)", "No Management (Fished)",
                              "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)")
                     , guide="none")+
  scale_fill_manual(values= c(`Pre-1987`="grey20", `NTZ.Normal`="#69BE28", `NTZ.Nothing`="#005594", `NTZ.Temp Closure`="#8AD2D8", `NTZ.Temp Closure and NTZ`="#53AF8B",
                              `Fished.Normal`="white", `Fished.Nothing`="white",`Fished.Temp Closure`="white", `Fished.Temp Closure and NTZ`="white"),
                    labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                             "Normal Scenario (Fished)", "No Management (Fished)",
                             "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)"), name="Area and Scenario")+
  scale_colour_manual(values = c("#69BE28", "#005594", "#8AD2D8", "#53AF8B", "#005594", "#69BE28", "grey 20", "grey 20", "#53AF8B", "#8AD2D8"), guide="none")+ 
  guides(fill = guide_legend(override.aes = list(shape = c("circle", "square filled", "triangle filled", "diamond filled", "triangle down filled", "square filled", "triangle filled", "diamond filled", "triangle down filled"),
                                                 colour=c("grey20", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B"))))+
  ylab(NULL)+
  xlab(NULL)+
  ggplot2::annotate("text", x=1.9, y=0.0148, label="(b) Sub-legal sized", size = 3, fontface=2)
line.sublegal.Fished

line.legal.NTZ <- ScenarioWholePop %>% 
  filter(Stage %in% c("Legal")) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre 1987", ifelse(Scenario %in% c("Normal") & Year>1985, "NTZs as normal", 
                                                             ifelse(Scenario %in% c("Temp Closure") & Year>1985, "Temporal Closure Only", 
                                                                    ifelse(Scenario %in% c("Nothing"), "None", "Temporal Closure and NTZs"))))) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  filter(Status %in% c("NTZ")) %>% 
  ggplot(.)+
  geom_point(aes(x=Year, y=Total, group=Scenario, colour=ShapeGroup, fill=ShapeGroup,  shape=ShapeGroup), size=2.5)+
  geom_line(aes(x=Year, y=Total, group=Scenario, colour=ColourGroup))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=10, colour="grey20")+
  geom_vline(xintercept=12.5, linetype="dashed", colour="grey20")+
  #geom_vline(xintercept=10.6, linetype="dotted", colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled", `NTZ.Temp Closure`="diamond filled", `NTZ.Temp Closure and NTZ`="triangle down filled" , `Pre-1987`="circle",
                               `Fished.Normal`="square filled", `Fished.Nothing`="triangle filled", `Fished.Temp Closure`="diamond filled", `Fished.Temp Closure and NTZ`="triangle down filled"), name="Area and scenario",
                     labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                              "Normal Scenario (Fished)", "No Management (Fished)",
                              "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)")
                     , guide="none")+
  scale_fill_manual(values= c(`Pre-1987`="grey20", `NTZ.Normal`="#69BE28", `NTZ.Nothing`="#005594", `NTZ.Temp Closure`="#8AD2D8", `NTZ.Temp Closure and NTZ`="#53AF8B",
                              `Fished.Normal`="white", `Fished.Nothing`="white",`Fished.Temp Closure`="white", `Fished.Temp Closure and NTZ`="white"),
                    labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                             "Normal Scenario (Fished)", "No Management (Fished)",
                             "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)"), name="Area and Scenario")+
  scale_colour_manual(values = c("#005594", "#69BE28", "#005594", "#8AD2D8", "#53AF8B", "#69BE28", "grey 20", "grey 20", "#53AF8B", "#8AD2D8"), guide="none")+ 
  guides(fill = guide_legend(override.aes = list(shape = c("circle", "square filled", "triangle filled", "diamond filled", "triangle down filled", "square filled", "triangle filled", "diamond filled", "triangle down filled"),
                                                 colour=c("grey20", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B"))))+
  ylab(NULL)+
  xlab(NULL)
  #ggplot2::annotate("text", x=1.7, y=0.3, label="(c) Legal sized", size = 3, fontface=2)
line.legal.NTZ

line.legal.Fished <- ScenarioWholePop %>% 
  filter(Stage %in% c("Legal")) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre 1987", ifelse(Scenario %in% c("Normal") & Year>1985, "NTZs as normal", 
                                                             ifelse(Scenario %in% c("Temp Closure") & Year>1985, "Temporal Closure Only", 
                                                                    ifelse(Scenario %in% c("Nothing"), "None", "Temporal Closure and NTZs"))))) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  filter(Status %in% c("Fished")) %>% 
  ggplot(.)+
  geom_point(aes(x=Year, y=Total, group=Scenario, colour=ShapeGroup, fill=ShapeGroup,  shape=ShapeGroup), size=2.5)+
  geom_line(aes(x=Year, y=Total, group=Scenario, colour=ColourGroup))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=10, colour="grey20")+
  geom_vline(xintercept=12.5, linetype="dashed", colour="grey20")+
  #geom_vline(xintercept=10.6, linetype="dotted", colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled", `NTZ.Temp Closure`="diamond filled", `NTZ.Temp Closure and NTZ`="triangle down filled" , `Pre-1987`="circle",
                               `Fished.Normal`="square filled", `Fished.Nothing`="triangle filled", `Fished.Temp Closure`="diamond filled", `Fished.Temp Closure and NTZ`="triangle down filled"), name="Area and scenario",
                     labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                              "Normal Scenario (Fished)", "No Management (Fished)",
                              "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)")
                     , guide="none")+
  scale_fill_manual(values= c(`Pre-1987`="grey20", `NTZ.Normal`="#69BE28", `NTZ.Nothing`="#005594", `NTZ.Temp Closure`="#8AD2D8", `NTZ.Temp Closure and NTZ`="#53AF8B",
                              `Fished.Normal`="white", `Fished.Nothing`="white",`Fished.Temp Closure`="white", `Fished.Temp Closure and NTZ`="white"),
                    labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                             "Normal Scenario (Fished)", "No Management (Fished)",
                             "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)"), name="Area and Scenario")+
  scale_colour_manual(values = c("#69BE28", "#005594", "#8AD2D8", "#53AF8B", "#005594", "#69BE28", "grey 20", "grey 20", "#53AF8B", "#8AD2D8"), guide="none")+ 
  guides(fill = guide_legend(override.aes = list(shape = c("circle", "square filled", "triangle filled", "diamond filled", "triangle down filled", "square filled", "triangle filled", "diamond filled", "triangle down filled"),
                                                 colour=c("grey20", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B"))))+
  ylab(NULL)+
  xlab(NULL)
  #ggplot2::annotate("text", x=1.7, y=0.016, label="(c) Legal sized", size = 3, fontface=2)
line.legal.Fished

line.biglegal.NTZ <- ScenarioWholePop %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "grey20", ifelse(Scenario %in% c("Normal") & Year>1985, "#69BE28", 
                                                             ifelse(Scenario %in% c("Temp Closure") & Year>1985, "#8AD2D8", 
                                                                    ifelse(Scenario %in% c("Nothing"), "#005594", "#53AF8B"))))) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  filter(Status %in% c("NTZ")) %>% 
  ggplot(.)+
  geom_point(aes(x=Year, y=Total, group=interaction(Status,Scenario), colour=ShapeGroup, fill=ShapeGroup,  shape=ShapeGroup), size=2.5,)+
  geom_line(aes(x=Year, y=Total, group=interaction(Status,Scenario), colour=ColourGroup))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=10, colour="grey20")+
  geom_vline(xintercept=12.5, linetype="dashed", colour="grey20")+
  #geom_vline(xintercept=10.6, linetype="dotted", colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled", `NTZ.Temp Closure`="diamond filled", `NTZ.Temp Closure and NTZ`="triangle down filled" , `Pre-1987`="circle",
                               `Fished.Normal`="square filled", `Fished.Nothing`="triangle filled", `Fished.Temp Closure`="diamond filled", `Fished.Temp Closure and NTZ`="triangle down filled"), name="Area and scenario",
                     labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                              "Normal Scenario (Fished)", "No Management (Fished)",
                              "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)")
                     , guide="none")+
  scale_fill_manual(values= c(`Pre-1987`="grey20", `NTZ.Normal`="#69BE28", `NTZ.Nothing`="#005594", `NTZ.Temp Closure`="#8AD2D8", `NTZ.Temp Closure and NTZ`="#53AF8B",
                              `Fished.Normal`="white", `Fished.Nothing`="white",`Fished.Temp Closure`="white", `Fished.Temp Closure and NTZ`="white"),
                    labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                             "Normal Scenario (Fished)", "No Management (Fished)",
                             "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)"), name="Area and Scenario")+
  scale_colour_manual(values = c("#005594", "#53AF8B", "#69BE28", "#8AD2D8", "grey 20", "#69BE28", "#005594", "#8AD2D8", "#53AF8B", "grey 20"), guide="none")+ 
  guides(fill = guide_legend(override.aes = list(shape = c("circle", "square filled", "triangle filled", "diamond filled", "triangle down filled", "square filled", "triangle filled", "diamond filled", "triangle down filled"),
                                                 colour=c("grey20", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B"))))+
  ylab(NULL)+
  xlab(NULL)+
  ggplot2::annotate("text", x=2, y=0.2, label="(d) Large legal sized", size = 3, fontface=2)
line.biglegal.NTZ

line.biglegal.Fished <- ScenarioWholePop %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "grey20", ifelse(Scenario %in% c("Normal") & Year>1985, "#69BE28", 
                                                           ifelse(Scenario %in% c("Temp Closure") & Year>1985, "#8AD2D8", 
                                                                  ifelse(Scenario %in% c("Nothing"), "#005594", "#53AF8B"))))) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ShapeGroup = ifelse(Year>1985, paste(Status, Scenario, sep="."), "Pre-1987")) %>% 
  mutate(ShapeGroup = as.factor(ShapeGroup)) %>% 
  filter(Status %in% c("Fished")) %>% 
  ggplot(.)+
  geom_point(aes(x=Year, y=Total, group=interaction(Status,Scenario), colour=ShapeGroup, fill=ShapeGroup,  shape=ShapeGroup), size=2.5,)+
  geom_line(aes(x=Year, y=Total, group=interaction(Status,Scenario), colour=ColourGroup))+
  theme_classic()+
  geom_vline(xintercept=6.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=10, colour="grey20")+
  geom_vline(xintercept=12.5, linetype="dashed", colour="grey20")+
 # geom_vline(xintercept=10.6, linetype="dotted", colour="grey20")+
  scale_shape_manual(values= c(`Pre-1987`="circle",`NTZ.Normal`="square filled", `NTZ.Nothing`="triangle filled", `NTZ.Temp Closure`="diamond filled", `NTZ.Temp Closure and NTZ`="triangle down filled" , `Pre-1987`="circle",
                               `Fished.Normal`="square filled", `Fished.Nothing`="triangle filled", `Fished.Temp Closure`="diamond filled", `Fished.Temp Closure and NTZ`="triangle down filled"), name="Area and scenario",
                     labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                              "Normal Scenario (Fished)", "No Management (Fished)",
                              "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)")
                     , guide="none")+
  scale_fill_manual(values= c(`Pre-1987`="grey20", `NTZ.Normal`="#69BE28", `NTZ.Nothing`="#005594", `NTZ.Temp Closure`="#8AD2D8", `NTZ.Temp Closure and NTZ`="#53AF8B",
                              `Fished.Normal`="white", `Fished.Nothing`="white",`Fished.Temp Closure`="white", `Fished.Temp Closure and NTZ`="white"),
                    labels=c("Pre-1987", "Normal Scenario (NTZ)", "No Management (NTZ)", "Temporal Closure (NTZ)", "Temporal Closure and NTZs (NTZ)",
                             "Normal Scenario (Fished)", "No Management (Fished)",
                             "Temporal Closure (Fished)" ,"Temporal Closure and NTZs (Fished)"), name="Area and Scenario")+
  scale_colour_manual(values = c("#005594", "#53AF8B", "#69BE28", "#8AD2D8", "#69BE28", "#005594", "#8AD2D8", "#53AF8B", "grey 20", "grey 20"), guide="none")+ 
  guides(fill = guide_legend(override.aes = list(shape = c("circle", "square filled", "triangle filled", "diamond filled", "triangle down filled", "square filled", "triangle filled", "diamond filled", "triangle down filled"),
                                                 colour=c("grey20", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B", "#69BE28",  "#005594", "#8AD2D8", "#53AF8B"))))+
  ylab(NULL)+
  xlab(NULL)+
  ggplot2::annotate("text", x=2, y=0.012, label="(d) Large legal sized", size = 3, fontface=2)
line.biglegal.Fished

#### PUT PLOTS TOGETHER FOR PUBLISHING ####
x.label <- textGrob("Year", gp=gpar(fontsize=14))
y.label <- textGrob("No. Fish per"~km^2, gp=gpar(fontsize=14), rot=90)
legend <- gtable_filter(ggplotGrob(line.biglegal), "guide-box")

LinePlotsxGroup.NTZ <-grid.arrange(arrangeGrob(line.recruits.NTZ + theme(legend.position="none"),
                                           line.sublegal.NTZ + theme(legend.position="none"),
                                           line.legal.NTZ + theme(legend.position="none"),
                                           line.biglegal.NTZ + theme(legend.position="none"),
                                           nrow = 2,
                                           left = y.label,
                                           bottom = x.label), 
                               legend, 
                               widths=unit.c(unit(1, "npc") - legend$width, legend$width), 
                               nrow=1)

LinePlotsxGroup.F <-grid.arrange(arrangeGrob(line.recruits.Fished + theme(legend.position="none"),
                                               line.sublegal.Fished + theme(legend.position="none"),
                                               line.legal.Fished + theme(legend.position="none"),
                                               line.biglegal.Fished + theme(legend.position="none"),
                                               nrow = 2,
                                               left = y.label,
                                               bottom = x.label), 
                                   legend, 
                                   widths=unit.c(unit(1, "npc") - legend$width, legend$width), 
                                   nrow=1)

LinePlotsxSubLegal <-grid.arrange(arrangeGrob(line.recruits.NTZ + theme(legend.position="none"),
                                              line.recruits.Fished + theme(legend.position="none"),
                                              line.sublegal.NTZ + theme(legend.position="none"),
                                              line.sublegal.Fished + theme(legend.position="none"),
                                              nrow = 2,
                                              left = y.label,
                                              bottom = x.label), 
                                  legend, 
                                  widths=unit.c(unit(1, "npc") - legend$width, legend$width), 
                                  nrow=1)

LinePlotsxLegal <-grid.arrange(arrangeGrob(line.legal.NTZ + theme(legend.position="none"),
                                              line.legal.Fished + theme(legend.position="none"),
                                              line.biglegal.NTZ + theme(legend.position="none"),
                                              line.biglegal.Fished + theme(legend.position="none"),
                                              nrow = 2,
                                              left = y.label,
                                              bottom = x.label), 
                                  legend, 
                                  widths=unit.c(unit(1, "npc") - legend$width, legend$width), 
                                  nrow=1)

#### KOBE STYLE PLOT ####
Years <- seq(5, 55, 5)
Years[12] <- 59

TotMatBio <- matrix(0, ncol=12, nrow=30)

for(Y in 1:12){
  
  Population <-  get(paste0("Year", Years[Y]))
  
  for (A in 1:30){
    adults <- Population[ ,10, ] %>% 
      colSums(.) # Gives us just females because they are the limiting factor for reproduction
    adults <- adults * 0.5
    
    MatBio<- lapply(1:dim(Population)[3], function(Age){
      SB <- adults[Age] * maturity[Age,1] * weight[(Age*12)+1] #Gives us spawning biomass in each age group at the end of the year, hence the x 12+1 as it starts at 1 not zero
      TotMatBio <- sum(SB) #Gives us total mature spawning biomass
    })
    MatBio <- do.call(rbind, MatBio)
    TotMatBio[ ,Y] <- MatBio
  }
}

TotMatBio <- as.data.frame(TotMatBio) %>% 
  `colnames<-`(c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018")) %>% 
  mutate(Age = seq(1, 30, 1)) %>% 
  pivot_longer(cols=-c(Age), names_to="Year", values_to="Number") %>% 
  mutate(Year = as.factor(Year)) %>% 
  group_by(Year) %>% 
  summarise(TotalBio=sum(Number)) %>% 
  mutate(TotalBio = TotalBio/1000)

YearlyFishing <- fishing %>% 
  colSums(fishing) %>% 
  colSums(.)

YearlyFishing <- YearlyFishing[Years, ]
YearlyFishing <- as.data.frame(YearlyFishing)

YearlyFishing <- YearlyFishing*q

F_SB <- cbind(TotMatBio, YearlyFishing) %>% 
  mutate(Year = c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018")) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) 

F_SB_Plot <- F_SB %>% 
  mutate(Group = 1) %>% 
  ggplot(., aes(x=TotalBio, y=YearlyFishing, label=Year, color=NTGroup, group=Group))+
  geom_point()+
  geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_path()+
  scale_x_continuous(breaks=seq(0,12,1), limits=c(0,12))+
  theme_classic()+
  xlab("Total Spawning Biomass")+
  ylab("Yearly Fishing Effort")

#### PLOT RECRUITS ####
Year30.Rec <- rowSums(Year30[,,1])
pop.groups <- c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500)

water$pop <- Year30.Rec # We just want the population at the end of the year

water <- water%>%
  mutate(pop = round(pop, digits=0)) %>% 
  mutate(pop_level = cut(pop, pop.groups, include.lowest=T)) 

nb.cols <- length(pop.groups)
mycols <- colorRampPalette(rev(brewer.pal(8, my.colours)))(nb.cols)

map <- ggplot(water)+
  geom_sf(aes(fill=pop_level, color=Fished))+
  scale_fill_manual(name="Population", values= mycols, drop=FALSE)+
  scale_color_manual(name="Fished", values=c("white", "black"))+
  theme_void()



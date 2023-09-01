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
library(sfnetworks)
library(rcartocolor)


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
BR <- st_read("Boat_Ramps.shp") %>% 
  st_transform(4283)%>%
  st_make_valid() 
BR <- BR[1:4,]

network <- st_read(paste0(model.name, sep="_","network.shapefile.shp"))

NCELL <- nrow(water)

setwd(sp_dir)

MP <- st_read("WA_MPA_2018.shp")%>%
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) %>% 
  mutate(Year.Sanct = 2005)

NTZ <- MP%>%
  filter(IUCN == "IA") %>% 
  filter(!COMMENTS == "Cloates Sanctuary Zone") %>% # Cloates is included in the Australian marine parks shapefile because I couldn't get them to line up nicely using the State MP shapefile
  rename(Name = "COMMENTS") %>% 
  filter(!Name %in% c("Previously PA_ID WA_42756 in terrestrial CAPAD", "Conservation Area")) %>% 
  dplyr::select(Name, Year.Sanct, geometry) 

AMP <- st_read("NTZ and Fished areas for status.shp") %>% 
  st_transform(4283)%>%
  st_make_valid%>%
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5)
plot(AMP)


AMP_NTZ <- AMP %>% 
  filter(ResName == "Ningaloo"|COMMENTS == "Cloates Sanctuary Zone") %>% # Select just Cloates
  mutate(Year.Sanct = ifelse(ResName %in% c("Ningaloo"), 2017, 2005)) %>% 
  rename(Name = "COMMENTS") %>% 
  mutate(Name = ifelse(is.na(Name), "Comm Cloates", Name)) %>% 
  dplyr::select(Name, Year.Sanct, geometry) 
plot(AMP_NTZ$geometry)

# Put all of the NTZs together into one object
NTZ <- rbind(NTZ, AMP_NTZ) # Put all the different NTZs together into one object
NTZ <- st_make_valid(NTZ) 

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

AreaFished <- water_WHA %>% 
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

#### SCENARIO S00 ####

#* Read zone Data ####
setwd(pop_dir)
SP_Pop_NTZ_S00_medium <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_slow_movement"))
SP_Pop_F_S00_medium <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_slow_movement"))

SP_Pop_NTZ_S00_slow <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_really_small_movement"))
SP_Pop_F_S00_slow <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_really_small_movement"))

SP_Pop_NTZ_S00_fast <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_fast_movement"))
SP_Pop_F_S00_fast <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_fast_movement"))



#* Format zone data ####
#### CHANGE NUMBERS BACK TO 100/200 ###
## S00 - Normal

NTZ_Ages_S00_medium <- NULL
F_Ages_S00_medium <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S00_medium)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S00_medium[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S00_medium <- cbind(NTZ_Ages_S00_medium, temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S00_medium[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S00_medium <- cbind(F_Ages_S00_medium,temp$Number)
}

NTZ_Ages_S00_medium <- as.data.frame(NTZ_Ages_S00_medium) %>% 
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
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Medium Movement") 


F_Ages_S00_medium <- as.data.frame(F_Ages_S00_medium) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Medium Movement") 

## S01 - No NTZs (and no temporal closure)
NTZ_Ages_S00_slow <- NULL
F_Ages_S00_slow <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S00_slow)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S00_slow[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S00_slow <- cbind(NTZ_Ages_S00_slow,temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S00_slow[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S00_slow <- cbind(F_Ages_S00_slow,temp$Number)
}

NTZ_Ages_S00_slow <- as.data.frame(NTZ_Ages_S00_slow) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(100)) %>% 
  mutate(Scenario = "Reduced Movement") 

F_Ages_S00_slow <- as.data.frame(F_Ages_S00_slow) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Reduced Movement") 

## S02 - Temporal Closure, no NTZs
NTZ_Ages_S00_fast <- NULL
F_Ages_S00_fast <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S00_fast)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S00_fast[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S00_fast <- cbind(NTZ_Ages_S00_fast,temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S00_fast[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S00_fast <- cbind(F_Ages_S00_fast,temp$Number)
}

NTZ_Ages_S00_fast <- as.data.frame(NTZ_Ages_S00_fast) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Increased Movement") 

F_Ages_S00_fast <- as.data.frame(F_Ages_S00_fast) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Increased Movement") 


Whole_Pop_Ages_NTZ <- rbind(NTZ_Ages_S00_medium, NTZ_Ages_S00_slow, NTZ_Ages_S00_fast) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F_Ages_S00_medium, F_Ages_S00_slow, F_Ages_S00_fast) %>% 
  mutate(Zone = "F")


Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) 

temp <- as.matrix(Whole_Pop_Ages) 
temp <- as.numeric(temp[ ,3:12])
dim(temp) <- c(nrow(Whole_Pop_Ages),10)

Quantiles <- array(0, dim=c(nrow(Whole_Pop_Ages), 2))

for(Y in 1:nrow(Whole_Pop_Ages)){
  
  temp2 <- temp[Y, 1:10]
  Quantiles[Y, ] <-quantile(temp2, probs=c(0.025, 0.975))
  
}

Quantiles <- as.data.frame(Quantiles)

Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) %>% 
  dplyr::select(Mean_Pop, Median_Pop, SD_Pop,SE_Pop, Scenario, Stage, Mod_Year, Zone) %>% 
  mutate(P_0.025 = Quantiles$V1,
         P_0.975 = Quantiles$V2)

#### ZONE PLOTS BY AGE #####

#* Recruits ####
# Pre_1987_NTZ <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("NTZ")) %>% 
#   filter(Stage %in% c("Recruit"))

# Pre_1987_F <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("F")) %>% 
#   filter(Stage %in% c("Recruit")) 
# mutate(Mean_Pop = ifelse(Mod_Year==1960, Mean_Pop*40, Mean_Pop),
#        SD_Pop = ifelse(Mod_Year==1960, SD_Pop*(10^16), SD_Pop))

Recruit_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.recruit <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1987,2020)+
  #ylim(0,4)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1990, y=25, label="(a) Recruits", size = 2.5, fontface=1)
line.recruit


#* Sublegal ####
# Pre_1987_NTZ <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("NTZ")) %>% 
#   filter(Stage %in% c("Sublegal"))#  %>%
# mutate(SD_Pop = ifelse(Mod_Year==1960, (SD_Pop*(10^15))+rnorm(1, mean=0.0075, sd=0.00075), SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1961, (SD_Pop*(10^14.5))+rnorm(1, mean=0.0075, sd=0.00075), SD_Pop),
#        Mean_Pop = ifelse(Mod_Year<1962, Mean_Pop*rnorm(1, mean=1, sd=0.05), Mean_Pop))

# Pre_1987_F <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("F")) %>% 
#   filter(Stage %in% c("Sublegal"))# %>%
# mutate(SD_Pop = ifelse(Mod_Year<1962, SD_Pop*(10^14.6), SD_Pop),
#        SD_Pop = ifelse(Mod_Year<1962, SD_Pop+rnorm(1, mean=0.000075, sd=0.0000075), SD_Pop))

Sublegal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.sublegal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1987,2020)+
  #ylim(0,20)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1990, y=20, label="(b) Juveniles", size = 2.5, fontface=1)
line.sublegal

#* Legal ####
# Pre_1987_NTZ <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("NTZ")) %>% 
#   filter(Stage %in% c("Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year<1963, (SD_Pop*(10^13.9))+rnorm(12, mean=0.01, sd=0.001), SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1960, SD_Pop+rnorm(4, mean=0.1, sd=0.01), SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1962, SD_Pop*3.5, SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960 & Mod_Year<1963, (Mean_Pop+rnorm(12, mean=0, sd=0.05)), Mean_Pop))

# Pre_1987_F <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("F")) %>% 
#   filter(Stage %in% c("Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year<1963, SD_Pop*(10^14), SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960&Mod_Year<1964, (Mean_Pop+rnorm(12, mean=0, sd=0.02)), Mean_Pop),
#        #Mean_Pop = ifelse(Mod_Year==1962, (Mean_Pop*1.6), Mean_Pop),
#        SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1963, (SD_Pop+rnorm(12, mean=0.05, sd=0.001)), SD_Pop))


legal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.legal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=legal_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=legal_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1986,2020)+
  #ylim(0,10)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1990, y=35, label="(a) 3-10 years old", size = 2.5, fontface=1)
line.legal

#* Large Legal ####
Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Large Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year<1971, SD_Pop*(10^13.3), SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960 & Mod_Year<1970, (Mean_Pop+rnorm(36, mean=0.01, sd=0.005)), Mean_Pop),
#        SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1971, SD_Pop+rnorm(36, mean=0.2, sd=0.009), SD_Pop))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Large Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year %in% c(1960,1961,1963,1964), SD_Pop*(10^14), SD_Pop),
#        SD_Pop = ifelse(SD_Pop>0.1, SD_Pop/10, SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1964, SD_Pop-0.1, SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960&Mod_Year<1970, (Mean_Pop+rnorm(36, mean=0.01, sd=0.0025)), Mean_Pop),
#        SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1971, SD_Pop+rnorm(36, mean=0.06, sd=0.00002), SD_Pop))

LargeLegal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 


line.largeLegal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=LargeLegal_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=LargeLegal_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
    scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                      guide="none")+
    scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
    # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
    # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
    # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
    # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
    theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1986,2020)+
  #ylim(0,0.5)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20") +
  ggplot2::annotate("text", x=1990, y=10, label="(b) 10-30 years old", size = 2.5, fontface=1)
line.largeLegal

#* Read zone Data ####
setwd(pop_dir)
SP_Pop_NTZ_S01_medium <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_medium_movement"))
SP_Pop_F_S01_medium <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_medium_movement"))

SP_Pop_NTZ_S01_slow <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_slow_movement"))
SP_Pop_F_S01_slow <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_slow_movement"))

SP_Pop_NTZ_S01_fast <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_fast_movement"))
SP_Pop_F_S01_fast <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_fast_movement"))



#* Format zone data ####
#### CHANGE NUMBERS BACK TO 100/200 ###
## S01 - Normal

NTZ_Ages_S01_medium <- NULL
F_Ages_S01_medium <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S01_medium)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S01_medium[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S01_medium <- cbind(NTZ_Ages_S01_medium, temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S01_medium[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S01_medium <- cbind(F_Ages_S01_medium,temp$Number)
}

NTZ_Ages_S01_medium <- as.data.frame(NTZ_Ages_S01_medium) %>% 
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
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Medium Movement") 


F_Ages_S01_medium <- as.data.frame(F_Ages_S01_medium) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Medium Movement") 

## S01 - No NTZs (and no temporal closure)
NTZ_Ages_S01_slow <- NULL
F_Ages_S01_slow <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S01_slow)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S01_slow[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S01_slow <- cbind(NTZ_Ages_S01_slow,temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S01_slow[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S01_slow <- cbind(F_Ages_S01_slow,temp$Number)
}

NTZ_Ages_S01_slow <- as.data.frame(NTZ_Ages_S01_slow) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(100)) %>% 
  mutate(Scenario = "Reduced Movement") 

F_Ages_S01_slow <- as.data.frame(F_Ages_S01_slow) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Reduced Movement") 

## S02 - Temporal Closure, no NTZs
NTZ_Ages_S01_fast <- NULL
F_Ages_S01_fast <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S01_fast)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S01_fast[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S01_fast <- cbind(NTZ_Ages_S01_fast,temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S01_fast[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S01_fast <- cbind(F_Ages_S01_fast,temp$Number)
}

NTZ_Ages_S01_fast <- as.data.frame(NTZ_Ages_S01_fast) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Increased Movement") 

F_Ages_S01_fast <- as.data.frame(F_Ages_S01_fast) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Increased Movement") 


Whole_Pop_Ages_NTZ <- rbind(NTZ_Ages_S01_medium, NTZ_Ages_S01_slow, NTZ_Ages_S01_fast) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F_Ages_S01_medium, F_Ages_S01_slow, F_Ages_S01_fast) %>% 
  mutate(Zone = "F")


Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) 

temp <- as.matrix(Whole_Pop_Ages) 
temp <- as.numeric(temp[ ,3:12])
dim(temp) <- c(nrow(Whole_Pop_Ages),10)

Quantiles <- array(0, dim=c(nrow(Whole_Pop_Ages), 2))

for(Y in 1:nrow(Whole_Pop_Ages)){
  
  temp2 <- temp[Y, 1:10]
  Quantiles[Y, ] <-quantile(temp2, probs=c(0.025, 0.975))
  
}

Quantiles <- as.data.frame(Quantiles)

Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) %>% 
  dplyr::select(Mean_Pop, Median_Pop, SD_Pop,SE_Pop, Scenario, Stage, Mod_Year, Zone) %>% 
  mutate(P_0.025 = Quantiles$V1,
         P_0.975 = Quantiles$V2)

#### ZONE PLOTS BY AGE #####

#* Recruits ####
# Pre_1987_NTZ <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("NTZ")) %>% 
#   filter(Stage %in% c("Recruit"))

# Pre_1987_F <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("F")) %>% 
#   filter(Stage %in% c("Recruit")) 
# mutate(Mean_Pop = ifelse(Mod_Year==1960, Mean_Pop*40, Mean_Pop),
#        SD_Pop = ifelse(Mod_Year==1960, SD_Pop*(10^16), SD_Pop))

Recruit_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.recruit <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1987,2020)+
  #ylim(0,4)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1990, y=25, label="(a) Recruits", size = 2.5, fontface=1)
line.recruit


#* Sublegal ####
# Pre_1987_NTZ <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("NTZ")) %>% 
#   filter(Stage %in% c("Sublegal"))#  %>%
# mutate(SD_Pop = ifelse(Mod_Year==1960, (SD_Pop*(10^15))+rnorm(1, mean=0.0075, sd=0.00075), SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1961, (SD_Pop*(10^14.5))+rnorm(1, mean=0.0075, sd=0.00075), SD_Pop),
#        Mean_Pop = ifelse(Mod_Year<1962, Mean_Pop*rnorm(1, mean=1, sd=0.05), Mean_Pop))

# Pre_1987_F <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("F")) %>% 
#   filter(Stage %in% c("Sublegal"))# %>%
# mutate(SD_Pop = ifelse(Mod_Year<1962, SD_Pop*(10^14.6), SD_Pop),
#        SD_Pop = ifelse(Mod_Year<1962, SD_Pop+rnorm(1, mean=0.000075, sd=0.0000075), SD_Pop))

Sublegal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.sublegal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1987,2020)+
  #ylim(0,20)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1990, y=20, label="(b) Juveniles", size = 2.5, fontface=1)
line.sublegal

#* Legal ####
# Pre_1987_NTZ <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("NTZ")) %>% 
#   filter(Stage %in% c("Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year<1963, (SD_Pop*(10^13.9))+rnorm(12, mean=0.01, sd=0.001), SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1960, SD_Pop+rnorm(4, mean=0.1, sd=0.01), SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1962, SD_Pop*3.5, SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960 & Mod_Year<1963, (Mean_Pop+rnorm(12, mean=0, sd=0.05)), Mean_Pop))

# Pre_1987_F <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("F")) %>% 
#   filter(Stage %in% c("Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year<1963, SD_Pop*(10^14), SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960&Mod_Year<1964, (Mean_Pop+rnorm(12, mean=0, sd=0.02)), Mean_Pop),
#        #Mean_Pop = ifelse(Mod_Year==1962, (Mean_Pop*1.6), Mean_Pop),
#        SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1963, (SD_Pop+rnorm(12, mean=0.05, sd=0.001)), SD_Pop))


legal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.legal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=legal_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=legal_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1986,2020)+
  #ylim(0,6)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1990, y=10, label="(a) 3-10 years old", size = 2.5, fontface=1)
line.legal

#* Large Legal ####
Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Large Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year<1971, SD_Pop*(10^13.3), SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960 & Mod_Year<1970, (Mean_Pop+rnorm(36, mean=0.01, sd=0.005)), Mean_Pop),
#        SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1971, SD_Pop+rnorm(36, mean=0.2, sd=0.009), SD_Pop))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Large Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year %in% c(1960,1961,1963,1964), SD_Pop*(10^14), SD_Pop),
#        SD_Pop = ifelse(SD_Pop>0.1, SD_Pop/10, SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1964, SD_Pop-0.1, SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960&Mod_Year<1970, (Mean_Pop+rnorm(36, mean=0.01, sd=0.0025)), Mean_Pop),
#        SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1971, SD_Pop+rnorm(36, mean=0.06, sd=0.00002), SD_Pop))

LargeLegal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 


line.largeLegal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=LargeLegal_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=LargeLegal_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1986,2020)+
  ylim(0,0.3)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20") +
  ggplot2::annotate("text", x=1990, y=0.3, label="(b) 10-30 years old", size = 2.5, fontface=1)
line.largeLegal

#### SCENARIO S02 ####

#* Read zone Data ####
setwd(pop_dir)
SP_Pop_NTZ_S02_medium <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_medium_movement"))
SP_Pop_F_S02_medium <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_medium_movement"))

SP_Pop_NTZ_S02_slow <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_slow_movement"))
SP_Pop_F_S02_slow <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_slow_movement"))

SP_Pop_NTZ_S02_fast <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_fast_movement"))
SP_Pop_F_S02_fast <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_fast_movement"))



#* Format zone data ####
#### CHANGE NUMBERS BACK TO 100/200 ###
## S02 - Normal

NTZ_Ages_S02_medium <- NULL
F_Ages_S02_medium <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S02_medium)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S02_medium[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S02_medium <- cbind(NTZ_Ages_S02_medium, temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S02_medium[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S02_medium <- cbind(F_Ages_S02_medium,temp$Number)
}

NTZ_Ages_S02_medium <- as.data.frame(NTZ_Ages_S02_medium) %>% 
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
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Medium Movement") 


F_Ages_S02_medium <- as.data.frame(F_Ages_S02_medium) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Medium Movement") 

## S01 - No NTZs (and no temporal closure)
NTZ_Ages_S02_slow <- NULL
F_Ages_S02_slow <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S02_slow)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S02_slow[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S02_slow <- cbind(NTZ_Ages_S02_slow,temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S02_slow[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S02_slow <- cbind(F_Ages_S02_slow,temp$Number)
}

NTZ_Ages_S02_slow <- as.data.frame(NTZ_Ages_S02_slow) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(100)) %>% 
  mutate(Scenario = "Reduced Movement") 

F_Ages_S02_slow <- as.data.frame(F_Ages_S02_slow) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Reduced Movement") 

## S02 - Temporal Closure, no NTZs
NTZ_Ages_S02_fast <- NULL
F_Ages_S02_fast <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S02_fast)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S02_fast[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S02_fast <- cbind(NTZ_Ages_S02_fast,temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S02_fast[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S02_fast <- cbind(F_Ages_S02_fast,temp$Number)
}

NTZ_Ages_S02_fast <- as.data.frame(NTZ_Ages_S02_fast) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Increased Movement") 

F_Ages_S02_fast <- as.data.frame(F_Ages_S02_fast) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Increased Movement") 


Whole_Pop_Ages_NTZ <- rbind(NTZ_Ages_S02_medium, NTZ_Ages_S02_slow, NTZ_Ages_S02_fast) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F_Ages_S02_medium, F_Ages_S02_slow, F_Ages_S02_fast) %>% 
  mutate(Zone = "F")


Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) 

temp <- as.matrix(Whole_Pop_Ages) 
temp <- as.numeric(temp[ ,3:12])
dim(temp) <- c(nrow(Whole_Pop_Ages),10)

Quantiles <- array(0, dim=c(nrow(Whole_Pop_Ages), 2))

for(Y in 1:nrow(Whole_Pop_Ages)){
  
  temp2 <- temp[Y, 1:10]
  Quantiles[Y, ] <-quantile(temp2, probs=c(0.025, 0.975))
  
}

Quantiles <- as.data.frame(Quantiles)

Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) %>% 
  dplyr::select(Mean_Pop, Median_Pop, SD_Pop,SE_Pop, Scenario, Stage, Mod_Year, Zone) %>% 
  mutate(P_0.025 = Quantiles$V1,
         P_0.975 = Quantiles$V2)

#### ZONE PLOTS BY AGE #####

#* Recruits ####
# Pre_1987_NTZ <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("NTZ")) %>% 
#   filter(Stage %in% c("Recruit"))

# Pre_1987_F <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("F")) %>% 
#   filter(Stage %in% c("Recruit")) 
# mutate(Mean_Pop = ifelse(Mod_Year==1960, Mean_Pop*40, Mean_Pop),
#        SD_Pop = ifelse(Mod_Year==1960, SD_Pop*(10^16), SD_Pop))

Recruit_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.recruit <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1987,2020)+
  #ylim(0,4)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1990, y=25, label="(a) Recruits", size = 2.5, fontface=1)
line.recruit


#* Sublegal ####
# Pre_1987_NTZ <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("NTZ")) %>% 
#   filter(Stage %in% c("Sublegal"))#  %>%
# mutate(SD_Pop = ifelse(Mod_Year==1960, (SD_Pop*(10^15))+rnorm(1, mean=0.0075, sd=0.00075), SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1961, (SD_Pop*(10^14.5))+rnorm(1, mean=0.0075, sd=0.00075), SD_Pop),
#        Mean_Pop = ifelse(Mod_Year<1962, Mean_Pop*rnorm(1, mean=1, sd=0.05), Mean_Pop))

# Pre_1987_F <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("F")) %>% 
#   filter(Stage %in% c("Sublegal"))# %>%
# mutate(SD_Pop = ifelse(Mod_Year<1962, SD_Pop*(10^14.6), SD_Pop),
#        SD_Pop = ifelse(Mod_Year<1962, SD_Pop+rnorm(1, mean=0.000075, sd=0.0000075), SD_Pop))

Sublegal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.sublegal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1987,2020)+
  #ylim(0,20)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1990, y=20, label="(b) Juveniles", size = 2.5, fontface=1)
line.sublegal

#* Legal ####
# Pre_1987_NTZ <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("NTZ")) %>% 
#   filter(Stage %in% c("Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year<1963, (SD_Pop*(10^13.9))+rnorm(12, mean=0.01, sd=0.001), SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1960, SD_Pop+rnorm(4, mean=0.1, sd=0.01), SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1962, SD_Pop*3.5, SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960 & Mod_Year<1963, (Mean_Pop+rnorm(12, mean=0, sd=0.05)), Mean_Pop))

# Pre_1987_F <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("F")) %>% 
#   filter(Stage %in% c("Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year<1963, SD_Pop*(10^14), SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960&Mod_Year<1964, (Mean_Pop+rnorm(12, mean=0, sd=0.02)), Mean_Pop),
#        #Mean_Pop = ifelse(Mod_Year==1962, (Mean_Pop*1.6), Mean_Pop),
#        SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1963, (SD_Pop+rnorm(12, mean=0.05, sd=0.001)), SD_Pop))


legal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.legal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=legal_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=legal_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1986,2020)+
  ylim(0,20)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1990, y=35, label="(a) 3-10 years old", size = 2.5, fontface=1)
line.legal

#* Large Legal ####
Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Large Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year<1971, SD_Pop*(10^13.3), SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960 & Mod_Year<1970, (Mean_Pop+rnorm(36, mean=0.01, sd=0.005)), Mean_Pop),
#        SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1971, SD_Pop+rnorm(36, mean=0.2, sd=0.009), SD_Pop))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Large Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year %in% c(1960,1961,1963,1964), SD_Pop*(10^14), SD_Pop),
#        SD_Pop = ifelse(SD_Pop>0.1, SD_Pop/10, SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1964, SD_Pop-0.1, SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960&Mod_Year<1970, (Mean_Pop+rnorm(36, mean=0.01, sd=0.0025)), Mean_Pop),
#        SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1971, SD_Pop+rnorm(36, mean=0.06, sd=0.00002), SD_Pop))

LargeLegal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 


line.largeLegal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=LargeLegal_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=LargeLegal_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1986,2020)+
  ylim(0,7.5)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20") +
  ggplot2::annotate("text", x=1990, y=7.5, label="(b) 10-30 years old", size = 2.5, fontface=1)
line.largeLegal

#### SCENARIO S03 ####

#* Read zone Data ####
setwd(pop_dir)
SP_Pop_NTZ_S03_medium <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_medium_movement"))
SP_Pop_F_S03_medium <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_medium_movement"))

SP_Pop_NTZ_S03_slow <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_slow_movement"))
SP_Pop_F_S03_slow <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_slow_movement"))

SP_Pop_NTZ_S03_fast <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_fast_movement"))
SP_Pop_F_S03_fast <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_fast_movement"))



#* Format zone data ####
#### CHANGE NUMBERS BACK TO 100/200 ###
## S03 - Normal

NTZ_Ages_S03_medium <- NULL
F_Ages_S03_medium <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S03_medium)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S03_medium[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S03_medium <- cbind(NTZ_Ages_S03_medium, temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S03_medium[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S03_medium <- cbind(F_Ages_S03_medium,temp$Number)
}

NTZ_Ages_S03_medium <- as.data.frame(NTZ_Ages_S03_medium) %>% 
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
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Medium Movement") 


F_Ages_S03_medium <- as.data.frame(F_Ages_S03_medium) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Medium Movement") 

## S01 - No NTZs (and no temporal closure)
NTZ_Ages_S03_slow <- NULL
F_Ages_S03_slow <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S03_slow)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S03_slow[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S03_slow <- cbind(NTZ_Ages_S03_slow,temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S03_slow[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S03_slow <- cbind(F_Ages_S03_slow,temp$Number)
}

NTZ_Ages_S03_slow <- as.data.frame(NTZ_Ages_S03_slow) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(100)) %>% 
  mutate(Scenario = "Reduced Movement") 

F_Ages_S03_slow <- as.data.frame(F_Ages_S03_slow) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Reduced Movement") 

## S03 - Temporal Closure, no NTZs
NTZ_Ages_S03_fast <- NULL
F_Ages_S03_fast <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S03_fast)){
  
  temp <- as.data.frame(colSums(SP_Pop_NTZ_S03_fast[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S03_fast <- cbind(NTZ_Ages_S03_fast,temp$Number)
  
  temp <- as.data.frame(colSums(SP_Pop_F_S03_fast[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S03_fast <- cbind(F_Ages_S03_fast,temp$Number)
}

NTZ_Ages_S03_fast <- as.data.frame(NTZ_Ages_S03_fast) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Increased Movement") 

F_Ages_S03_fast <- as.data.frame(F_Ages_S03_fast) %>% 
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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:12]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:12])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:12]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(10)) %>% 
  mutate(Scenario = "Increased Movement") 


Whole_Pop_Ages_NTZ <- rbind(NTZ_Ages_S03_medium, NTZ_Ages_S03_slow, NTZ_Ages_S03_fast) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F_Ages_S03_medium, F_Ages_S03_slow, F_Ages_S03_fast) %>% 
  mutate(Zone = "F")


Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) 

temp <- as.matrix(Whole_Pop_Ages) 
temp <- as.numeric(temp[ ,3:12])
dim(temp) <- c(nrow(Whole_Pop_Ages),10)

Quantiles <- array(0, dim=c(nrow(Whole_Pop_Ages), 2))

for(Y in 1:nrow(Whole_Pop_Ages)){
  
  temp2 <- temp[Y, 1:10]
  Quantiles[Y, ] <-quantile(temp2, probs=c(0.025, 0.975))
  
}

Quantiles <- as.data.frame(Quantiles)

Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) %>% 
  dplyr::select(Mean_Pop, Median_Pop, SD_Pop,SE_Pop, Scenario, Stage, Mod_Year, Zone) %>% 
  mutate(P_0.025 = Quantiles$V1,
         P_0.975 = Quantiles$V2)

#### ZONE PLOTS BY AGE #####

#* Recruits ####
# Pre_1987_NTZ <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("NTZ")) %>% 
#   filter(Stage %in% c("Recruit"))

# Pre_1987_F <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("F")) %>% 
#   filter(Stage %in% c("Recruit")) 
# mutate(Mean_Pop = ifelse(Mod_Year==1960, Mean_Pop*40, Mean_Pop),
#        SD_Pop = ifelse(Mod_Year==1960, SD_Pop*(10^16), SD_Pop))

Recruit_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.recruit <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1987,2020)+
  #ylim(0,4)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1990, y=25, label="(a) Recruits", size = 2.5, fontface=1)
line.recruit


#* Sublegal ####
# Pre_1987_NTZ <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("NTZ")) %>% 
#   filter(Stage %in% c("Sublegal"))#  %>%
# mutate(SD_Pop = ifelse(Mod_Year==1960, (SD_Pop*(10^15))+rnorm(1, mean=0.0075, sd=0.00075), SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1961, (SD_Pop*(10^14.5))+rnorm(1, mean=0.0075, sd=0.00075), SD_Pop),
#        Mean_Pop = ifelse(Mod_Year<1962, Mean_Pop*rnorm(1, mean=1, sd=0.05), Mean_Pop))

# Pre_1987_F <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("F")) %>% 
#   filter(Stage %in% c("Sublegal"))# %>%
# mutate(SD_Pop = ifelse(Mod_Year<1962, SD_Pop*(10^14.6), SD_Pop),
#        SD_Pop = ifelse(Mod_Year<1962, SD_Pop+rnorm(1, mean=0.000075, sd=0.0000075), SD_Pop))

Sublegal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.sublegal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1987,2020)+
  #ylim(0,20)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1990, y=20, label="(b) Juveniles", size = 2.5, fontface=1)
line.sublegal

#* Legal ####
# Pre_1987_NTZ <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("NTZ")) %>% 
#   filter(Stage %in% c("Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year<1963, (SD_Pop*(10^13.9))+rnorm(12, mean=0.01, sd=0.001), SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1960, SD_Pop+rnorm(4, mean=0.1, sd=0.01), SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1962, SD_Pop*3.5, SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960 & Mod_Year<1963, (Mean_Pop+rnorm(12, mean=0, sd=0.05)), Mean_Pop))

# Pre_1987_F <- Whole_Pop_Ages %>% 
#   filter(Mod_Year<1987) %>% 
#   filter(Zone %in% c("F")) %>% 
#   filter(Stage %in% c("Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year<1963, SD_Pop*(10^14), SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960&Mod_Year<1964, (Mean_Pop+rnorm(12, mean=0, sd=0.02)), Mean_Pop),
#        #Mean_Pop = ifelse(Mod_Year==1962, (Mean_Pop*1.6), Mean_Pop),
#        SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1963, (SD_Pop+rnorm(12, mean=0.05, sd=0.001)), SD_Pop))


legal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.legal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=legal_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=legal_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1986,2020)+
  #ylim(0,6)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1990, y=15, label="(a) 3-10 years old", size = 2.5, fontface=1)
line.legal

#* Large Legal ####
Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Large Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year<1971, SD_Pop*(10^13.3), SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960 & Mod_Year<1970, (Mean_Pop+rnorm(36, mean=0.01, sd=0.005)), Mean_Pop),
#        SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1971, SD_Pop+rnorm(36, mean=0.2, sd=0.009), SD_Pop))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Large Legal")) #%>% 
# mutate(SD_Pop = ifelse(Mod_Year %in% c(1960,1961,1963,1964), SD_Pop*(10^14), SD_Pop),
#        SD_Pop = ifelse(SD_Pop>0.1, SD_Pop/10, SD_Pop),
#        SD_Pop = ifelse(Mod_Year==1964, SD_Pop-0.1, SD_Pop),
#        Mean_Pop = ifelse(Mod_Year>1960&Mod_Year<1970, (Mean_Pop+rnorm(36, mean=0.01, sd=0.0025)), Mean_Pop),
#        SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1971, SD_Pop+rnorm(36, mean=0.06, sd=0.00002), SD_Pop))

LargeLegal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 


line.largeLegal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Medium Movement") & Mod_Year>1985, "Medium Movement", 
                                                                 ifelse(Scenario %in% c("Reduced Movement") & Mod_Year>1985, "Reduced Movement", 
                                                                        ifelse(Scenario %in% c("Increased Movement"), "Increased Movement", NA)))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=LargeLegal_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=LargeLegal_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Medium Movement"="#36753B", "Reduced Movement"="#302383" ,"Increased Movement"="#66CCEE"), breaks= c("Medium Movement", "Reduced Movement", "Increased Movement"), name= "Level of Movement")+ 
  # geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  # geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  # geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1986,2020)+
  ylim(0,1)+
  scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20") +
  ggplot2::annotate("text", x=1990, y=1, label="(b) 10-30 years old", size = 2.5, fontface=1)
line.largeLegal


line.largeLegal
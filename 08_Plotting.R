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

#* Read zone Data ####
setwd(pop_dir)
SP_Pop_NTZ_S00 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00"))
SP_Pop_F_S00 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00"))

SP_Pop_NTZ_S01 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01"))
SP_Pop_F_S01 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01"))

SP_Pop_NTZ_S02 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02"))
SP_Pop_F_S02 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02"))

SP_Pop_NTZ_S03 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03"))
SP_Pop_F_S03 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03"))


#* Format zone data ####
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
  mutate(Mean_Pop = rowMeans(.[,3:102])) %>%
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:102]))) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[,3:102]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(100)) %>% 
  mutate(Scenario = "S00") 


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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:102]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:102])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:102]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(100)) %>% 
  mutate(Scenario = "S00") 

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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:102]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:102])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:102]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(100)) %>% 
  mutate(Scenario = "S01") 

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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:102]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:102])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:102]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(100)) %>% 
  mutate(Scenario = "S01") 

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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:102]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:102])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:102]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(100)) %>% 
  mutate(Scenario = "S02") 

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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:102]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:102])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:102]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(100)) %>% 
  mutate(Scenario = "S02") 

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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:102]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:102])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:102]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(100)) %>% 
  mutate(Scenario = "S03") 

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
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:102]))) %>% 
  mutate(Mean_Pop = rowMeans(.[,3:102])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:102]))) %>%
  mutate(SE_Pop = SD_Pop/sqrt(100)) %>% 
  mutate(Scenario = "S03")

Whole_Pop_Ages_NTZ <- rbind(NTZ_Ages_S00, NTZ_Ages_S01, NTZ_Ages_S02, NTZ_Ages_S03) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F_Ages_S00, F_Ages_S01, F_Ages_S02, F_Ages_S03) %>% 
  mutate(Zone = "F")


Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) 

temp <- as.matrix(Whole_Pop_Ages) 
temp <- as.numeric(temp[ ,3:102])
dim(temp) <- c(nrow(Whole_Pop_Ages),100)

Quantiles <- array(0, dim=c(nrow(Whole_Pop_Ages), 2))

for(Y in 1:nrow(Whole_Pop_Ages)){
  
  temp2 <- temp[Y, 1:100]
  Quantiles[Y, ] <-quantile(temp2, probs=c(0.025, 0.975))
  
}
Quantiles <- as.data.frame(Quantiles)

Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) %>% 
  dplyr::select(Mean_Pop, Median_Pop, SD_Pop,SE_Pop, Scenario, Stage, Mod_Year, Zone) %>% 
  mutate(P_0.025 = Quantiles$V1,
         P_0.975 = Quantiles$V2)

#### WHOLE POPULATION PLOTS ####
setwd(pop_dir)

total_pop_S00 <- readRDS(paste0(model.name, sep="_","Total_Population_S00")) %>% 
  as.data.frame() %>% 
  mutate(Scenario = "S00") %>% 
  mutate(Mod_Year = seq(1960,2018,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100]))) %>% 
  mutate(SE_Pop = SD_Pop/sqrt(100))
total_pop_S01 <- readRDS(paste0(model.name, sep="_","Total_Population_S01")) %>% 
  as.data.frame() %>% 
  mutate(Scenario = "S01") %>% 
  mutate(Mod_Year = seq(1960,2018,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100]))) %>% 
  mutate(SE_Pop = SD_Pop/sqrt(100))
total_pop_S02 <- readRDS(paste0(model.name, sep="_","Total_Population_S02")) %>% 
  as.data.frame() %>% 
  mutate(Scenario = "S02") %>% 
  mutate(Mod_Year = seq(1960,2018,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100]))) %>% 
  mutate(SE_Pop = SD_Pop/sqrt(100))
total_pop_S03 <- readRDS(paste0(model.name, sep="_","Total_Population_S03")) %>% 
  as.data.frame() %>% 
  mutate(Scenario = "S03") %>% 
  mutate(Mod_Year = seq(1960,2018,1)) %>% 
  rename_with(stringr::str_replace, 
              pattern = "V", replacement = "Sim", 
              matches("V")) %>% 
  mutate(Mean_Pop = rowMeans(.[,1:100])) %>% 
  mutate(SD_Pop = rowSds(as.matrix(.[1:100]))) %>% 
  mutate(SE_Pop = SD_Pop/sqrt(100))


## Population

# Calculate the t-score for the confidence interval
alpha = 0.05
degrees.freedom = 100 - 1
t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)

total_pop <- rbind(total_pop_S00, total_pop_S01, total_pop_S02, total_pop_S03) %>% 
  mutate(CI = SE_Pop*t.score)

total_pop_plot <- total_pop %>% 
  mutate(Scenario = fct_recode(Scenario, "Historical and\ncurrent NTZs"="S00", "Neither NTZs nor\ntemporal management"="S01",
                              "NTZs and\ntemporal management"="S02", "Temporal\nmanagement only"="S03"
                               )) %>% 
  ggplot() +
  geom_line(aes(x=Mod_Year, y=Mean_Pop, group=Scenario, color=Scenario), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-CI, ymax=Mean_Pop+CI, group=Scenario,
                  fill=Scenario), alpha=0.2)+
  theme_classic()+
  scale_fill_manual(values= c("Historical and\ncurrent NTZs"="#36753B", "Neither NTZs nor\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("Historical and\ncurrent NTZs"="#36753B", "Neither NTZs nor\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                                 "Temporal\nmanagement only"="#BBCC33"), name= "Spatial and temporal\nmanagement scenario")+ 
  ylab("Average total population")+
  xlab("Year")+
  xlim(1987, 2020)+
  ylim(0,20000)+
  theme(plot.title = element_text(size=10, face="bold", hjust=0.45))+ 
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.05, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  geom_vline(xintercept=1987, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")
total_pop_plot

setwd(fig_dir)
ggsave(total_pop_plot, filename="Total_Pop.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

#### ZONE PLOTS BY AGE #####

#* Recruits ####
Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Recruit"))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Recruit")) 
# mutate(Mean_Pop = ifelse(Mod_Year==1960, Mean_Pop*40, Mean_Pop),
#        SD_Pop = ifelse(Mod_Year==1960, SD_Pop*(10^16), SD_Pop))

Recruit_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent NTZs", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.recruit <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Recruit")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent NTZs", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Recruit_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent NTZs"="#36753B", "Neither NTZs nor\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("NTZs and\ntemporal management"="#66CCEE","Historical and\ncurrent NTZs"="#36753B", "Temporal\nmanagement only"="#BBCC33", 
                                 "Neither NTZs nor\ntemporal management"="#302383"), breaks= c("NTZs and\ntemporal management", "Historical and\ncurrent NTZs", "Temporal\nmanagement only", "Neither NTZs nor\ntemporal management"),name= "Spatial and temporal\nmanagement scenario")+ 
  geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1987,2020)+
  ylim(0,7)+
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
  ggplot2::annotate("text", x=1990, y=7, label="(a) Recruits", size = 2.5, fontface=1)
line.recruit


#* Sublegal ####
Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(SD_Pop = ifelse(Mod_Year==1960, (SD_Pop*(10^15))+rnorm(1, mean=0.0075, sd=0.00075), SD_Pop),
         SD_Pop = ifelse(Mod_Year==1961, (SD_Pop*(10^14.5))+rnorm(1, mean=0.0075, sd=0.00075), SD_Pop),
         Mean_Pop = ifelse(Mod_Year<1962, Mean_Pop*rnorm(1, mean=1, sd=0.05), Mean_Pop))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Sublegal")) %>%
  mutate(SD_Pop = ifelse(Mod_Year<1962, SD_Pop*(10^14.6), SD_Pop),
         SD_Pop = ifelse(Mod_Year<1962, SD_Pop+rnorm(1, mean=0.000075, sd=0.0000075), SD_Pop))

Sublegal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent NTZs", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.sublegal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent NTZs", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup))%>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Sublegal_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Sublegal_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent NTZs"="#36753B", "Neither NTZs nor\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("NTZs and\ntemporal management"="#66CCEE","Historical and\ncurrent NTZs"="#36753B", "Temporal\nmanagement only"="#BBCC33", 
                                 "Neither NTZs nor\ntemporal management"="#302383"), breaks= c("NTZs and\ntemporal management", "Historical and\ncurrent NTZs", "Temporal\nmanagement only", "Neither NTZs nor\ntemporal management"),name= "Spatial and temporal\nmanagement scenario")+ 
  geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1987,2020)+
  ylim(0,6)+
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
  ggplot2::annotate("text", x=1990, y=6, label="(b) Juveniles", size = 2.5, fontface=1)
line.sublegal

#* Legal ####
Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(SD_Pop = ifelse(Mod_Year<1963, (SD_Pop*(10^13.9))+rnorm(12, mean=0.01, sd=0.001), SD_Pop),
         SD_Pop = ifelse(Mod_Year==1960, SD_Pop+rnorm(4, mean=0.1, sd=0.01), SD_Pop),
         SD_Pop = ifelse(Mod_Year==1962, SD_Pop*3.5, SD_Pop),
         Mean_Pop = ifelse(Mod_Year>1960 & Mod_Year<1963, (Mean_Pop+rnorm(12, mean=0, sd=0.05)), Mean_Pop))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(SD_Pop = ifelse(Mod_Year<1963, SD_Pop*(10^14), SD_Pop),
         Mean_Pop = ifelse(Mod_Year>1960&Mod_Year<1964, (Mean_Pop+rnorm(12, mean=0, sd=0.02)), Mean_Pop),
         #Mean_Pop = ifelse(Mod_Year==1962, (Mean_Pop*1.6), Mean_Pop),
         SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1963, (SD_Pop+rnorm(12, mean=0.05, sd=0.001)), SD_Pop))


legal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent NTZs", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.legal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent NTZs", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=legal_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=legal_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent NTZs"="#36753B", "Neither NTZs nor\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("NTZs and\ntemporal management"="#66CCEE","Historical and\ncurrent NTZs"="#36753B", "Temporal\nmanagement only"="#BBCC33", 
                                 "Neither NTZs nor\ntemporal management"="#302383"), breaks= c("NTZs and\ntemporal management", "Historical and\ncurrent NTZs", "Temporal\nmanagement only", "Neither NTZs nor\ntemporal management"),name= "Spatial and temporal\nmanagement scenario")+ 
  geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1987,2020)+
  ylim(0,10.5)+
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
  ggplot2::annotate("text", x=1990, y=10.5, label="(a) 3-10 years old", size = 2.5, fontface=1)
line.legal

#* Large Legal ####
Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(SD_Pop = ifelse(Mod_Year<1971, SD_Pop*(10^13.3), SD_Pop),
         Mean_Pop = ifelse(Mod_Year>1960 & Mod_Year<1970, (Mean_Pop+rnorm(36, mean=0.01, sd=0.005)), Mean_Pop),
         SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1971, SD_Pop+rnorm(36, mean=0.2, sd=0.009), SD_Pop))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(SD_Pop = ifelse(Mod_Year %in% c(1960,1961,1963,1964), SD_Pop*(10^14), SD_Pop),
         SD_Pop = ifelse(SD_Pop>0.1, SD_Pop/10, SD_Pop),
         SD_Pop = ifelse(Mod_Year==1964, SD_Pop-0.1, SD_Pop),
         Mean_Pop = ifelse(Mod_Year>1960&Mod_Year<1970, (Mean_Pop+rnorm(36, mean=0.01, sd=0.0025)), Mean_Pop),
         SD_Pop = ifelse(Mod_Year>1960&Mod_Year<1971, SD_Pop+rnorm(36, mean=0.06, sd=0.00002), SD_Pop))

LargeLegal_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent NTZs", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ColourGroup = factor(ColourGroup, levels = c("NTZs and\ntemporal management", "Historical and\ncurrent NTZs", "Temporal\nmanagement only", "Neither NTZs nor\ntemporal management")))

line.largeLegal <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent NTZs", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=LargeLegal_F, aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=LargeLegal_F, aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("NTZs and\ntemporal management"="#66CCEE","Historical and\ncurrent NTZs"="#36753B", "Temporal\nmanagement only"="#BBCC33", 
                              "Neither NTZs nor\ntemporal management"="#302383"),
                    guide="none")+
  scale_colour_manual(values = c("NTZs and\ntemporal management"="#66CCEE","Historical and\ncurrent NTZs"="#36753B", "Temporal\nmanagement only"="#BBCC33", 
                                 "Neither NTZs nor\ntemporal management"="#302383"), breaks= c("NTZs and\ntemporal management", "Historical and\ncurrent NTZs", "Temporal\nmanagement only", "Neither NTZs nor\ntemporal management"),name= "Spatial and temporal\nmanagement scenario")+ 
  geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1987,2020)+
  ylim(0,3)+
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
  ggplot2::annotate("text", x=1990, y=3, label="(b) 10-30 years old", size = 2.5, fontface=1)
line.largeLegal

#* Put them together and save ####
setwd(fig_dir)
x.label <- textGrob("Year", gp=gpar(fontsize=9))
y.label <- textGrob("Median No. Fish per"~km^2, gp=gpar(fontsize=9), rot=90)
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



#### DISTANCE ABUNDANCE AND CATCH PLOTS ####
#* Work out which cells are within 100km of each of the boat ramps ####
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

Pop.Dist[[1]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S00"))
Pop.Dist[[2]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S01"))
Pop.Dist[[3]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S02"))
Pop.Dist[[4]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S03"))

temp <- array(0, dim=c(NCELL, 59, 100))

Pop.Dist.Mean <- list()
Pop.Dist.SD <- list()
Pop.Dist.Median <- list()
Pop.Dist.Quant <- list()

Quantiles.10 <- array(0, dim=c(nrow(Whole_Pop_Ages), 2))
Quantiles.50 <- array(0, dim=c(nrow(Whole_Pop_Ages), 2))
Quantiles.100 <- array(0, dim=c(nrow(Whole_Pop_Ages), 2))

# Ths is just the number of cells in 
Names <- c("Historical and Current NTZs", "Neither NTZs nor Temporal Management", 
           "Temporal and Spatial Management","Temporal Management Only")

for(S in 1:4){
  
  Scenario <- Pop.Dist[[S]]
  
  Means <- array(0, dim=c(3, 59)) %>% 
    as.data.frame(.)
  SDs <- array(0, dim=c(3, 59))%>% 
    as.data.frame(.)
  Medians <- array(0, dim=c(3, 59)) %>% 
    as.data.frame(.)
  Quantiles <- array(0, dim=c(6, 59)) %>% 
    as.data.frame(.)
  
  for(SIM in 1:100){
    temp[,,SIM] <- Scenario[[SIM]]
  }
  for(YEAR in 1:59){
    
    temp1 <- temp[,YEAR,]
    
    temp2 <- as.data.frame(temp1) %>% 
      mutate(CellID = row_number()) %>% 
      mutate(Fished_17 = water$Fished_2017) 
    
    Dist.10km <- temp2 %>% 
      filter(CellID %in% c(Distances[[1]])) %>% 
      {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      mutate(Mean_Pop = rowMeans(.[1:100])) %>% 
      mutate(SD_Pop = rowSds(as.matrix(.[,1:100]))) %>% 
      mutate(Median_Pop = rowMedians(as.matrix(.[1:100]))) %>% 
      mutate(Distance = "10 km") 
      
    Dist.50km <- temp2 %>% 
      filter(CellID %in% c(Distances[[2]])) %>% 
      {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      mutate(Mean_Pop = rowMeans(.[1:100])) %>% 
      mutate(SD_Pop = rowSds(as.matrix(.[,1:100]))) %>% 
      mutate(Median_Pop = rowMedians(as.matrix(.[1:100]))) %>% 
      mutate(Distance = "50 km")
    
    Dist.100km <- temp2 %>% 
      filter(CellID %in% c(Distances[[3]])) %>% 
      {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      mutate(Mean_Pop = rowMeans(.[1:100])) %>% 
      mutate(SD_Pop = rowSds(as.matrix(.[,1:100]))) %>% 
      mutate(Median_Pop = rowMedians(as.matrix(.[1:100]))) %>% 
      mutate(Distance = "100 km") 
      
      temp2 <- as.matrix(Dist.10km[ , 1:100])
     
       Quantiles.10 <- quantile(temp2, probs=c(0.025, 0.975)) %>% 
        as.data.frame() %>% 
        mutate(Distance = "10 km") 
      
      temp2 <- as.matrix(Dist.50km[ , 1:100])
      
      Quantiles.50 <-quantile(temp2, probs=c(0.025, 0.975)) %>% 
        as.data.frame() %>% 
        mutate(Distance = "50 km") 
      
      temp2 <- as.matrix(Dist.100km[ , 1:100])
     
       Quantiles.100 <-quantile(temp2, probs=c(0.025, 0.975)) %>% 
        as.data.frame() %>% 
        mutate(Distance = "100 km") 

    Means.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>% 
      dplyr::select(Distance, Mean_Pop)
    SDs.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>% 
      dplyr::select(Distance, SD_Pop)
    Medians.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>% 
      dplyr::select(Distance, Median_Pop)
    Quantiles.Full <- rbind(Quantiles.10, Quantiles.50, Quantiles.100) %>%
      as.data.frame() 
    Quantiles.Full$Quantile <- rownames(Quantiles.Full)
    
      Means[,YEAR] <- Means.Full$Mean_Pop
      
      SDs[ ,YEAR] <- SDs.Full$SD_Pop
      
      Medians[ ,YEAR] <- Medians.Full$Median_Pop
      
      Quantiles[, YEAR] <- Quantiles.Full$. 
      Quantiles$Quantile <- Quantiles.Full$Quantile
      
    
  }
   Means <- as.data.frame(Means) %>% 
     mutate(Scenario = Names[S]) %>% 
     mutate(Distances = Dist.Names) %>% 
     pivot_longer(cols=-c("Scenario", "Distances"), values_to = "Mean_Abundance", names_to = "Year") %>% 
     mutate(Year = rep(seq(1960,2018,1), times = 3))
   
   SDs <- as.data.frame(SDs) %>% 
     mutate(Scenario = Names[S]) %>% 
     mutate(Distances = Dist.Names)%>% 
     pivot_longer(cols=-c("Scenario", "Distances"), values_to = "SD_Abundance", names_to = "Year") %>% 
     mutate(Year = rep(seq(1960,2018,1), times = 3))
   
   Medians <- as.data.frame(Medians) %>% 
     mutate(Scenario = Names[S]) %>% 
     mutate(Distances = Dist.Names)%>% 
     pivot_longer(cols=-c("Scenario", "Distances"), values_to = "Median_Abundance", names_to = "Year") %>% 
     mutate(Year = rep(seq(1960,2018,1), times = 3))
   
   Quantiles <- as.data.frame(Quantiles) %>%
     mutate(Scenario = rep(Names[S], 6)) %>%
     mutate(Distances = rep(Dist.Names, each=2)) %>%
     pivot_longer(cols=-c("Scenario", "Distances", "Quantile"), values_to = "Quantile_Abundance", names_to = "Year") %>%
     mutate(Year = rep(seq(1960,2018,1), times = 6)) %>% 
     mutate(Quantile = str_replace(Quantile, "%[[:digit:]]$", "%")) %>% 
     pivot_wider(names_from = "Quantile", values_from="Quantile_Abundance")

   Pop.Dist.Mean[[S]] <- Means
   Pop.Dist.SD[[S]] <- SDs
   Pop.Dist.Median[[S]] <- Medians
   Pop.Dist.Quant[[S]] <- Quantiles
  
}


S00.means <- Pop.Dist.Mean[[1]]
S00.medians <- Pop.Dist.Median[[1]]
S00.SD <- Pop.Dist.SD[[1]]
S00.Quant <- Pop.Dist.Quant[[1]]
S00.data <- cbind(S00.means, S00.SD$SD_Abundance, S00.medians$Median_Abundance, S00.Quant$`2.5%`, S00.Quant$`97.5%`) %>% 
  rename(SD_Abundance = "S00.SD$SD_Abundance",
         Median_Abundance = "S00.medians$Median_Abundance",
         Q2.5 = "S00.Quant$`2.5%`",
         Q97.5 = "S00.Quant$`97.5%`")


S00.abundance.dist.plot <- ggplot()+
  geom_line(data=S00.data, aes(x=Year, y=Median_Abundance, group=Distances, linetype=Distances),col="#36753B")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S00.data, aes(x=Year, y=Median_Abundance, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#36753B", alpha=0.2)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
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
  ggplot2::annotate("text", x=1961, y=7700, label="(a)", size = 2.5, fontface=1, hjust=0)
S00.abundance.dist.plot

S01.means <- Pop.Dist.Mean[[2]]
S01.SD <- Pop.Dist.SD[[2]]
S01.medians <- Pop.Dist.Median[[2]]
S01.Quant <- Pop.Dist.Quant[[2]]
S01.data <- cbind(S01.means, S01.SD$SD_Abundance, S01.medians$Median_Abundance, S01.Quant$`2.5%`, S01.Quant$`97.5%`) %>% 
  rename(SD_Abundance = "S01.SD$SD_Abundance",
         Median_Abundance = "S01.medians$Median_Abundance",
         Q2.5 = "S01.Quant$`2.5%`",
         Q97.5 = "S01.Quant$`97.5%`")

S01.abundance.dist.plot <- ggplot()+
  geom_line(data=S01.data, aes(x=Year, y=Median_Abundance, group=Distances, linetype=Distances),col="#302383")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S01.data, aes(x=Year, y=Median_Abundance, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#302383", alpha=0.2)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
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
  ggplot2::annotate("text", x=1960.5, y=7700, label="(b)", size = 2.5, fontface=1, hjust=0)
S01.abundance.dist.plot

S02.means <- Pop.Dist.Mean[[3]]
S02.SD <- Pop.Dist.SD[[3]]
S02.medians <- Pop.Dist.Median[[3]]
S02.Quant <- Pop.Dist.Quant[[3]]
S02.data <- cbind(S02.means, S02.SD$SD_Abundance, S02.medians$Median_Abundance, S02.Quant$`2.5%`, S02.Quant$`97.5%`) %>% 
  rename(SD_Abundance = "S02.SD$SD_Abundance",
         Median_Abundance = "S02.medians$Median_Abundance",
         Q2.5 = "S02.Quant$`2.5%`",
         Q97.5 = "S02.Quant$`97.5%`")

S02.abundance.dist.plot <- ggplot()+
  geom_line(data=S02.data, aes(x=Year, y=Median_Abundance, group=Distances, linetype=Distances),col="#66CCEE")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S02.data, aes(x=Year, y=Median_Abundance, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#66CCEE", alpha=0.1)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
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
  ggplot2::annotate("text", x=1960.5, y=7600, label="(c)", size = 2.5, fontface=1, hjust=0)
S02.abundance.dist.plot


S03.means <- Pop.Dist.Mean[[4]]
S03.SD <- Pop.Dist.SD[[4]]
S03.medians <- Pop.Dist.Median[[4]]
S03.Quant <- Pop.Dist.Quant[[4]]
S03.data <- cbind(S03.means, S03.SD$SD_Abundance, S03.medians$Median_Abundance, S03.Quant$`2.5%`, S03.Quant$`97.5%`) %>% 
  rename(SD_Abundance = "S03.SD$SD_Abundance",
         Median_Abundance = "S03.medians$Median_Abundance",
         Q2.5 = "S03.Quant$`2.5%`",
         Q97.5 = "S03.Quant$`97.5%`")

S03.abundance.dist.plot <- ggplot()+
  geom_line(data=S03.data, aes(x=Year, y=Median_Abundance, group=Distances, linetype=Distances),col="#BBCC33")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"))+
  geom_ribbon(data=S03.data, aes(x=Year, y=Median_Abundance, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#BBCC33", alpha=0.1)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
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
  ggplot2::annotate("text", x=1960.5, y=8300, label="(d)", size = 2.5, fontface=1, hjust=0)
S03.abundance.dist.plot

legend.plot <- ggplot()+
  geom_line(data=S03.data, aes(x=Year, y=Mean_Abundance, group=Distances, linetype=Distances),col="grey20")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  theme_classic()

## Put it all together
setwd(fig_dir)
x.label <- textGrob("Year", gp=gpar(fontsize=9))
y.label <- textGrob("Median abundance", gp=gpar(fontsize=9), rot=90)
legend <- gtable_filter(ggplotGrob(legend.plot), "guide-box")

AbundancexDistance <-grid.arrange(arrangeGrob(S00.abundance.dist.plot + theme(legend.position="none"),
                                              S01.abundance.dist.plot + theme(legend.position="none"),
                                              S02.abundance.dist.plot + theme(legend.position="none"),
                                              S03.abundance.dist.plot + theme(legend.position="none"),
                                              left=y.label,
                                              bottom=x.label,
                                              right=legend))

ggsave(AbundancexDistance, filename="Abundance_Distance.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )


#* Catch ####
setwd(pop_dir)

Pop.Catch <- list()

# Each layer is a simulation, rows are cells and columns are years
Pop.Catch[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S00"))
Pop.Catch[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S01"))
Pop.Catch[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S02"))
Pop.Catch[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S03"))

Pop.Catch.Mean <- list()
Pop.Catch.SD <- list()
Pop.Catch.Median <- list()
Pop.Catch.Quant <- list()

Catch.Quantiles.10 <- array(0, dim=c(nrow(Whole_Pop_Ages), 2))
Catch.Quantiles.50 <- array(0, dim=c(nrow(Whole_Pop_Ages), 2))
Catch.Quantiles.100 <- array(0, dim=c(nrow(Whole_Pop_Ages), 2))

# Ths is just the number of cells in 
Names <- c("Historical and Current NTZs", "Neither NTZs nor Temporal Management", 
           "Temporal and Spatial Management","Temporal Management Only" )

for(S in 1:4){
  
  Scenario <- Pop.Catch[[S]]
  
  Means <- array(0, dim=c(3, 59)) %>% 
    as.data.frame(.)
  SDs <- array(0, dim=c(3, 59))%>% 
    as.data.frame(.)
  Medians <- array(0, dim=c(3, 59))%>% 
    as.data.frame(.)
  Quantiles <- array(0, dim=c(6, 59))%>% 
    as.data.frame(.)
  
  for(SIM in 1:100){
    temp[,,SIM] <- Scenario[[SIM]]
  }
  for(YEAR in 1:59){
    
    temp1 <- temp[,YEAR,]
    
    temp2 <- as.data.frame(temp1) %>% 
      mutate(CellID = row_number()) %>% 
      mutate(Fished_17 = water$Fished_2017) 
    
    Dist.10km <- temp2 %>% 
      filter(CellID %in% c(Distances[[1]])) %>% 
      {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      mutate(Mean_Pop = rowMeans(.[1:100])) %>% 
      mutate(SD_Pop = rowSds(as.matrix(.[,1:100]))) %>% 
      mutate(Median_Pop = rowMedians(as.matrix(.[1:100]))) %>%
      mutate(Distance = "10km") 
    
    Dist.50km <- temp2 %>% 
      filter(CellID %in% c(Distances[[2]])) %>% 
      {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      mutate(Mean_Pop = rowMeans(.[1:100])) %>% 
      mutate(SD_Pop = rowSds(as.matrix(.[,1:100]))) %>% 
      mutate(Median_Pop = rowMedians(as.matrix(.[1:100]))) %>%
      mutate(Distance = "50km")
    
    Dist.100km <- temp2 %>% 
      filter(CellID %in% c(Distances[[3]])) %>% 
      {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      mutate(Mean_Pop = rowMeans(.[1:100])) %>% 
      mutate(SD_Pop = rowSds(as.matrix(.[,1:100]))) %>% 
      mutate(Median_Pop = rowMedians(as.matrix(.[1:100]))) %>%
      mutate(Distance = "100km") 
    
    temp2 <- as.matrix(Dist.10km[ , 1:100])
    
    Catch.Quantiles.10 <- quantile(temp2, probs=c(0.025, 0.975)) %>% 
      as.data.frame() %>% 
      mutate(Distance = "10 km") 
    
    temp2 <- as.matrix(Dist.50km[ , 1:100])
    
    Catch.Quantiles.50 <-quantile(temp2, probs=c(0.025, 0.975)) %>% 
      as.data.frame() %>% 
      mutate(Distance = "50 km") 
    
    temp2 <- as.matrix(Dist.100km[ , 1:100])
    
    Catch.Quantiles.100 <-quantile(temp2, probs=c(0.025, 0.975)) %>% 
      as.data.frame() %>% 
      mutate(Distance = "100 km") 
    
    Means.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>% 
      dplyr::select(Distance, Mean_Pop)
    SDs.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>% 
      dplyr::select(Distance, SD_Pop)
    Medians.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>% 
      dplyr::select(Distance, Median_Pop)
    Quantiles.Full <- rbind(Catch.Quantiles.10, Catch.Quantiles.50, Catch.Quantiles.100) %>%
      as.data.frame() 
    Quantiles.Full$Quantile <- rownames(Quantiles.Full)
    
    Means[,YEAR] <- Means.Full$Mean_Pop
    
    SDs[ ,YEAR] <- SDs.Full$SD_Pop
    
    Medians[ ,YEAR] <- Medians.Full$Median_Pop
    
    Quantiles[, YEAR] <- Quantiles.Full$. 
    Quantiles$Quantile <- Quantiles.Full$Quantile
    
    
  }
  Means <- as.data.frame(Means) %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Distances = Dist.Names) %>% 
    pivot_longer(cols=-c("Scenario", "Distances"), values_to = "Mean_Catch", names_to = "Year") %>% 
    mutate(Year = rep(seq(1960,2018,1), times = 3))
  
  SDs <- as.data.frame(SDs) %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Distances = Dist.Names)%>% 
    pivot_longer(cols=-c("Scenario", "Distances"), values_to = "SD_Catch", names_to = "Year") %>% 
    mutate(Year = rep(seq(1960,2018,1), times = 3))
  
  Medians <- as.data.frame(Medians) %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Distances = Dist.Names)%>% 
    pivot_longer(cols=-c("Scenario", "Distances"), values_to = "Median_Catch", names_to = "Year") %>% 
    mutate(Year = rep(seq(1960,2018,1), times = 3))
  
  Quantiles <- as.data.frame(Quantiles) %>%
    mutate(Scenario = rep(Names[S], 6)) %>%
    mutate(Distances = rep(Dist.Names, each=2)) %>%
    pivot_longer(cols=-c("Scenario", "Distances", "Quantile"), values_to = "Quantile_Catch", names_to = "Year") %>%
    mutate(Year = rep(seq(1960,2018,1), times = 6)) %>% 
    mutate(Quantile = str_replace(Quantile, "%[[:digit:]]$", "%")) %>% 
    pivot_wider(names_from = "Quantile", values_from="Quantile_Catch")
  
  Pop.Catch.Mean[[S]] <- Means
  Pop.Catch.SD[[S]] <- SDs
  Pop.Catch.Median[[S]] <- Medians
  Pop.Catch.Quant[[S]] <- Quantiles
  
}

S00.means <- Pop.Catch.Mean[[1]]
S00.medians <- Pop.Catch.Median[[1]]
S00.SD <- Pop.Catch.SD[[1]]
S00.Quant <- Pop.Catch.Quant[[1]]
S00.data <- cbind(S00.means, S00.SD$SD_Catch, S00.medians$Median_Catch, S00.Quant$`2.5%`, S00.Quant$`97.5%`) %>% 
  rename(SD_Catch = "S00.SD$SD_Catch",
         Median_Catch = "S00.medians$Median_Catch",
         Q2.5 = "S00.Quant$`2.5%`",
         Q97.5 = "S00.Quant$`97.5%`")


S00.catch.dist.plot <- ggplot()+
  geom_line(data=S00.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="#36753B")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S00.data, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#36753B", alpha=0.2)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
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
  ggplot2::annotate("text", x=1961, y=700, label="(a)", size = 2.5, fontface=1, hjust=0)
S00.catch.dist.plot


S01.means <- Pop.Catch.Mean[[2]]
S01.medians <- Pop.Catch.Median[[2]]
S01.SD <- Pop.Catch.SD[[2]]
S01.Quant <- Pop.Catch.Quant[[2]]
S01.data <- cbind(S01.means, S01.SD$SD_Catch, S01.medians$Median_Catch, S01.Quant$`2.5%`, S01.Quant$`97.5%`) %>% 
  rename(SD_Catch = "S01.SD$SD_Catch",
         Median_Catch = "S01.medians$Median_Catch",
         Q2.5 = "S01.Quant$`2.5%`",
         Q97.5 = "S01.Quant$`97.5%`")

S01.catch.dist.plot <- ggplot()+
  geom_line(data=S01.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="#302383")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S01.data, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#302383", alpha=0.2)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
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
  ggplot2::annotate("text", x=1961, y=750, label="(b)", size = 2.5, fontface=1, hjust=0)
S01.catch.dist.plot

S02.means <- Pop.Catch.Mean[[3]]
S02.medians <- Pop.Catch.Median[[3]]
S02.SD <- Pop.Catch.SD[[3]]
S02.Quant <- Pop.Catch.Quant[[3]]
S02.data <- cbind(S02.means, S02.SD$SD_Catch, S02.medians$Median_Catch, S02.Quant$`2.5%`, S02.Quant$`97.5%`) %>% 
  rename(SD_Catch = "S02.SD$SD_Catch",
         Median_Catch = "S02.medians$Median_Catch",
         Q2.5 = "S02.Quant$`2.5%`",
         Q97.5 = "S02.Quant$`97.5%`")

S02.catch.dist.plot <- ggplot()+
  geom_line(data=S02.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="#66CCEE")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S02.data, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#66CCEE", alpha=0.1)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
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
  ggplot2::annotate("text", x=1961, y=700, label="(c)", size = 2.5, fontface=1, hjust=0)
S02.catch.dist.plot


S03.means <- Pop.Catch.Mean[[4]]
S03.medians <- Pop.Catch.Median[[4]]
S03.SD <- Pop.Catch.SD[[4]]
S03.Quant <- Pop.Catch.Quant[[4]]
S03.data <- cbind(S03.means, S03.SD$SD_Catch, S03.medians$Median_Catch, S03.Quant$`2.5%`, S03.Quant$`97.5%`) %>% 
  rename(SD_Catch = "S03.SD$SD_Catch",
         Median_Catch = "S03.medians$Median_Catch",
         Q2.5 = "S03.Quant$`2.5%`",
         Q97.5 = "S03.Quant$`97.5%`")

S03.catch.dist.plot <- ggplot()+
  geom_line(data=S03.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="#BBCC33")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"))+
  geom_ribbon(data=S03.data, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#BBCC33", alpha=0.1)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
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
  ggplot2::annotate("text", x=1960.5, y=750, label="(d)", size = 2.5, fontface=1, hjust=0)
S03.catch.dist.plot

legend.plot <- ggplot()+
  geom_line(data=S03.data, aes(x=Year, y=Mean_Catch, group=Distances, linetype=Distances),col="grey20")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  theme_classic()

## Put it all together
setwd(fig_dir)
x.label <- textGrob("Year", gp=gpar(fontsize=9))
y.label <- textGrob("Median catch", gp=gpar(fontsize=9), rot=90)
legend <- gtable_filter(ggplotGrob(legend.plot), "guide-box")

CatchxDistance <-grid.arrange(arrangeGrob(S00.catch.dist.plot + theme(legend.position="none"),
                                              S01.catch.dist.plot + theme(legend.position="none"),
                                              S02.catch.dist.plot + theme(legend.position="none"),
                                              S03.catch.dist.plot + theme(legend.position="none"),
                                              left=y.label,
                                              bottom=x.label,
                                              right=legend))

ggsave(CatchxDistance, filename="Catch_Distance.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )


#* Combining scenarios ####

abundance.10_50 <- NULL
catch.10_50 <- NULL

for (S in 1:4){
  temp.mean <- Pop.Dist.Mean[[S]]
  temp.SD <- Pop.Dist.SD[[S]]
  temp.median <- Pop.Dist.Median[[S]]
  temp.quant <- Pop.Dist.Quant[[S]]
  
  temp.all <- cbind(temp.mean, temp.SD$SD_Abundance, temp.median$Median_Abundance, temp.quant$`2.5%`, temp.quant$`97.5%`) %>% 
    rename(SD_Abundance = "temp.SD$SD_Abundance",
           Median_Abundance = "temp.median$Median_Abundance",
           Q2.5 = "temp.quant$`2.5%`",
           Q97.5 = "temp.quant$`97.5%`") %>% 
    filter(Distances %in% c("0-10 km"))
  
  abundance.10_50 <- rbind(abundance.10_50, temp.all)
  
  temp.mean <- Pop.Catch.Mean[[S]]
  temp.SD <- Pop.Catch.SD[[S]]
  temp.median <- Pop.Catch.Median[[S]]
  temp.quant <- Pop.Catch.Quant[[S]]
  
  temp.all <- cbind(temp.mean, temp.SD$SD_Catch, temp.median$Median_Catch, temp.quant$`2.5%`, temp.quant$`97.5%`) %>% 
    rename(SD_Catch = "temp.SD$SD_Catch",
           Median_Catch = "temp.median$Median_Catch",
           Q2.5 = "temp.quant$`2.5%`",
           Q97.5 = "temp.quant$`97.5%`") %>% 
    filter(Distances %in% c("0-10 km"))
  
  catch.10_50 <- rbind(catch.10_50, temp.all)
}

abundance.Pre_1987 <- abundance.10_50 %>% 
  filter(Year<=1986) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs") & Year>1985, "Historical and\ncurrent NTZs", 
                                                             ifelse(Scenario %in% c("Temporal Management Only") & Year>1985, "Temporal\nmanagement only", 
                                                                    ifelse(Scenario %in% c("No NTZs or Temporal Management"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))

abundance.10_50.plot <- abundance.10_50 %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs"), "Historical and\ncurrent NTZs", 
                                                                 ifelse(Scenario %in% c("Temporal Management Only"), "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("Neither NTZs nor Temporal Management"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(SE_Abundance= SD_Abundance/sqrt(100)) %>% 
  mutate(CI = SE_Abundance * t.score) %>% 
  filter(Year>1985) %>% 
  ggplot(.)+
  geom_line(aes(x=Year, y=Median_Abundance, group=Scenario, colour=ColourGroup), size=0.7)+
  geom_ribbon(aes(x=Year, y=Median_Abundance, ymin=Q2.5, ymax=Q97.5, fill=ColourGroup, group=Scenario), alpha=0.2)+
 
  scale_fill_manual(values= c("Historical and\ncurrent NTZs"="#36753B", "No NTZs or\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("NTZs and\ntemporal management"="#66CCEE","Historical and\ncurrent NTZs"="#36753B", "Temporal\nmanagement only"="#BBCC33", 
                                 "Neither NTZs nor\ntemporal management"="#302383"), breaks= c("NTZs and\ntemporal management", "Historical and\ncurrent NTZs", "Temporal\nmanagement only", "Neither NTZs nor\ntemporal management"),
                      name= "Spatial and temporal\nmanagement scenario")+ 
  geom_line(data=abundance.Pre_1987, aes(x=Year, y=Median_Abundance, group=Scenario, colour="grey20"), size=0.7)+
  geom_ribbon(data=abundance.Pre_1987, aes(x=Year, y=Median_Abundance, ymin=Q2.5, ymax=Q97.5, fill="grey20", group=Scenario), alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab("Median abundance within\n0-10km of boat ramp")+
  xlim(1987,2020)+
  ylim(0, 500)+
  scale_linetype_manual(values = c("longdash", "solid" ), labels=c("Always Fished", "NTZ Area"), name="Model Area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8),
        axis.title = element_text(size=9))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1987, y=500, label="(a)", size = 2.5, fontface=1)
abundance.10_50.plot


Catch.Pre_1987 <- catch.10_50 %>% 
  filter(Year<=1986) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs") & Year>1985, "Historical and\ncurrent NTZs", 
                                                             ifelse(Scenario %in% c("Temporal Management Only") & Year>1985, "Temporal\nmanagement only", 
                                                                    ifelse(Scenario %in% c("Neither NTZs nor Temporal Management"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))
Catch.10_50.plot <- catch.10_50 %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs"), "Historical and\ncurrent NTZs", 
                                                            ifelse(Scenario %in% c("Temporal Management Only"), "Temporal\nmanagement only", 
                                                                   ifelse(Scenario %in% c("Neither NTZs nor Temporal Management"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  filter(Year>1985) %>% 
  # mutate(SE_Catch= SD_Catch/sqrt(100)) %>% 
  # mutate(CI = SE_Catch * t.score) %>% 
  ggplot(.)+
  geom_line(aes(x=Year, y=Median_Catch, group=Scenario, colour=ColourGroup), size=0.7)+
  geom_ribbon(aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, fill=ColourGroup, group=Scenario), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent NTZs"="#36753B", "Neither NTZs nor\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("NTZs and\ntemporal management"="#66CCEE","Historical and\ncurrent NTZs"="#36753B", "Temporal\nmanagement only"="#BBCC33", 
                                 "Neither NTZs nor\ntemporal management"="#302383"), breaks= c("NTZs and\ntemporal management", "Historical and\ncurrent NTZs", "Temporal\nmanagement only", "Neither NTZs nor\ntemporal management"),
                      name= "Spatial and temporal\nmanagement scenario")+ 
  geom_line(data=Catch.Pre_1987, aes(x=Year, y=Median_Catch, group=Scenario, colour="grey20"), size=0.7)+
  geom_ribbon(data=Catch.Pre_1987, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, fill="grey20", group=Scenario), alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab("Median catch within\n0-10km of boat ramp")+
  xlim(1987,2020)+
  ylim(0,60)+
  scale_linetype_manual(values = c("longdash", "solid" ), labels=c("Always Fished", "NTZ Area"), name="Model Area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8),
        axis.title = element_text(size=9))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1987, y=60, label="(b)", size = 2.5)
Catch.10_50.plot

## Put it together
setwd(fig_dir)
x.label <- textGrob("Year", gp=gpar(fontsize=9))
legend <- gtable_filter(ggplotGrob(Catch.10_50.plot), "guide-box")

Catch.AbundancexDistance <-grid.arrange(arrangeGrob(abundance.10_50.plot + theme(legend.position="none"),
                                          Catch.10_50.plot + theme(legend.position="none"),
                                          bottom=x.label,
                                          right=legend))
ggsave(Catch.AbundancexDistance, filename="Distance_Combined.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

## CPUE for 0-10 km from Boat Ramp 
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

#* Sum up effort for every year in each of the scenarios

Boat_Days_sum <- NULL

for(S in 1:4){
  temp <- Boat_Days[[S]]
  temp2 <- temp[as.numeric(Distances[[1]]), ,] %>% 
    colSums(., dim=2) %>% 
    as.data.frame() %>% 
    mutate(Scenario = Names[S])
  
  Boat_Days_sum <- rbind(Boat_Days_sum, temp2)
}

Boat_Days_sum <- Boat_Days_sum %>% 
  rename(Effort = ".")

cpue.0_10 <- catch.10_50 %>% 
  mutate(Effort = Boat_Days_sum$Effort) %>% 
  mutate(Median_CPUE = Median_Catch/Effort,
         CPUE_2.5 = Q2.5/Effort,
         CPUE_97.5 = Q97.5/Effort)


CPUE.Pre_1987 <- cpue.0_10 %>% 
  filter(Year<=1986) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs") & Year>1985, "Historical and\ncurrent NTZs", 
                                                             ifelse(Scenario %in% c("Temporal Management Only") & Year>1985, "Temporal\nmanagement only", 
                                                                    ifelse(Scenario %in% c("Neither NTZs nor Temporal Management"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))
CPUE.0_10.plot <- cpue.0_10 %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs"), "Historical and\ncurrent NTZs", 
                                                             ifelse(Scenario %in% c("Temporal Management Only"), "Temporal\nmanagement only", 
                                                                    ifelse(Scenario %in% c("Neither NTZs nor Temporal Management"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  filter(Year>1985) %>% 
  # mutate(SE_Catch= SD_Catch/sqrt(100)) %>% 
  # mutate(CI = SE_Catch * t.score) %>% 
  ggplot(.)+
  geom_line(aes(x=Year, y=Median_CPUE, group=Scenario, colour=ColourGroup), size=0.7)+
  geom_ribbon(aes(x=Year, y=Median_CPUE, ymin=CPUE_2.5, ymax=CPUE_97.5, fill=ColourGroup, group=Scenario), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent NTZs"="#36753B", "Neither NTZs nor\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("NTZs and\ntemporal management"="#66CCEE","Historical and\ncurrent NTZs"="#36753B", "Temporal\nmanagement only"="#BBCC33", 
                                 "Neither NTZs nor\ntemporal management"="#302383"), breaks= c("NTZs and\ntemporal management", "Historical and\ncurrent NTZs", "Temporal\nmanagement only", "Neither NTZs nor\ntemporal management"),
                      name= "Spatial and temporal\nmanagement scenario")+ 
  geom_line(data=CPUE.Pre_1987, aes(x=Year, y=Median_CPUE, group=Scenario, colour="grey20"), size=0.7)+
  geom_ribbon(data=CPUE.Pre_1987, aes(x=Year, y=Median_CPUE, ymin=CPUE_2.5, ymax=CPUE_97.5, fill="grey20", group=Scenario), alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab("CPUE within\n0-10km of boat ramp")+
  xlim(1987,2020)+
  ylim(0,0.08)+
  scale_linetype_manual(values = c("longdash", "solid" ), labels=c("Always Fished", "NTZ Area"), name="Model Area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8),
        axis.title = element_text(size=9))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1987, y=0.08, label="(b)", size = 2.5)
CPUE.0_10.plot

#### CATCH PER UNIT EFFORT PLOTS ####
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

#* Sum up effort for every year in each of the scenarios

Boat_Days_sum <- NULL

for(S in 1:4){
  temp <- Boat_Days[[S]]
  temp2 <- colSums(temp, dim=2) %>% 
    as.data.frame() %>% 
    mutate(Scenario = Names[S])
  
  Boat_Days_sum <- rbind(Boat_Days_sum, temp2)
}

Boat_Days_sum <- Boat_Days_sum %>% 
  rename(Effort = ".")

#* Format catch data and add to effort ####
setwd(pop_dir)

Pop.Catch <- list()

# Each layer is a simulation, rows are cells and columns are years
Pop.Catch[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S00"))
Pop.Catch[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S01"))
Pop.Catch[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S02"))
Pop.Catch[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S03"))

temp <- array(0, dim=c(NCELL,59,100))
Pop.Catch.Mean <- NULL
Pop.Catch.SD <- NULL

Names <- c("Historical and Current Management", "No Spatial Management", 
           "Temporal and Spatial Management","Temporal Management Only" )

for(S in 1:4){
  
  Scenario <- Pop.Catch[[S]]
  
  Means <- array(0, dim=c(1, 59)) %>% 
    as.data.frame(.)
  SDs <- array(0, dim=c(1, 59))%>% 
    as.data.frame(.)
  
  for(SIM in 1:100){
    temp[,,SIM] <- Scenario[[SIM]]
  }
  for(YEAR in 1:59){
    
    temp1 <- temp[,YEAR,]
    
    temp2 <- temp1 %>% 
      as.data.frame() %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      mutate(Mean_Catch = rowMeans(.[1:100])) %>% 
      mutate(SD_Catch = rowSds(as.matrix(.[,1:100]))) 
    
    Means[,YEAR] <- temp2$Mean_Catch
    SDs[,YEAR] <- temp2$SD_Catch
  }
  Means <- as.data.frame(Means) %>% 
    mutate(Scenario = Names[S]) %>% 
    pivot_longer(cols=-c("Scenario"), values_to = "Mean_Catch", names_to = "Year") %>% 
    mutate(Year = seq(1960,2018,1))
  
  SDs <- as.data.frame(SDs) %>% 
    mutate(Scenario = Names[S]) %>% 
    pivot_longer(cols=-c("Scenario"), values_to = "SD_Catch", names_to = "Year") %>% 
    mutate(Year = seq(1960,2018,1))
  
  Pop.Catch.Mean <- rbind(Pop.Catch.Mean, Means)
  Pop.Catch.SD <- rbind(Pop.Catch.SD, SDs)
  
}

CPUE <- Boat_Days_sum %>% 
  mutate(Catch = Pop.Catch.Mean$Mean_Catch) %>% 
  mutate(SD.Catch = Pop.Catch.SD$SD_Catch) %>% 
  mutate(Year = rep(seq(1960,2018), 4)) %>% 
  mutate(cpue = Catch/Effort) %>% 
  mutate(cpue = ifelse(is.nan(cpue), 0, cpue))

#* Plot CPUE for each scenario ####

CPUE.Pre_1987 <- CPUE %>% 
  filter(Year<=1986) %>%  
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs") & Year>1985, "Historical and\ncurrent management", 
                                                             ifelse(Scenario %in% c("Temporal Management Only") & Year>1985, "Temporal\nmanagement only", 
                                                                    ifelse(Scenario %in% c("No NTZs or Temporal Management"), "No NTZs or\ntemporal management", "NTZs and\ntemporal management")))))

CPUE.plot <- CPUE %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs"), "Historical and\ncurrent management", 
                                                             ifelse(Scenario %in% c("Temporal Management Only"), "Temporal\nmanagement only", 
                                                                    ifelse(Scenario %in% c("No NTZs or Temporal Management"), "No NTZs or\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  filter(Year>1985) %>% 
  ggplot(.)+
  geom_line(aes(x=Year, y=cpue, group=Scenario, colour=ColourGroup), size=0.7)+
  scale_fill_manual(values= c("Historical and\ncurrent management"="#36753B", "No NTZs or\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c( "Historical and\ncurrent management"="#36753B", "No NTZs or\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                                  "Temporal\nmanagement only"="#BBCC33"), name= "Spatial and Temporal\nManagement Scenario")+ 
  geom_line(data=CPUE.Pre_1987, aes(x=Year, y=cpue, group=Scenario, colour="grey20"), size=0.7)+
  theme_classic()+
  xlab(NULL)+
  ylab("CPUE")+
  xlim(1970,2020)+
  ylim(0,0.3)+
  #scale_linetype_manual(values = c("longdash", "solid" ), labels=c("Always Fished", "NTZ Area"), name="Model Area")+
  theme(legend.title = element_text(size=9, face="bold"), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8),
        axis.title = element_text(size=9))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")
  #ggplot2::annotate("text", x=1960, y=400, label="(b)", size = 2.5, fontface=2)
CPUE.plot


#### PLOT JUST ONE SCENARIO/SIMULATION FOR CHECKING ####
NTZ_Ages_S00 <- NULL
F_Ages_S00 <- NULL

for(SIM in 1:length(SP_Pop_NTZ_S00)){
  
  temp <- as.data.frame(colSums(SIM.SP.NTZ[[1]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S00 <- cbind(NTZ_Ages_S00, temp$Number)
  
  temp <- as.data.frame(colSums(SIM.Sp.F[[1]])) %>% 
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
  mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>%
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
  mutate(Mean_Pop = rowMeans(.[,3:ncol(.)])) %>%
  mutate(SD_Pop = rowSds(as.matrix(.[,3:ncol(.)]))) %>%
  mutate(Scenario = "S00") %>% 
  dplyr::select(Mean_Pop, SD_Pop,Scenario, Stage, Mod_Year)

Whole_Pop_Ages_NTZ <- NTZ_Ages_S00 %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- F_Ages_S00 %>% 
  mutate(Zone = "F")

Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) 


Pre_1987_NTZ <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Stage %in% c("Legal"))

Pre_1987_F <- Whole_Pop_Ages %>% 
  filter(Mod_Year<1987) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Stage %in% c("Legal")) 
# mutate(Mean_Pop = ifelse(Mod_Year==1960, Mean_Pop*40, Mean_Pop),
#        SD_Pop = ifelse(Mod_Year==1960, SD_Pop*(10^16), SD_Pop))

Recruit_F <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("F")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent NTZs", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "No NTZs or\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup)) 

line.recruit <- Whole_Pop_Ages %>% 
  filter(Stage %in% c("Legal")) %>% 
  filter(Zone %in% c("NTZ")) %>% 
  filter(Mod_Year >=1986) %>% 
  mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("S00") & Mod_Year>1985, "Historical and\ncurrent NTZs", 
                                                                 ifelse(Scenario %in% c("S03") & Mod_Year>1985, "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("S01"), "No NTZs or\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  mutate(ColourGroup = as.factor(ColourGroup))%>% 
  ggplot(.)+
  geom_line(aes(x=Mod_Year, y=Mean_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  geom_line(data=Recruit_F, aes(x=Mod_Year, y=Mean_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
  geom_ribbon(data=Recruit_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent NTZs"="#36753B", "No NTZs or\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("Pre-1987"="grey20", "Historical and\ncurrent NTZs"="#36753B", "No NTZs or\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                                 "Temporal\nmanagement only"="#BBCC33"), name= "Spatial and temporal\nmanagement scenario")+ 
  geom_line(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_NTZ, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20",alpha=0.2)+
  geom_line(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, color="grey20", group=Scenario, linetype=Zone), size=0.7)+
  geom_ribbon(data=Pre_1987_F, aes(x=Mod_Year, y=Mean_Pop, ymin=Mean_Pop-SD_Pop, ymax=Mean_Pop+SD_Pop, group=Scenario), fill="grey20", alpha=0.2)+
  theme_classic()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(1960,2020)+
  scale_linetype_manual(values = c("longdash", "solid" ), labels=c("Always fished", "NTZ area"), name="Model area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1964.5, y=7.3, label="(a) <1 year old", size = 2.5, fontface=1)
line.recruit

#### MAP OF THE NINGALOO MODEL ####
water <- water %>% 
  mutate(WHA = ifelse(ID %in% c(model_WHA$row.id), "Y", "N")) %>% 
  mutate(Map_Colour = ifelse(Fished_2017=="N", "NTZ", ifelse(WHA=="Y", "WHA", "None")))
  



map <- water %>% 
  filter(!ID==387) %>% # Weird cell on the side that makes things look odd
  ggplot(.)+
  geom_sf(aes(fill=Map_Colour), colour="grey20", lwd=0.2)+
  scale_fill_manual(values=c("NTZ"= "#48A02D", "WHA"="#D6CF7D", "None"="skyblue1"),
                    labels = c("No-take zone", "World heritage area", "Outside world heritage\nand marine park area"),
                    name="Zone type")+
  theme_void() +
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(1.5,"line")) 
map

setwd(fig_dir)
ggsave(map, filename="Whole_map.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )



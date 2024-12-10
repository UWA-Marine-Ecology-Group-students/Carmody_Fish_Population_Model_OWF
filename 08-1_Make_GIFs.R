#######################################################

# Script for creating little GIFs of the catch and
# abundance over time 

#######################################################
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
gif_dir <-  paste(working.dir, "GIFs", sep="/")

a4.width = 160
model.name <- "ningaloo"
## Functions
setwd(working.dir)
source("X_Functions.R")

Names <- c("Current NTZs", "No temporal management or NTZs", 
           "Temporal management and NTZs","Temporal management")
#### READ IN DATA ####
setwd(sg_dir)

weight <- readRDS("Weight")
mature <- readRDS("maturity")

NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
water <- readRDS(paste0(model.name, sep="_","water")) %>% 
  mutate(cell_area = as.vector(cell_area/1000000))
NCELL <- nrow(water)

setwd(sp_dir)
bathy <- raster("ga_bathy_ningaloocrop.tif")
WHA <- st_read("2013_02_WorldHeritageMarineProgramme.shp") %>% 
  st_transform(4283)%>%
  st_make_valid %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) 

old_MP <- st_read("1987_NMP_boundary.shp") %>% 
  st_transform(4283)%>%
  st_make_valid %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) %>% 
  dplyr::filter(DESCRIPTOR == "sanctuary") %>% # Only pick the areas that are sanctuaries 
  mutate(Year.Sanct = 1987) %>% 
  rename(Name = "NAME") %>% 
  dplyr::select(Name, Year.Sanct, geometry) 

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
NTZ <- rbind(old_MP, NTZ, AMP_NTZ) # Put all the different NTZs together into one object
NTZ <- st_make_valid(NTZ) 

NTZ_cropped <- NTZ %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-23.45, ymax=-20.5)
plot(NTZ_cropped$geometry)

setwd(sg_dir)
NTZ_ID <- readRDS(paste0(model.name, sep="_", "NTZ_Cell_ID"))
F_ID <- readRDS(paste0(model.name, sep="_", "F_Cell_ID"))
Shallow_ID <- c(NTZ_ID, F_ID)

# crop to the WHA area as well
water_points <- st_centroid_within_poly(water) 

water_bathy <- raster::extract(bathy, water_points, fun=mean, df=TRUE)

water_bathy <- water_bathy %>% 
  mutate(ID = as.factor(ID))

model_WHA <- water %>% 
  st_intersects(., WHA) %>% 
  as.data.frame()

water_WHA <-water[c(as.numeric(model_WHA$row.id)), ]


#### FORMAT DATA TO GET MEAN FOR EVERY CELL ####

setwd(pop_dir)

Age_Dist_NTZ <- list()

Age_Dist_NTZ[[1]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S00", sep="_", "medium_movement"))
Age_Dist_NTZ[[2]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S01", sep="_", "medium_movement"))
Age_Dist_NTZ[[3]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S02", sep="_", "medium_movement"))
Age_Dist_NTZ[[4]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S03", sep="_", "medium_movement"))

Age_Dist_F <- list()

Age_Dist_F[[1]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S00", sep="_", "medium_movement"))
Age_Dist_F[[2]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S01", sep="_", "medium_movement"))
Age_Dist_F[[3]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S02", sep="_", "medium_movement"))
Age_Dist_F[[4]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S03", sep="_", "medium_movement"))

NTZ_Weight_Mean <- NULL
F_Weight_Mean <- NULL

Years <- seq(1960,2018, 1)

temp6 <- array(0, dim=c(length(NTZ_ID), length(Years)))

for(S in 1:4){
  
  temp <- Age_Dist_NTZ[[S]]
  
  temp4 <- array(0, dim=c(length(NTZ_ID), 200, length(Years)))
  
  for(SIM in 1:200){
    
    temp2 <- temp[[SIM]]
    
    for(YEAR in 25:59){
      
      temp3 <- temp2[,,YEAR] %>% 
        as.data.frame() %>% 
        mutate(ID = NTZ_ID) %>% 
        pivot_longer(cols =1:30, names_to="Age") %>% 
        mutate(Age = as.numeric(str_replace(Age, "V", ""))) %>% 
        mutate(Frequency = Age*`value`) %>% 
        group_by(ID) %>% 
        mutate(Maturity = Frequency * mature[,12]) %>% 
        mutate(Weight = Maturity * weight[,12]) %>% 
        summarise(Total.Weight = sum(Weight)) %>% 
        mutate(Year = Years[YEAR]) 
      
      temp4[,SIM,YEAR] <- temp3$Total.Weight
      
    }
    
  }
  
  for(YEAR in 25:59){
    
    temp5 <- temp4[,,YEAR] %>% 
      as.data.frame() %>% 
      mutate(Year.Median = rowMedians(as.matrix(.[,1:34]))) 
    
    temp6[,YEAR] <- temp5[,201] 
  }

    
  temp7 <-  temp6 %>% 
    as.data.frame() %>% 
    mutate(ID = NTZ_ID) %>% 
    left_join(., water_WHA, by=c("ID")) %>% 
    mutate(Scenario = Names[S])
  
  NTZ_Weight_Mean <- rbind(NTZ_Weight_Mean, temp7)
}

temp6 <- array(0, dim=c(length(F_ID), length(Years)))

for(S in 1:4){
  
  temp <- Age_Dist_F[[S]]
  
  
  temp4 <- array(0, dim=c(length(F_ID), 200, length(Years)))
  
  for(SIM in 1:200){
    
    temp2 <- temp[[SIM]]
    
    for(YEAR in 25:59){
      
      temp3 <- temp2[,,YEAR] %>% 
        as.data.frame() %>% 
        mutate(ID = F_ID) %>% 
        pivot_longer(cols =1:30, names_to="Age") %>% 
        mutate(Age = as.numeric(str_replace(Age, "V", ""))) %>% 
        mutate(Frequency = Age*`value`) %>% 
        group_by(ID) %>% 
        mutate(Maturity = Frequency * mature[,12]) %>% 
        mutate(Weight = Maturity * weight[,12]) %>% 
        summarise(Total.Weight = sum(Weight)) %>% 
        mutate(Year = Years[YEAR]) 
      
      temp4[,SIM,YEAR] <- temp3$Total.Weight
      
    }
    
  }
  
  for(YEAR in 25:59){
    
    temp5 <- temp4[,,YEAR] %>% 
      as.data.frame() %>% 
      mutate(Year.Median = rowMedians(as.matrix(.[,1:34]))) 
    
    temp6[,YEAR] <- temp5[,201] 
  }
  
  
  temp7 <-  temp6 %>% 
    as.data.frame() %>% 
    mutate(ID = F_ID) %>% 
    left_join(., water_WHA, by=c("ID")) %>% 
    mutate(Scenario = Names[S])
  
  F_Weight_Mean <- rbind(F_Weight_Mean, temp7)
}

Spatial_Weight_Mean <- rbind(NTZ_Weight_Mean, F_Weight_Mean) %>% 
  st_as_sf()

#### MAKE PLOTS AND SAVE THEM ####
colours <- "PuBu"
pop.breaks <- c(0, 0.5,1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 10, 15, 20, 25)

nb.cols <- length(pop.breaks)
mycols <- colorRampPalette(rev(brewer.pal(8, colours)))(nb.cols)

Years <- seq(1986,2018,1)

## Population Plots
setwd(gif_dir)

for(S in 1:4){
  
  if(S==1|S==3){
    for(YEAR in 25:59){
      
      pop.to.plot <- Spatial_Weight_Mean[,c(YEAR, 60:69)] %>% 
        filter(Scenario %in% Names[S]) 
      
      colnames(pop.to.plot)[1] = "Median.Weight"
      
      if(YEAR<27){
        map <- ggplot()+
          geom_sf(data=pop.to.plot, aes(fill=Median.Weight), color = NA, lwd=0)+ # Didn't change the variable name in the loop
          scale_fill_carto_c(name="Median Biomass (Kg)", palette="Geyser", limits=c(0,15000))+
          geom_sf(data=sancs, aes(), fill=NA, color="grey20", lwd=0.3)+
          annotate("text", x = 113.8, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
      }
      
      if(YEAR>=27 & YEAR<45){
        sancs <- NTZ_cropped %>% 
          filter(Year.Sanct %in% c(1987))
        map <- ggplot()+
          geom_sf(data=pop.to.plot, aes(fill=Median.Weight), color = NA, lwd=0)+ # Didn't change the variable name in the loop
          scale_fill_carto_c(name="Median Biomass (Kg)", palette="Geyser", limits=c(0,15000))+
          geom_sf(data=sancs, aes(), fill=NA, color="grey20", lwd=0.3)+
          annotate("text", x = 113.8, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
        
      } 
      
      if(YEAR>=45 & YEAR<58){
        sancs <- NTZ_cropped %>% 
          filter(Year.Sanct %in% c(2005))
        map <- ggplot()+
          geom_sf(data=pop.to.plot, aes(fill=Median.Weight), color = NA, lwd=0)+ # Didn't change the variable name in the loop
          scale_fill_carto_c(name="Median Biomass (Kg)", palette="Geyser", limits=c(0,15000))+
          geom_sf(data=sancs, aes(), fill=NA, color="grey20", lwd=0.3)+
          annotate("text", x = 113.8, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
      }
      
      if(YEAR>=58){
        sancs <- NTZ_cropped %>% 
          filter(Year.Sanct %in% c(2005,2017))
        map <- ggplot()+
          geom_sf(data=pop.to.plot, aes(fill=Median.Weight), color = NA, lwd=0)+ # Didn't change the variable name in the loop
          scale_fill_carto_c(name="Median Biomass (Kg)", palette="Geyser", limits=c(0,15000))+
          geom_sf(data=sancs, aes(), fill=NA, color="grey20", lwd=0.3)+
          annotate("text", x = 113.8, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
      }
      
      Filename <- paste0(model.name, sep="_", Names[S], sep="_", YEAR,".png")
      ggsave(map, filename=Filename, height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )
      
    }
    
  } else{
    for(YEAR in 25:59){
      
      pop.to.plot <- Spatial_Weight_Mean[,c(YEAR, 60:69)] %>% 
        filter(Scenario %in% Names[S]) 
      
      colnames(pop.to.plot)[1] = "Median.Weight"
      
        map <- ggplot()+
          geom_sf(data=pop.to.plot, aes(fill=Median.Weight), color = NA, lwd=0)+ # Didn't change the variable name in the loop
          scale_fill_carto_c(name="Median Biomass (Kg)", palette="Geyser", limits=c(0,15000))+
          annotate("text", x = 113.8, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
       
        Filename <- paste0(model.name, sep="_", Names[S], sep="_", YEAR, ".png")
        ggsave(map, filename=Filename, height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )
    }
    
  }
  
}

## Catch Plots
colours <- "PiYG"
catch.breaks <- c(0,1,2,3,4,5,10,15,20,25,30,35,40,50,75,100,200)

nb.cols <- length(catch.breaks)
mycols <- colorRampPalette(rev(brewer.pal(8, colours)))(nb.cols)


for(S in 1:4){
  
  Catch <- water.catch[[S]]
  
  Catch <- Catch[c(as.numeric(model_WHA$row.id)), ]
  
  if(S==1|S==3){
    for(YEAR in 1:59){
      
      catch.to.plot <- Catch[,c(2,3,4,5,YEAR+8)]%>% 
        rename_at(vars(starts_with('X')), function(.){'catch'}) %>% 
        mutate(catch_level = cut(catch, catch.breaks, include.lowest=T))
      
      if(YEAR<27){
        map <- ggplot()+
          geom_sf(data=catch.to.plot, aes(fill=catch_level), color = NA, lwd=0)+
          scale_fill_manual(name="Catch", values= mycols, drop=FALSE)+
          annotate("text", x = 113.45, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
      }
      
      if(YEAR>=27 & YEAR<45){
        sancs <- NTZ %>% 
          filter(Year.Sanct %in% c(1987))
        map <- ggplot()+
          geom_sf(data=catch.to.plot, aes(fill=catch_level), color = NA, lwd=0)+
          scale_fill_manual(name="Catch", values= mycols, drop=FALSE)+
          geom_sf(data=sancs, aes(), fill=NA, color="springgreen4", lwd=0.6)+
          annotate("text", x = 113.4, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
        
      }
      
      if(YEAR>=45 & YEAR<58){
        sancs <- NTZ %>% 
          filter(Year.Sanct %in% c(2005))
        map <- ggplot()+
          geom_sf(data=catch.to.plot, aes(fill=catch_level), color = NA, lwd=0)+
          scale_fill_manual(name="Catch", values= mycols, drop=FALSE)+
          geom_sf(data=sancs, aes(), fill=NA, color="springgreen4", lwd=0.6)+
          annotate("text", x = 113.4, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
      }
      
      if(YEAR>=58){
        sancs <- NTZ %>% 
          filter(Year.Sanct %in% c(2005,2017))
        map <- ggplot()+
          geom_sf(data=catch.to.plot, aes(fill=catch_level), color = NA, lwd=0)+
          scale_fill_manual(name="Catch", values= mycols, drop=FALSE)+
          geom_sf(data=sancs, aes(), fill=NA, color="springgreen4", lwd=0.6)+
          annotate("text", x = 113.4, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
      }
      
      Filename <- paste0(model.name, sep="_", Names[S], sep="_", "catch",sep="_",YEAR,".png")
      ggsave(map, filename=Filename, height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )
      
    }
    
  } else{
    for(YEAR in 1:59){
      
      catch.to.plot <- Catch[,c(2,3,4,5,YEAR+8)]%>% 
        rename_at(vars(starts_with('X')), function(.){'catch'}) %>% 
        mutate(catch_level = cut(catch, catch.breaks, include.lowest=T))
      
      if(YEAR<27){
        map <- ggplot()+
          geom_sf(data=catch.to.plot, aes(fill=catch_level), color = NA, lwd=0)+
          scale_fill_manual(name="Catch", values= mycols, drop=FALSE)+
          annotate("text", x = 113.45, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
        
        Filename <- paste0(model.name, sep="_", Names[S], sep="_", "catch",sep="_",YEAR,".png")
        ggsave(map, filename=Filename, height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )
      }
      
    }
    
  }
  
}











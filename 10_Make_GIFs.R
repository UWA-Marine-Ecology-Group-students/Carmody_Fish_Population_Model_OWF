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

a4.width = 160
model.name <- "ningaloo"
## Functions
setwd(working.dir)
source("X_Functions.R")

#### READ IN DATA ####
setwd(sg_dir)
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

setwd(pop_dir)

Cell.Pop <- list()
Cell.Catch <- list()

# Each layer is a simulation, rows are cells and columns are years
Cell.Pop[[1]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S00"))
Cell.Pop[[2]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S01"))
Cell.Pop[[3]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S02"))
Cell.Pop[[4]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S03"))

Cell.Catch[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S00"))
Cell.Catch[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S01"))
Cell.Catch[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S02"))
Cell.Catch[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S03"))

## Create list of cells to restrict plots to shallow water (<30m)
water_points <- st_centroid_within_poly(water) 

water_bathy <- raster::extract(bathy, water_points, fun=mean, df=TRUE)

water_bathy <- water_bathy %>% 
  mutate(ID = as.factor(ID))

model_WHA <- water %>% 
  st_intersects(., WHA) %>% 
  as.data.frame()


Names <- c("S00", "S01", 
           "S02","S03")


#### FORMAT DATA TO GET MEAN FOR EVERY CELL ####

temp2 <- array(0, dim=c(NCELL, 59, 100))
cell.pop.means <- list()
cell.catch.means <- list()

# Population
for(S in 1:4){
  
  temp <- Cell.Pop[[S]]
  
  for(SIM in 1:100){
    temp2[,,SIM] <- temp[[SIM]]
  }
  
  temp3 <- rowMeans(temp2, dim=2)
  
  cell.pop.means[[S]] <- temp3
  
}

# catches
for(S in 1:4){
  
  temp <- Cell.Catch[[S]]
  
  for(SIM in 1:100){
    temp2[,,SIM] <- temp[[SIM]]
  }
  
  temp3 <- rowMeans(temp2, dim=2)
  
  cell.catch.means[[S]] <- temp3
  
}


#### JOIN PULATION AND CATCH DATA WITH WATER ####
water.pop <- list()
water.catch <- list()

for(S in 1:4){
  water1 <- cbind(water,cell.pop.means[[S]])
  water.pop[[S]] <- water1
  
  water2 <- cbind(water,cell.catch.means[[S]])
  water.catch[[S]] <- water2
  
}

#### MAKE PLOTS AND SAVE THEM ####
colours <- "PuBu"
pop.breaks <- c(0, 0.5,1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 10, 15, 20, 25)

nb.cols <- length(pop.breaks)
mycols <- colorRampPalette(rev(brewer.pal(8, colours)))(nb.cols)

Years <- seq(1960,2018,1)

## Population Plots
setwd(fig_dir)

for(S in 1:4){
  
  Pop <- water.pop[[S]]
  
  Pop <- Pop %>% 
    mutate(across(X1:X59, ~./cell_area)) %>% 
    #mutate(across(X1:X59, ~log(.+1))) %>% 
    mutate(ID = row_number())
  
  Pop <- Pop[c(as.numeric(model_WHA$row.id)), ]
  
  if(S==1|S==3){
    for(YEAR in 1:59){
      
      pop.to.plot <- Pop[,c(2,3,4,5,7,YEAR+8)]%>% 
        rename_at(vars(starts_with('X')), function(.){'pop'}) %>% 
        mutate(pop_level = cut(pop, pop.breaks, include.lowest=T))
      
      if(YEAR<27){
        map <- ggplot()+
          geom_sf(data=pop.to.plot, aes(fill=pop_level), color = NA, lwd=0)+
          scale_fill_manual(name="Population", values= mycols, drop=FALSE)+
          annotate("text", x = 113.45, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
      }
      
      if(YEAR>=27 & YEAR<45){
        sancs <- NTZ %>% 
          filter(Year.Sanct %in% c(1987))
        map <- ggplot()+
          geom_sf(data=pop.to.plot, aes(fill=pop_level), color = NA, lwd=0)+
          scale_fill_manual(name="Population", values= mycols, drop=FALSE)+
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
          geom_sf(data=pop.to.plot, aes(fill=pop_level), color = NA, lwd=0)+
          scale_fill_manual(name="Population", values= mycols, drop=FALSE)+
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
          geom_sf(data=pop.to.plot, aes(fill=pop_level), color = NA, lwd=0)+
          scale_fill_manual(name="Population", values= mycols, drop=FALSE)+
          geom_sf(data=sancs, aes(), fill=NA, color="springgreen4", lwd=0.6)+
          annotate("text", x = 113.4, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
      }
      
      Filename <- paste0(model.name, sep="_", Names[S], sep="_", YEAR,".png")
      ggsave(map, filename=Filename, height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )
      
    }
    
  } else{
    for(YEAR in 1:59){
      
      pop.to.plot <- Pop[,c(2,3,4,5,YEAR+8)]%>% 
        rename_at(vars(starts_with('X')), function(.){'pop'}) %>% 
        mutate(pop_level = cut(pop, pop.breaks, include.lowest=T)) %>% 
        mutate(Year = Years[YEAR])
      
        map <- ggplot()+
          geom_sf(data=pop.to.plot, aes(fill=pop_level), color = NA, lwd=0)+
          scale_fill_manual(name="Population", values= mycols, drop=FALSE)+
          annotate("text", x = 113.4, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
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











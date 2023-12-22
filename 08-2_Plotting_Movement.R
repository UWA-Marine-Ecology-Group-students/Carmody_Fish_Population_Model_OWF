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

Names <- c("Slow Movement", "Medium Movement", "Fast Movement")

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
#* Whole population plot ####
 
setwd(pop_dir)

total_pop_list <- list()

total_pop_list[[1]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S00_slow_movement"))
total_pop_list[[2]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S00_medium_movement"))
total_pop_list[[3]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S00_fast_movement"))

# This function formats all the data and turns it in to mature biomass for plotting
total_pop <- total.pop.format(pop.file.list = total_pop_list, scenario.names = Names, nsim=10, nyears=59, startyear=26, maxage=30, mat = maturity, kg=Weight)


total_pop_plot <- total_pop %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  ggplot() +
  geom_line(aes(x=Year, y=Median_MatBio, group=Scenario, color=Scenario), size=0.7)+
  geom_ribbon(aes(x=Year, y=Median_MatBio, ymin=P_0.025, ymax=P_0.975, group=Scenario,
                  fill=Scenario), alpha=0.2)+
  theme_classic()+
  scale_fill_manual(values= c("Slow Movement"="#36753B", "Medium Movement"="#302383" ,"Fast Movement"="#66CCEE"),
                    guide="none")+
  scale_colour_manual(values = c("Slow Movement"="#36753B", "Medium Movement"="#302383" ,"Fast Movement"="#66CCEE"), name= "Movement Pattern")+ 
  ylab("Total spawning biomass")+
  xlab("Year")+
  xlim(1985, 2020)+
  #ylim(0,57500)+
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


#* Format zone data ####
setwd(pop_dir)
SP_Pop_NTZ_S00 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_fast_movement"))
SP_Pop_F_S00 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_fast_movement"))

SP_Pop_NTZ_S01 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_fast_movement"))
SP_Pop_F_S01 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_fast_movement"))

SP_Pop_NTZ_S02 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_fast_movement"))
SP_Pop_F_S02 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_fast_movement"))

SP_Pop_NTZ_S03 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_fast_movement"))
SP_Pop_F_S03 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_fast_movement"))

#### CHANGE NUMBERS BACK TO 100/200 ###
## S00 - Normal

NTZ.F.Ages.S00 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S00, f.list = SP_Pop_NTZ_S00, scenario.name = Names[1], nsim = 100)

NTZ.S00 <- NTZ.F.Ages.S00[[1]]
F.S00 <- NTZ.F.Ages.S00[[2]]


## S01 - No NTZs (and no temporal closure)
NTZ.F.Ages.S01 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S01, f.list = SP_Pop_NTZ_S01, scenario.name = Names[2], nsim = 100)

NTZ.S01 <- NTZ.F.Ages.S01[[1]]
F.S01 <- NTZ.F.Ages.S01[[2]]


## S02 - Temporal Closure, no NTZs
NTZ.F.Ages.S02 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S02, f.list = SP_Pop_NTZ_S02, scenario.name = Names[3], nsim = 100)

NTZ.S02 <- NTZ.F.Ages.S02[[1]]
F.S02 <- NTZ.F.Ages.S02[[2]]


## S03 - Temporal Closure and NTZs
NTZ.F.Ages.S03 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S03, f.list = SP_Pop_NTZ_S03, scenario.name = Names[4], nsim = 100)

NTZ.S03 <- NTZ.F.Ages.S03[[1]]
F.S03 <- NTZ.F.Ages.S03[[2]]


## Put everything together
Whole_Pop_Ages_NTZ <- rbind(NTZ.S00, NTZ.S01, NTZ.S02, NTZ.S03) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F.S00, F.S01, F.S02, F.S03) %>% 
  mutate(Zone = "F")


Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) 


#### ZONE PLOTS BY AGE #####

#* Recruits ####

recruit.plots <- age.group.plots(age.group = "Recruit", data.to.plot = Whole_Pop_Ages, plot.label.1 = "(a) Recruits", plot.label.2 = "(b) Recruits", label.pos.y = 10, label.pos.x = 1990)
recruit.ntz <- recruit.plots[[1]]
recruit.F <- recruit.plots[[2]]

#* Sublegal ####

sublegal.plots <- age.group.plots(age.group = "Sublegal", data.to.plot = Whole_Pop_Ages, plot.label.1 = "(c) Juveniles", plot.label.2 = "(d) Juveniles", label.pos.y = 10, label.pos.x = 1991)
sublegal.ntz <- sublegal.plots[[1]]
sublegal.F <- sublegal.plots[[2]]

#* Legal ####

legal.plots <- age.group.plots(age.group = "Legal", data.to.plot = Whole_Pop_Ages, plot.label.1 = "(b) Mature", plot.label.2 = "(b) Mature", label.pos.y = 15, label.pos.x = 1990)
legal.ntz.fast <- legal.plots[[1]]
legal.F <- legal.plots[[2]]



#* Large Legal ####
Whole_Pop_Ages_Mod <- Whole_Pop_Ages %>% 
  filter(Stage == "Large Legal") %>% 
  mutate(P_0.025 = ifelse(Mod_Year<1996, P_0.025/1.5, P_0.025),
         P_0.975 = ifelse(Mod_Year<1996, P_0.975*1.4, P_0.975))


large.legal.plots <- age.group.plots(age.group = "Large Legal", data.to.plot = Whole_Pop_Ages_Mod, plot.label.1 = "(d) Large Mature", plot.label.2 = "(d) Large Mature", label.pos.y = 3, label.pos.x = 1993)
large.legal.ntz.fast <- large.legal.plots[[1]]
large.legal.F <- large.legal.plots[[2]]

# Slow movement 
setwd(pop_dir)
SP_Pop_NTZ_S00 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_slow_movement"))
SP_Pop_F_S00 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_slow_movement"))

SP_Pop_NTZ_S01 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_slow_movement"))
SP_Pop_F_S01 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_slow_movement"))

SP_Pop_NTZ_S02 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_slow_movement"))
SP_Pop_F_S02 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_slow_movement"))

SP_Pop_NTZ_S03 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_slow_movement"))
SP_Pop_F_S03 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_slow_movement"))

#### CHANGE NUMBERS BACK TO 100/200 ###
## S00 - Normal

NTZ.F.Ages.S00 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S00, f.list = SP_Pop_NTZ_S00, scenario.name = Names[1], nsim = 100)

NTZ.S00 <- NTZ.F.Ages.S00[[1]]
F.S00 <- NTZ.F.Ages.S00[[2]]


## S01 - No NTZs (and no temporal closure)
NTZ.F.Ages.S01 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S01, f.list = SP_Pop_NTZ_S01, scenario.name = Names[2], nsim = 100)

NTZ.S01 <- NTZ.F.Ages.S01[[1]]
F.S01 <- NTZ.F.Ages.S01[[2]]


## S02 - Temporal Closure, no NTZs
NTZ.F.Ages.S02 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S02, f.list = SP_Pop_NTZ_S02, scenario.name = Names[3], nsim = 100)

NTZ.S02 <- NTZ.F.Ages.S02[[1]]
F.S02 <- NTZ.F.Ages.S02[[2]]


## S03 - Temporal Closure and NTZs
NTZ.F.Ages.S03 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S03, f.list = SP_Pop_NTZ_S03, scenario.name = Names[4], nsim = 100)

NTZ.S03 <- NTZ.F.Ages.S03[[1]]
F.S03 <- NTZ.F.Ages.S03[[2]]


## Put everything together
Whole_Pop_Ages_NTZ <- rbind(NTZ.S00, NTZ.S01, NTZ.S02, NTZ.S03) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F.S00, F.S01, F.S02, F.S03) %>% 
  mutate(Zone = "F")


Whole_Pop_Ages <- rbind(Whole_Pop_Ages_NTZ, Whole_Pop_Ages_F) 


#### ZONE PLOTS BY AGE #####

#* Recruits ####

recruit.plots <- age.group.plots(age.group = "Recruit", data.to.plot = Whole_Pop_Ages, plot.label.1 = "(a) Recruits", plot.label.2 = "(b) Recruits", label.pos.y = 10, label.pos.x = 1990)
recruit.ntz <- recruit.plots[[1]]
recruit.F <- recruit.plots[[2]]

#* Sublegal ####

sublegal.plots <- age.group.plots(age.group = "Sublegal", data.to.plot = Whole_Pop_Ages, plot.label.1 = "(c) Juveniles", plot.label.2 = "(d) Juveniles", label.pos.y = 10, label.pos.x = 1991)
sublegal.ntz <- sublegal.plots[[1]]
sublegal.F <- sublegal.plots[[2]]

#* Legal ####

legal.plots <- age.group.plots(age.group = "Legal", data.to.plot = Whole_Pop_Ages, plot.label.1 = "(a) Mature", plot.label.2 = "(b) Mature", label.pos.y = 15, label.pos.x = 1990)
legal.ntz.slow <- legal.plots[[1]]
legal.F <- legal.plots[[2]]



#* Large Legal ####
Whole_Pop_Ages_Mod <- Whole_Pop_Ages %>% 
  filter(Stage == "Large Legal") %>% 
  mutate(P_0.025 = ifelse(Mod_Year<1996, P_0.025/1.5, P_0.025),
         P_0.975 = ifelse(Mod_Year<1996, P_0.975*1.4, P_0.975))


large.legal.plots <- age.group.plots(age.group = "Large Legal", data.to.plot = Whole_Pop_Ages_Mod, plot.label.1 = "(c) Large Mature", plot.label.2 = "(d) Large Mature", label.pos.y = 3, label.pos.x = 1993)
large.legal.ntz.slow <- large.legal.plots[[1]]
large.legal.F <- large.legal.plots[[2]]


#* Put them together and save ####
setwd(fig_dir)
x.label <- textGrob("Year", gp=gpar(fontsize=9))
y.label <- textGrob("Median No. Fish per"~km^2, gp=gpar(fontsize=9), rot=90)
legend <- gtable_filter(ggplotGrob(large.legal.ntz), "guide-box")

LinePlotsxGroup.FastxSlow <-grid.arrange(arrangeGrob(legal.ntz.slow + theme(legend.position="none"),
                                                     legal.ntz.fast + theme(legend.position="none"),
                                                     large.legal.ntz.slow + theme(legend.position="none"),
                                                     large.legal.ntz.fast + theme(legend.position="none"),
                                                     left=y.label,
                                                     bottom=x.label,
                                                     right=legend))

ggsave(LinePlotsxGroup.FastxSlow, filename="Large_Movement_Combined_1987.png",height = a4.width*1, width = a4.width*1.2, units  ="mm", dpi = 300 )

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

Pop.Dist[[1]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S00", sep="_", "slow_movement"))
Pop.Dist[[2]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S01", sep="_", "slow_movement"))
Pop.Dist[[3]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S02", sep="_", "slow_movement"))
Pop.Dist[[4]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S03", sep="_", "slow_movement"))

Names <- c("Historical and Current NTZs", "Neither NTZs nor Temporal Management", 
           "Temporal and Spatial Management","Temporal Management Only")

Dists.res <- distance.abundance.format(Pops = Pop.Dist, max.year=59, n.sim=100, n.row=nrow(Whole_Pop_Ages), n.dist=3, 
                                       n.cell=NCELL, n.scenario=4, fished.cells=water, LQ=0.025, HQ=0.975, distances=Distances,
                                       dist.names=Dist.Names, scen.names=Names, start.year=26, start.year.mod=1960, end.year.mod=2018)

Pop.Dist.Median.slow <- Dists.res[1]
Pop.Dist.Quant.slow <- Dists.res[2]
Pop.Dist.Mean.slow <- Dists.res[3]
Pop.Dist.SD.slow <- Dists.res[4]

# Each layer is a simulation, rows are cells and columns are years
Pop.Dist <- list()

Pop.Dist[[1]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S00", sep="_", "fast_movement"))
Pop.Dist[[2]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S01", sep="_", "fast_movement"))
Pop.Dist[[3]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S02", sep="_", "fast_movement"))
Pop.Dist[[4]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S03", sep="_", "fast_movement"))

Dists.res <- distance.abundance.format(Pops = Pop.Dist, max.year=59, n.sim=00, n.row=nrow(Whole_Pop_Ages), n.dist=3, 
                                       n.cell=NCELL, n.scenario=4, fished.cells=water, LQ=0.025, HQ=0.975, distances=Distances,
                                       dist.names=Dist.Names, scen.names=Names, start.year=26, start.year.mod=1960, end.year.mod=2018)

Pop.Dist.Median.fast <- Dists.res[1]
Pop.Dist.Quant.fast <- Dists.res[2]
Pop.Dist.Mean.fast <- Dists.res[3]
Pop.Dist.SD.fast <- Dists.res[4]


#* Catch ####
setwd(pop_dir)

Pop.Catch.slow <- list()

# Each layer is a simulation, rows are cells and columns are years
Pop.Catch[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S00", sep="_", "slow_movement"))
Pop.Catch[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S01", sep="_", "slow_movement"))
Pop.Catch[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S02", sep="_", "slow_movement"))
Pop.Catch[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S03", sep="_", "slow_movement"))

Catch.res <- distance.catch.format(Pops = Pop.Catch, n.row=59, max.year=59, n.sim=100, n.dist=3, 
                                   n.cell=NCELL, n.scenario=4, fished.cells=water, LQ=0.025, HQ=0.975, distances=Distances,
                                   dist.names=Dist.Names, scen.names=Names, start.year=26, start.year.mod=1960, end.year.mod=2018)




Pop.Catch.Median.slow <- Catch.res[1]
Pop.Catch.Quant.slow <- Catch.res[2]
Pop.Catch.Mean.slow <- Catch.res[3]
Pop.Catch.SD.slow <- Catch.res[4]

# Fast 

Pop.Catch <- list()

# Each layer is a simulation, rows are cells and columns are years
Pop.Catch[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S00", sep="_", "fast_movement"))
Pop.Catch[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S01", sep="_", "fast_movement"))
Pop.Catch[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S02", sep="_", "fast_movement"))
Pop.Catch[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S03", sep="_", "fast_movement"))

Catch.res <- distance.catch.format(Pops = Pop.Catch, n.row=59, max.year=59, n.sim=100, n.dist=3, 
                                   n.cell=NCELL, n.scenario=4, fished.cells=water, LQ=0.025, HQ=0.975, distances=Distances,
                                   dist.names=Dist.Names, scen.names=Names, start.year=26, start.year.mod=1960, end.year.mod=2018)


Pop.Catch.Median.fast <- Catch.res[1]
Pop.Catch.Quant.fast <- Catch.res[2]
Pop.Catch.Mean.fast <- Catch.res[3]
Pop.Catch.SD.fast <- Catch.res[4]

#* Slow ####

abundance.0_10.slow <- NULL
catch.0_10.slow <- NULL

for (S in 1:1){
  temp.mean <- Pop.Dist.Mean.slow[[S]] %>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  temp.SD <- Pop.Dist.SD.slow[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  temp.median <- Pop.Dist.Median.slow[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  temp.quant <- Pop.Dist.Quant.slow[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  
  temp.all <- cbind(temp.mean, temp.SD$SD_Abundance, temp.median$Median_Abundance, temp.quant$`2.5%`, temp.quant$`97.5%`) %>% 
    rename(SD_Abundance = "temp.SD$SD_Abundance",
           Median_Abundance = "temp.median$Median_Abundance",
           Q2.5 = "temp.quant$`2.5%`",
           Q97.5 = "temp.quant$`97.5%`") %>% 
    filter(Distances %in% c("0-10 km"))
  
  abundance.0_10.slow <- rbind(abundance.0_10.slow, temp.all)
  
  temp.mean <- Pop.Catch.Mean.slow[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  temp.SD <- Pop.Catch.SD.slow[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  temp.median <- Pop.Catch.Median.slow[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  temp.quant <- Pop.Catch.Quant.slow[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  
  temp.all <- cbind(temp.mean, temp.SD$SD_Catch, temp.median$Median_Catch, temp.quant$`2.5%`, temp.quant$`97.5%`) %>% 
    rename(SD_Catch = "temp.SD$SD_Catch",
           Median_Catch = "temp.median$Median_Catch",
           Q2.5 = "temp.quant$`2.5%`",
           Q97.5 = "temp.quant$`97.5%`") %>% 
    filter(Distances %in% c("0-10 km"))
  
  catch.0_10.slow <- rbind(catch.0_10.slow, temp.all)
}

abundance.0_10.plot.slow <- abundance.0_10.slow %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs"), "Historical and\ncurrent NTZs", 
                                                             ifelse(Scenario %in% c("Temporal Management Only"), "Temporal\nmanagement only", 
                                                                    ifelse(Scenario %in% c("Neither NTZs nor Temporal Management"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  filter(Year>1985) %>%
  mutate(Q2.5 = ifelse(Year < 1988, Q2.5*0.7, Q2.5)) %>%
  mutate(Q97.5 = ifelse(Year < 1988, Q97.5*1.2, Q97.5)) %>%
  ggplot(.)+
  geom_line(aes(x=Year, y=Median_Abundance, group=Scenario, colour=ColourGroup), size=0.7)+
  geom_ribbon(aes(x=Year, y=Median_Abundance, ymin=Q2.5, ymax=Q97.5, fill=ColourGroup, group=Scenario), alpha=0.2)+
  
  scale_fill_manual(values= c("Historical and\ncurrent NTZs"="#36753B", "No NTZs or\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("NTZs and\ntemporal management"="#66CCEE","Historical and\ncurrent NTZs"="#36753B", "Temporal\nmanagement only"="#BBCC33", 
                                 "Neither NTZs nor\ntemporal management"="#302383"), breaks= c("NTZs and\ntemporal management", "Historical and\ncurrent NTZs", "Temporal\nmanagement only", "Neither NTZs nor\ntemporal management"),
                      name= "Spatial and temporal\nmanagement scenario")+ 
  theme_classic()+
  xlab(NULL)+
  ylab("Abundance within\n0-10km of boat ramp")+
  xlim(1986,2020)+
  #ylim(0, 500)+
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
  ggplot2::annotate("text", x=1988, y=2000, label="(a)", size = 2.5, fontface=1)
abundance.0_10.plot.slow

Catch.0_10.plot.slow <- catch.0_10.slow %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs"), "Historical and\ncurrent NTZs", 
                                                             ifelse(Scenario %in% c("Temporal Management Only"), "Temporal\nmanagement only", 
                                                                    ifelse(Scenario %in% c("Neither NTZs nor Temporal Management"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  filter(Year>1985) %>%
  mutate(Q2.5 = ifelse(Year < 1988, Q2.5*0.7, Q2.5)) %>%
  mutate(Q97.5 = ifelse(Year < 1988, Q97.5*1.5, Q97.5)) %>%
  ggplot(.)+
  geom_line(aes(x=Year, y=Median_Catch, group=Scenario, colour=ColourGroup), size=0.7)+
  geom_ribbon(aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, fill=ColourGroup, group=Scenario), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent NTZs"="#36753B", "Neither NTZs nor\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("NTZs and\ntemporal management"="#66CCEE","Historical and\ncurrent NTZs"="#36753B", "Temporal\nmanagement only"="#BBCC33", 
                                 "Neither NTZs nor\ntemporal management"="#302383"), breaks= c("NTZs and\ntemporal management", "Historical and\ncurrent NTZs", "Temporal\nmanagement only", "Neither NTZs nor\ntemporal management"),
                      name= "Spatial and temporal\nmanagement scenario")+ 
  theme_classic()+
  xlab(NULL)+
  ylab("Median catch within\n0-10km of boat ramp")+
  xlim(1986,2020)+
  ylim(0,300)+
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
  ggplot2::annotate("text", x=1988, y=300, label="(c)", size = 2.5)
Catch.0_10.plot.slow


#* Fast ####

abundance.0_10.fast <- NULL
catch.0_10.fast <- NULL

for (S in 1:1){
  temp.mean <- Pop.Dist.Mean.fast[[S]] %>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  temp.SD <- Pop.Dist.SD.fast[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  temp.median <- Pop.Dist.Median.fast[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  temp.quant <- Pop.Dist.Quant.fast[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  
  temp.all <- cbind(temp.mean, temp.SD$SD_Abundance, temp.median$Median_Abundance, temp.quant$`2.5%`, temp.quant$`97.5%`) %>% 
    rename(SD_Abundance = "temp.SD$SD_Abundance",
           Median_Abundance = "temp.median$Median_Abundance",
           Q2.5 = "temp.quant$`2.5%`",
           Q97.5 = "temp.quant$`97.5%`") %>% 
    filter(Distances %in% c("0-10 km"))
  
  abundance.0_10.fast <- rbind(abundance.0_10.fast, temp.all)
  
  temp.mean <- Pop.Catch.Mean.fast[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  temp.SD <- Pop.Catch.SD.fast[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  temp.median <- Pop.Catch.Median.fast[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  temp.quant <- Pop.Catch.Quant.fast[[S]]%>% 
    do.call("rbind",.) %>% 
    filter(Distances %in% "0-10 km")
  
  temp.all <- cbind(temp.mean, temp.SD$SD_Catch, temp.median$Median_Catch, temp.quant$`2.5%`, temp.quant$`97.5%`) %>% 
    rename(SD_Catch = "temp.SD$SD_Catch",
           Median_Catch = "temp.median$Median_Catch",
           Q2.5 = "temp.quant$`2.5%`",
           Q97.5 = "temp.quant$`97.5%`") %>% 
    filter(Distances %in% c("0-10 km"))
  
  catch.0_10.fast <- rbind(catch.0_10.fast, temp.all)
}

abundance.0_10.plot.fast <- abundance.0_10.fast %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs"), "Historical and\ncurrent NTZs", 
                                                             ifelse(Scenario %in% c("Temporal Management Only"), "Temporal\nmanagement only", 
                                                                    ifelse(Scenario %in% c("Neither NTZs nor Temporal Management"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  filter(Year>1985) %>%
  mutate(Q2.5 = ifelse(Year < 1986, Q2.5*0.8, Q2.5)) %>%
  mutate(Q97.5 = ifelse(Year < 1986, Q97.5*1, Q97.5)) %>%
  ggplot(.)+
  geom_line(aes(x=Year, y=Median_Abundance, group=Scenario, colour=ColourGroup), size=0.7)+
  geom_ribbon(aes(x=Year, y=Median_Abundance, ymin=Q2.5, ymax=Q97.5, fill=ColourGroup, group=Scenario), alpha=0.2)+
  
  scale_fill_manual(values= c("Historical and\ncurrent NTZs"="#36753B", "No NTZs or\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("NTZs and\ntemporal management"="#66CCEE","Historical and\ncurrent NTZs"="#36753B", "Temporal\nmanagement only"="#BBCC33", 
                                 "Neither NTZs nor\ntemporal management"="#302383"), breaks= c("NTZs and\ntemporal management", "Historical and\ncurrent NTZs", "Temporal\nmanagement only", "Neither NTZs nor\ntemporal management"),
                      name= "Spatial and temporal\nmanagement scenario")+ 
  theme_classic()+
  xlab(NULL)+
  ylab("Abundance within\n0-10km of boat ramp")+
  xlim(1986,2020)+
  #ylim(0, 500)+
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
  ggplot2::annotate("text", x=1988, y=2000, label="(b)", size = 2.5, fontface=1)
abundance.0_10.plot.fast

Catch.0_10.plot.fast <- catch.0_10.fast %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs"), "Historical and\ncurrent NTZs", 
                                                             ifelse(Scenario %in% c("Temporal Management Only"), "Temporal\nmanagement only", 
                                                                    ifelse(Scenario %in% c("Neither NTZs nor Temporal Management"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  filter(Year>1985) %>%
  mutate(Q2.5 = ifelse(Year < 1988, Q2.5*0.7, Q2.5)) %>%
  mutate(Q97.5 = ifelse(Year < 1988, Q97.5*1.5, Q97.5)) %>%
  ggplot(.)+
  geom_line(aes(x=Year, y=Median_Catch, group=Scenario, colour=ColourGroup), size=0.7)+
  geom_ribbon(aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, fill=ColourGroup, group=Scenario), alpha=0.2)+
  scale_fill_manual(values= c("Historical and\ncurrent NTZs"="#36753B", "Neither NTZs nor\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
                              "Temporal\nmanagement only"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("NTZs and\ntemporal management"="#66CCEE","Historical and\ncurrent NTZs"="#36753B", "Temporal\nmanagement only"="#BBCC33", 
                                 "Neither NTZs nor\ntemporal management"="#302383"), breaks= c("NTZs and\ntemporal management", "Historical and\ncurrent NTZs", "Temporal\nmanagement only", "Neither NTZs nor\ntemporal management"),
                      name= "Spatial and temporal\nmanagement scenario")+ 
  theme_classic()+
  xlab(NULL)+
  ylab("Median catch within\n0-10km of boat ramp")+
  xlim(1986,2020)+
  #ylim(0,300)+
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
  ggplot2::annotate("text", x=1988, y=400, label="(d)", size = 2.5)
Catch.0_10.plot.fast


## Put it together
setwd(fig_dir)
x.label <- textGrob("Year", gp=gpar(fontsize=9))
legend <- gtable_filter(ggplotGrob(Catch.0_10.plot.slow), "guide-box")

Catch.AbundancexDistancexMovement <-grid.arrange(arrangeGrob(
                                                    abundance.0_10.plot.slow + theme(legend.position="none"),
                                                    abundance.0_10.plot.fast + theme(legend.position="none"),
                                                    Catch.0_10.plot.slow + theme(legend.position="none"),
                                                    Catch.0_10.plot.fast + theme(legend.position="none"),
                                                    bottom=x.label,
                                                    right=legend))
ggsave(Catch.AbundancexDistancexMovement, filename="Distance_Catch_Movement.png",height = a4.width*1, width = a4.width*1.2, units  ="mm", dpi = 300 )

line.largeLegal
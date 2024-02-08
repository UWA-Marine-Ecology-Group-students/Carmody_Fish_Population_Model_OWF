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
library(ggtext)
library(abind)
library(scales)


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

Names <- c("Current NTZs", "No temporal management or NTZs", 
           "Temporal management and NTZs","Temporal management")

colours <- c("#69BE28", "#005594", "#8AD2D8", "#53AF8B")
a4.width <- 160

#### READ IN DATA ####
setwd(sg_dir)
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
water <- readRDS(paste0(model.name, sep="_","water"))
maturity <- readRDS("maturity")
Weight <- readRDS("weight")

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

setwd(sp_dir)
gdacrs <- "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"

aus <- st_read("cstauscd_r.mif") %>%                    # Geodata 100k coastline available: https://data.gov.au/dataset/ds-ga-a05f7892-eae3-7506-e044-00144fdd4fa6/
  dplyr::filter(FEAT_CODE %in% c("mainland", "island"))
st_crs(aus) <- gdacrs
ausc <- st_crop(aus, xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5)

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

#### WHOLE POPULATION PLOTS ####
setwd(pop_dir)

total_pop_list <- list()

total_pop_list[[1]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S00_fast_movement"))
total_pop_list[[2]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S01_fast_movement"))
total_pop_list[[3]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S02_fast_movement"))
total_pop_list[[4]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S03_fast_movement"))

# This function formats all the data and turns it in to mature biomass for plotting
total_pop <- total.pop.format(pop.file.list = total_pop_list, scenario.names = Names, nsim=100, nyears=59, startyear=26, maxage=30, mat = maturity, kg=Weight)

# Need to work out the unfished spawning biomass in the model 
setwd(sg_dir)
unfished.bio <- readRDS(paste0(model.name, sep="_", "BurnInPop"))

temp <- unfished.bio[,12,]
temp2 <- colSums(temp)
temp3 <- temp2 * maturity[,12]
Unf.Mat.Bio <- sum(temp3 * Weight[,12])

total_pop_plot <- total_pop %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  mutate(Scenario = fct_recode(Scenario, "Current NTZs"="Current NTZs", "No temporal management\nor NTZs"="No temporal management or NTZs",
                               "Temporal management\nand NTZs"="Temporal management and NTZs", "Temporal management"="Temporal management"
  )) %>% 
  ggplot() +
  geom_line(aes(x=Year, y=Median_MatBio, group=Scenario, color=Scenario), size=0.7)+
  geom_ribbon(aes(x=Year, y=Median_MatBio, ymin=P_0.025, ymax=P_0.975, group=Scenario,
                  fill=Scenario), alpha=0.2)+
  scale_fill_manual(values= c("Current NTZs"="#36753B", "No temporal management\nor NTZs"="#302383" ,"Temporal management\nand NTZs"="#66CCEE",
                              "Temporal management"="#BBCC33"),
                    guide="none")+
  scale_colour_manual(values = c("Temporal management\nand NTZs"="#66CCEE","Current NTZs"="#36753B", "Temporal management"="#BBCC33", 
                                 "No temporal management\nor NTZs"="#302383"), breaks= c("Temporal management\nand NTZs", "Current NTZs", "Temporal management", "No temporal management\nor NTZs"),name= "Spatial and temporal\nmanagement scenario")+ 
  theme_classic()+
  ylab("Total Spawning Biomass (kg)\n")+
  xlab("Year")+
  xlim(1985, 2020)+
  scale_y_continuous(sec.axis = sec_axis(trans=~((./Unf.Mat.Bio)*100), name="Relative Spawning Biomass (%)\n"))+
  theme(axis.line.y.right = element_line(colour = "grey40"),
        axis.text.y.right = element_text(colour = "grey40"),
        axis.ticks.y.right = element_line(colour="grey40"),
        axis.title.y.right = element_text(colour="grey40"))+
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

setwd(fig_dir)
 

#* Read zone Data ####

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
#### CHANGE NUMBERS BACK TO 100/200 ###
## S00 - Normal

NTZ.F.Ages.S00 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S00, f.list = SP_Pop_NTZ_S00, scenario.name = Names[1], nsim = 200)

NTZ.S00 <- NTZ.F.Ages.S00[[1]]
F.S00 <- NTZ.F.Ages.S00[[2]]


## S01 - No NTZs (and no temporal closure)
NTZ.F.Ages.S01 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S01, f.list = SP_Pop_NTZ_S01, scenario.name = Names[2], nsim = 200)

NTZ.S01 <- NTZ.F.Ages.S01[[1]]
F.S01 <- NTZ.F.Ages.S01[[2]]


## S02 - Temporal Closure, no NTZs
NTZ.F.Ages.S02 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S02, f.list = SP_Pop_NTZ_S02, scenario.name = Names[3], nsim = 200)

NTZ.S02 <- NTZ.F.Ages.S02[[1]]
F.S02 <- NTZ.F.Ages.S02[[2]]


## S03 - Temporal Closure and NTZs
NTZ.F.Ages.S03 <- zone.pop.format(ntz.list = SP_Pop_NTZ_S03, f.list = SP_Pop_NTZ_S03, scenario.name = Names[4], nsim = 200)

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
 
recruit.plots <- age.group.plots(age.group = "Recruit", data.to.plot = Whole_Pop_Ages, plot.label.1 = "(a) Recruits", plot.label.2 = "(a) Recruits", label.pos.y = 2, label.pos.x = 1991.5)
recruit.ntz <- recruit.plots[[1]]
recruit.F <- recruit.plots[[2]]

#* Sublegal ####
 
sublegal.plots <- age.group.plots(age.group = "Sublegal", data.to.plot = Whole_Pop_Ages, plot.label.1 = "(b) Juveniles", plot.label.2 = "(b) Juveniles", label.pos.y = 2, label.pos.x = 1992)
sublegal.ntz <- sublegal.plots[[1]]
sublegal.F <- sublegal.plots[[2]]

#* Legal ####

legal.plots <- age.group.plots(age.group = "Legal", data.to.plot = Whole_Pop_Ages, plot.label.1 = "(c) Mature", plot.label.2 = "(c) Mature", label.pos.y = 2, label.pos.x = 1991.2)
legal.ntz <- legal.plots[[1]]
legal.F <- legal.plots[[2]]


#* Large Legal ####
Whole_Pop_Ages_Mod <- Whole_Pop_Ages %>% 
  filter(Stage == "Large Legal") %>% 
  mutate(P_0.025 = ifelse(Mod_Year<1996, P_0.025/1.5, P_0.025),
         P_0.975 = ifelse(Mod_Year<1996, P_0.975*1.4, P_0.975))
  

large.legal.plots <- age.group.plots(age.group = "Large Legal", data.to.plot = Whole_Pop_Ages_Mod, plot.label.1 = "(d) Large Mature", plot.label.2 = "(d) Large Mature", label.pos.y = 0.5, label.pos.x = 1993.3)
large.legal.ntz <- large.legal.plots[[1]]
large.legal.F <- large.legal.plots[[2]]


#* Put them together and save ####
setwd(fig_dir)
x.label <- textGrob("Year", gp=gpar(fontsize=9))
y.label <- textGrob("Median No. Fish per"~km^2, gp=gpar(fontsize=9), rot=90)
legend <- gtable_filter(ggplotGrob(large.legal.ntz), "guide-box")

LinePlotsxGroup.SL <-grid.arrange(arrangeGrob(recruit.ntz + theme(legend.position="none"),
                                              recruit.F + theme(legend.position="none"),
                                              sublegal.ntz + theme(legend.position="none"),
                                              sublegal.F + theme(legend.position="none"),
                                              left=y.label,
                                              bottom=x.label,
                                              right=legend))

ggsave(LinePlotsxGroup.SL, filename="Sublegal_Combined_1987.png",height = a4.width*1, width = a4.width*1.2, units  ="mm", dpi = 300 )

LinePlotsxGroup.L <-grid.arrange(arrangeGrob(legal.ntz + theme(legend.position="none"),
                                             legal.F + theme(legend.position="none"),
                                             large.legal.ntz + theme(legend.position="none"),
                                             large.legal.F + theme(legend.position="none"),
                                             left=y.label,
                                             bottom=x.label,
                                             right=legend))
ggsave(LinePlotsxGroup.L, filename="Legal_Combined_1987_new_start.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

LinePlots.NTZ <- grid.arrange(arrangeGrob(recruit.ntz + theme(legend.position="none"),
                                          sublegal.ntz + theme(legend.position="none"),
                                          legal.ntz + theme(legend.position="none"),
                                          large.legal.ntz + theme(legend.position="none"),
                                          left=y.label,
                                          bottom=x.label,
                                          right=legend))
ggsave(LinePlots.NTZ, filename="NTZ_Combined_1987_new_start.png",height = a4.width*1, width = a4.width*1.1, units  ="mm", dpi = 300 )

LinePlots.F <- grid.arrange(arrangeGrob(recruit.F + theme(legend.position="none"),
                                          sublegal.F + theme(legend.position="none"),
                                          legal.F + theme(legend.position="none"),
                                          large.legal.F + theme(legend.position="none"),
                                          left=y.label,
                                          bottom=x.label,
                                          right=legend))
ggsave(LinePlots.F, filename="F_Combined_1987_new_start.png",height = a4.width*1, width = a4.width*1.1, units  ="mm", dpi = 300 )


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

Pop.Dist[[1]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S00", sep="_", "medium_movement"))
Pop.Dist[[2]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S01", sep="_", "medium_movement"))
Pop.Dist[[3]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S02", sep="_", "medium_movement"))
Pop.Dist[[4]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S03", sep="_", "medium_movement"))

Names <- c("Historical and Current NTZs", "Neither NTZs nor Temporal Management", 
           "Temporal and Spatial Management","Temporal Management Only")

Dists.res <- distance.abundance.format(Pops = Pop.Dist, max.year=59, n.sim=200, n.row=nrow(Whole_Pop_Ages), n.dist=3, 
                                       n.cell=NCELL, n.scenario=4, fished.cells=water, LQ=0.025, HQ=0.975, distances=Distances,
                                       dist.names=Dist.Names, scen.names=Names, start.year=26, start.year.mod=1960, end.year.mod=2018)

Pop.Dist.Median <- Dists.res[[1]]
Pop.Dist.Quant <- Dists.res[[2]]
Pop.Dist.Mean <- Dists.res[[3]]
Pop.Dist.SD <- Dists.res[[4]]


S00.means <- Pop.Dist.Mean[[1]]
S00.medians <- Pop.Dist.Median[[1]]
S00.SD <- Pop.Dist.SD[[1]]
S00.Quant <- Pop.Dist.Quant[[1]]
S00.data <- cbind(S00.means$Mean_Abundance, S00.SD$SD_Abundance, S00.medians$Median_Abundance, S00.Quant$`2.5%`, S00.Quant$`97.5%`) %>% 
  as.data.frame() %>% 
  rename(Mean_Abundance = "V1",
         SD_Abundance = "V2",
         Median_Abundance = "V3",
         Q2.5 = "V4",
         Q97.5 = "V5") %>% 
  mutate(Distance = S00.means$Distances,
         Year = S00.means$Year)


S00.abundance.dist.plot <- ggplot()+
  geom_line(data=S00.data, aes(x=Year, y=Median_Abundance, group=Distance, linetype=Distance),col="#36753B")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S00.data, aes(x=Year, y=Median_Abundance, ymin=Q2.5, ymax=Q97.5, group=Distance), fill="#36753B", alpha=0.2)+
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
  xlim(1986,2020)+
  ylim(0,8000)+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1987, y=7700, label="(b)", size = 2.5, fontface=1, hjust=0)
S00.abundance.dist.plot

S01.means <- Pop.Dist.Mean[[2]]
S01.SD <- Pop.Dist.SD[[2]]
S01.medians <- Pop.Dist.Median[[2]]
S01.Quant <- Pop.Dist.Quant[[2]]
S01.data <- cbind(S01.means$Mean_Abundance, S01.SD$SD_Abundance, S01.medians$Median_Abundance, S01.Quant$`2.5%`, S01.Quant$`97.5%`) %>% 
  as.data.frame() %>% 
  rename(Mean_Abundance = "V1",
         SD_Abundance = "V2",
         Median_Abundance = "V3",
         Q2.5 = "V4",
         Q97.5 = "V5") %>%
  mutate(Distances = S01.means$Distances,
         Year = S01.means$Year)

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
  xlim(1986,2020)+
  ylim(0,8000)+
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8))+
  theme(axis.text=element_text(size=8))+
  geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1987, y=7700, label="(d)", size = 2.5, fontface=1, hjust=0)
S01.abundance.dist.plot

S02.means <- Pop.Dist.Mean[[3]]
S02.SD <- Pop.Dist.SD[[3]]
S02.medians <- Pop.Dist.Median[[3]]
S02.Quant <- Pop.Dist.Quant[[3]]
S02.data <- cbind(S02.means$Mean_Abundance, S02.SD$SD_Abundance, S02.medians$Median_Abundance, S02.Quant$`2.5%`, S02.Quant$`97.5%`) %>% 
  as.data.frame() %>% 
  rename(Mean_Abundance = "V1",
         SD_Abundance = "V2",
         Median_Abundance = "V3",
         Q2.5 = "V4",
         Q97.5 = "V5") %>%
  mutate(Distances = S02.means$Distances,
         Year = S02.means$Year)

S02.abundance.dist.plot <- ggplot()+
  geom_line(data=S02.data, aes(x=Year, y=Median_Abundance, group=Distances, linetype=Distances),col="#66CCEE")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S02.data, aes(x=Year, y=Median_Abundance, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#66CCEE", alpha=0.2)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1986,2020)+
  ylim(0,8000)+
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
  ggplot2::annotate("text", x=1987, y=7700, label="(a)", size = 2.5, fontface=1, hjust=0)
S02.abundance.dist.plot


S03.means <- Pop.Dist.Mean[[4]]
S03.SD <- Pop.Dist.SD[[4]]
S03.medians <- Pop.Dist.Median[[4]]
S03.Quant <- Pop.Dist.Quant[[4]]
S03.data <- cbind(S03.means$Mean_Abundance, S03.SD$SD_Abundance, S03.medians$Median_Abundance, S03.Quant$`2.5%`, S03.Quant$`97.5%`) %>% 
  as.data.frame() %>% 
  rename(Mean_Abundance = "V1",
         SD_Abundance = "V2",
         Median_Abundance = "V3",
         Q2.5 = "V4",
         Q97.5 = "V5") %>%
  mutate(Distances = S03.means$Distances,
         Year = S03.means$Year)

S03.abundance.dist.plot <- ggplot()+
  geom_line(data=S03.data, aes(x=Year, y=Median_Abundance, group=Distances, linetype=Distances),col="#BBCC33")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"))+
  geom_ribbon(data=S03.data, aes(x=Year, y=Median_Abundance, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#BBCC33", alpha=0.2)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1986,2020)+
  ylim(0,8000)+
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
  ggplot2::annotate("text", x=1987, y=7700, label="(c)", size = 2.5, fontface=1, hjust=0)
S03.abundance.dist.plot

legend.plot <- ggplot()+
  geom_line(data=S03.data, aes(x=Year, y=Mean_Abundance, group=Distances, linetype=Distances),col="grey20")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  theme_classic()

## Put it all together
setwd(fig_dir)
x.label <- textGrob("Year", gp=gpar(fontsize=9))
y.label <- textGrob("Median abundance of greater than age at MLL (no. fish)", gp=gpar(fontsize=9), rot=90)
legend <- gtable_filter(ggplotGrob(legend.plot), "guide-box")

AbundancexDistance <-grid.arrange(arrangeGrob(S02.abundance.dist.plot + theme(legend.position="none"),
                                              S00.abundance.dist.plot + theme(legend.position="none"),
                                              S03.abundance.dist.plot + theme(legend.position="none"),
                                              S01.abundance.dist.plot + theme(legend.position="none"),
                                              left=y.label,
                                              bottom=x.label,
                                              right=legend))

ggsave(AbundancexDistance, filename="Abundance_Distance.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )


#* Catch ####
setwd(pop_dir)

Pop.Catch <- list()

# Each layer is a simulation, rows are cells and columns are years
Pop.Catch[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S00", sep="_", "medium_movement"))
Pop.Catch[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S01", sep="_", "medium_movement"))
Pop.Catch[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S02", sep="_", "medium_movement"))
Pop.Catch[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S03", sep="_", "medium_movement"))

Catch.res <- distance.catch.format(Pops = Pop.Catch, n.row=59, max.year=59, n.sim=200, n.dist=3, 
                                       n.cell=NCELL, n.scenario=4, fished.cells=water, LQ=0.025, HQ=0.975, distances=Distances,
                                       dist.names=Dist.Names, scen.names=Names, start.year=26, start.year.mod=1960, end.year.mod=2018)




Pop.Catch.Median <- Catch.res[[1]]
Pop.Catch.Quant <- Catch.res[[2]]
Pop.Catch.Mean <- Catch.res[[3]]
Pop.Catch.SD <- Catch.res[[4]]

S00.means <- Pop.Catch.Mean[[1]]
S00.medians <- Pop.Catch.Median[[1]]
S00.SD <- Pop.Catch.SD[[1]]
S00.Quant <- Pop.Catch.Quant[[1]]
S00.data <- cbind(S00.means$Mean_Catch, S00.SD$SD_Catch, S00.medians$Median_Catch, S00.Quant$`2.5%`, S00.Quant$`97.5%`) %>% 
  as.data.frame() %>% 
  rename(Mean_Catch = "V1",
         SD_Catch = "V2",
         Median_Catch = "V3",
         Q2.5 = "V4",
         Q97.5 = "V5") %>% 
  mutate(Distance = S00.means$Distances,
         Year = S00.means$Year) %>% 
  mutate(Q2.5 = ifelse(Year < 1988, Q2.5*0.7, Q2.5)) %>% 
  mutate(Q97.5 = ifelse(Year < 1988, Q97.5*1.5, Q97.5)) 


S00.catch.dist.plot <- ggplot()+
  geom_line(data=S00.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="#36753B")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S00.data, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#36753B", alpha=0.2)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1986, 2020)+
  ylim(0,1000)+
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
  ggplot2::annotate("text", x=1988, y=1000, label="(a)", size = 2.5, fontface=1, hjust=0)
S00.catch.dist.plot


S01.means <- Pop.Catch.Mean[[2]]
S01.medians <- Pop.Catch.Median[[2]]
S01.SD <- Pop.Catch.SD[[2]]
S01.Quant <- Pop.Catch.Quant[[2]]
S01.data <- cbind(S01.means$Mean_Catch, S01.SD$SD_Catch, S01.medians$Median_Catch, S01.Quant$`2.5%`, S01.Quant$`97.5%`) %>% 
  as.data.frame() %>% 
  rename(Mean_Catch = "V1",
         SD_Catch = "V2",
         Median_Catch = "V3",
         Q2.5 = "V4",
         Q97.5 = "V5") %>% 
  mutate(Distance = S01.means$Distances,
         Year = S01.means$Year) %>% 
  mutate(Q2.5 = ifelse(Year < 1988, Q2.5*0.7, Q2.5)) %>% 
  mutate(Q97.5 = ifelse(Year < 1988, Q97.5*1.5, Q97.5)) 


S01.catch.dist.plot <- ggplot()+
  geom_line(data=S01.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="#302383")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S01.data, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#302383", alpha=0.2)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1986,2020)+
  ylim(0,1000)+
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
  ggplot2::annotate("text", x=1988, y=1000, label="(b)", size = 2.5, fontface=1, hjust=0)
S01.catch.dist.plot

S02.means <- Pop.Catch.Mean[[3]]
S02.medians <- Pop.Catch.Median[[3]]
S02.SD <- Pop.Catch.SD[[3]]
S02.Quant <- Pop.Catch.Quant[[3]]
S02.data <- cbind(S02.means$Mean_Catch, S02.SD$SD_Catch, S02.medians$Median_Catch, S02.Quant$`2.5%`, S02.Quant$`97.5%`) %>% 
  as.data.frame() %>% 
  rename(Mean_Catch = "V1",
         SD_Catch = "V2",
         Median_Catch = "V3",
         Q2.5 = "V4",
         Q97.5 = "V5") %>% 
  mutate(Distance = S02.means$Distances,
         Year = S02.means$Year) %>% 
  mutate(Q2.5 = ifelse(Year < 1988, Q2.5*0.7, Q2.5)) %>% 
  mutate(Q97.5 = ifelse(Year < 1988, Q97.5*1.5, Q97.5)) 

S02.catch.dist.plot <- ggplot()+
  geom_line(data=S02.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="#66CCEE")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S02.data, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#66CCEE", alpha=0.1)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1986, 2020)+
  ylim(0,1000)+
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
  ggplot2::annotate("text", x=1988, y=1000, label="(c)", size = 2.5, fontface=1, hjust=0)
S02.catch.dist.plot


S03.means <- Pop.Catch.Mean[[4]]
S03.medians <- Pop.Catch.Median[[4]]
S03.SD <- Pop.Catch.SD[[4]]
S03.Quant <- Pop.Catch.Quant[[4]]
S03.data <- cbind(S03.means$Mean_Catch, S03.SD$SD_Catch, S03.medians$Median_Catch, S03.Quant$`2.5%`, S03.Quant$`97.5%`) %>% 
  as.data.frame() %>% 
  rename(Mean_Catch = "V1",
         SD_Catch = "V2",
         Median_Catch = "V3",
         Q2.5 = "V4",
         Q97.5 = "V5") %>% 
  mutate(Distance = S03.means$Distances,
         Year = S03.means$Year) %>% 
  mutate(Q2.5 = ifelse(Year < 1988, Q2.5*0.7, Q2.5)) %>% 
  mutate(Q97.5 = ifelse(Year < 1988, Q97.5*1.5, Q97.5)) 


S03.catch.dist.plot <- ggplot()+
  geom_line(data=S03.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="#BBCC33")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"))+
  geom_ribbon(data=S03.data, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#BBCC33", alpha=0.1)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1986, 2020)+
  ylim(0,1000)+
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
  ggplot2::annotate("text", x=1988, y=1000, label="(d)", size = 2.5, fontface=1, hjust=0)
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

abundance.0_10 <- NULL
catch.0_10 <- NULL

for (S in 1:4){
  temp.mean <- Pop.Dist.Mean[[S]]
  temp.SD <- Pop.Dist.SD[[S]]
  temp.median <- Pop.Dist.Median[[S]]
  temp.quant <- Pop.Dist.Quant[[S]]
  
  temp.all <- cbind(temp.mean, temp.SD$SD_Abundance, temp.median$Median_Abundance, temp.quant$`2.5%`, temp.quant$`97.5%`) %>% 
    as.data.frame() %>% 
    rename(SD_Abundance = "temp.SD$SD_Abundance",
           Median_Abundance = "temp.median$Median_Abundance",
           Q2.5 = "temp.quant$`2.5%`",
           Q97.5 = "temp.quant$`97.5%`") %>% 
    filter(Distances %in% c("0-10 km"))
  
  abundance.0_10 <- rbind(abundance.0_10, temp.all)
  
  temp.mean <- Pop.Catch.Mean[[S]]
  temp.SD <- Pop.Catch.SD[[S]]
  temp.median <- Pop.Catch.Median[[S]]
  temp.quant <- Pop.Catch.Quant[[S]]
  
  temp.all <- cbind(temp.mean, temp.SD$SD_Catch, temp.median$Median_Catch, temp.quant$`2.5%`, temp.quant$`97.5%`) %>% 
    as.data.frame() %>% 
    rename(SD_Catch = "temp.SD$SD_Catch",
           Median_Catch = "temp.median$Median_Catch",
           Q2.5 = "temp.quant$`2.5%`",
           Q97.5 = "temp.quant$`97.5%`") %>% 
    filter(Distances %in% c("0-10 km"))
  
  catch.0_10 <- rbind(catch.0_10, temp.all)
}

abundance.0_10.plot <- abundance.0_10 %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs"), "Historical and\ncurrent NTZs", 
                                                                 ifelse(Scenario %in% c("Temporal Management Only"), "Temporal\nmanagement only", 
                                                                        ifelse(Scenario %in% c("Neither NTZs nor Temporal Management"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))%>% 
  filter(Year>1985) %>%
  # mutate(Q2.5 = ifelse(Year < 1988, Q2.5*0.1, Q2.5)) %>%
  # mutate(Q97.5 = ifelse(Year < 1988, Q97.5*1.2, Q97.5)) %>%
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
  ggplot2::annotate("text", x=1987, y=2000, label="(a)", size = 2.5, fontface=1)
abundance.0_10.plot


Catch.Pre_1987 <- catch.0_10 %>% 
  filter(Year>1984) %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs") & Year>1985, "Historical and\ncurrent NTZs", 
                                                             ifelse(Scenario %in% c("Temporal Management Only") & Year>1985, "Temporal\nmanagement only", 
                                                                    ifelse(Scenario %in% c("Neither NTZs nor Temporal Management"), "Neither NTZs nor\ntemporal management", "NTZs and\ntemporal management")))))
Catch.10_50.plot <- catch.0_10%>% 
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
  ggplot2::annotate("text", x=1987, y=300, label="(b)", size = 2.5)
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



#### CPUE PLOTS FOR DISTANCE FROM BOAT RAMP ####

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


## 0-10km all scenarios
cpue.0_10 <- catch.0_10 %>% 
  mutate(Effort = Effort_Dist[[1]]$Effort) %>% 
  mutate(Median_CPUE = Median_Catch/Effort,
         CPUE_2.5 = Q2.5/Effort,
         CPUE_97.5 = Q97.5/Effort)


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
  theme_classic()+
  xlab(NULL)+
  ylab(bquote(CPUE~(Fish~Boat~days^-1~~km^-1)))+
  xlim(1987,2020)+
  #ylim(0,0.08)+
  scale_linetype_manual(values = c("longdash", "solid" ), labels=c("Always Fished", "NTZ Area"), name="Model Area")+
  theme(legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=8), #change legend text font size
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(2,"line")) +
  guides(color = guide_legend(byrow = TRUE))+
  theme(axis.text=element_text(size=8),
        axis.title.y = element_text(size=9))+
  geom_vline(xintercept=1987, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  geom_vline(xintercept=2017, linetype="dotted", colour="grey20")+
  ggplot2::annotate("text", x=1988, y=0.8, label="(a)", size = 2.5)
CPUE.0_10.plot


## Save the plot
setwd(fig_dir)
x.label <- textGrob("Year", gp=gpar(fontsize=9))
y.label <- textGrob(bquote(CPUE~(Fish~Boat~days^-1~~km^-1)),gp=gpar(fontsize=9), rot=90)
legend <- gtable_filter(ggplotGrob(Catch.10_50.plot), "guide-box")

CPUE.AbundancexDistance <-grid.arrange(arrangeGrob(abundance.10_50.plot + theme(legend.position="none"),
                                                   CPUE.0_10.plot + theme(legend.position="none"),
                                                    bottom=x.label,
                                                    right=legend))

ggsave(CPUE.AbundancexDistance, filename="CPUE_Distance_Combined.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

CPUE.Distance.Scenarios <-grid.arrange(arrangeGrob(CPUE.0_10.plot + theme(legend.position="none")+ylab(NULL),
                                                   CPUE.10_50.plot + theme(legend.position="none")+ylab(NULL),
                                                   CPUE.50_100.plot + theme(legend.position="none")+ylab(NULL),
                                                   bottom=x.label,
                                                   left=y.label,
                                                   right=legend))

ggsave(CPUE.Distance.Scenarios, filename="CPUE_Distance_All_Scenarios.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

#### CPUE PLOTS FOR ALL DISTANCES BY SCENARIO ####

setwd(pop_dir)

Pop.Catch <- list()

# Each layer is a simulation, rows are cells and columns are years
Pop.Catch[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S00", sep="_", "medium_movement"))
Pop.Catch[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S01", sep="_", "medium_movement"))
Pop.Catch[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S02", sep="_", "medium_movement"))
Pop.Catch[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S03", sep="_", "medium_movement"))


for(S in 1:4){
  
  Scenario <- Pop.Catch[[S]]

  Medians <- array(0, dim=c(3, 59))%>% 
    as.data.frame(.)
  Quantiles <- array(0, dim=c(6, 59))%>% 
    as.data.frame(.)
  temp <- array(0, dim=c(1834,59,200))
  
  for(SIM in 1:200){
    temp[,,SIM] <- Scenario[[SIM]]
  }
  for(YEAR in 26:59){
    
    temp1 <- temp[,YEAR,]
    tempE <- Effort_Dist[[1]] %>% 
      filter(Scenario %in% Names[S])
    
    temp1.1 <- temp1/tempE$Effort
    
    temp2 <- as.data.frame(temp1.1) %>% 
      mutate(CellID = row_number()) %>% 
      mutate(Fished_17 = water$Fished_2017) 
    
    
    Dist.10km <- temp2 %>% 
      filter(CellID %in% c(Distances[[1]])) %>% 
      #{if (S>0) filter(., Fished_17 =="Y") else .} %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      # mutate(Mean_Pop = rowMeans(.[1:200])) %>% 
      # mutate(SD_Pop = rowSds(as.matrix(.[,1:200]))) %>% 
      mutate(Median_Pop = rowMedians(as.matrix(.[1:200]))) %>%
      mutate(Distance = "10km") 
    
    ## 10-50km
    tempE <- Effort_Dist[[2]] %>% 
      filter(Scenario %in% Names[S])
    
    temp1.2 <- temp1/tempE$Effort
    
    temp2 <- as.data.frame(temp1.2) %>% 
      mutate(CellID = row_number()) %>% 
      mutate(Fished_17 = water$Fished_2017)
    
    Dist.50km <- temp2 %>% 
      filter(CellID %in% c(Distances[[2]])) %>% 
      #{if (S>0) filter(., Fished_17 =="Y") else .} %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      # mutate(Mean_Pop = rowMeans(.[1:200])) %>% 
      # mutate(SD_Pop = rowSds(as.matrix(.[,1:200]))) %>% 
      mutate(Median_Pop = rowMedians(as.matrix(.[1:200]))) %>%
      mutate(Distance = "50km")
    
    ## 50-100km
    
    tempE <- Effort_Dist[[3]] %>% 
      filter(Scenario %in% Names[S])
    
    temp1.3 <- temp1/tempE$Effort
    
    temp2 <- as.data.frame(temp1.3) %>% 
      mutate(CellID = row_number()) %>% 
      mutate(Fished_17 = water$Fished_2017)
    
    Dist.100km <- temp2 %>% 
      filter(CellID %in% c(Distances[[3]])) %>% 
      #{if (S>0) filter(., Fished_17 =="Y") else .} %>% 
      summarise(across(where(is.numeric), sum)) %>% 
      # mutate(Mean_Pop = rowMeans(.[1:200])) %>% 
      # mutate(SD_Pop = rowSds(as.matrix(.[,1:200]))) %>% 
      mutate(Median_Pop = rowMedians(as.matrix(.[1:200]))) %>%
      mutate(Distance = "100km") 
    
    temp2 <- as.matrix(Dist.10km[ , 1:200])
    
    Catch.Quantiles.10 <- quantile(temp2, probs=c(0.025, 0.975)) %>% 
      as.data.frame() %>% 
      mutate(Distance = "10 km") 
    
    temp2 <- as.matrix(Dist.50km[ , 1:200])
    
    Catch.Quantiles.50 <-quantile(temp2, probs=c(0.025, 0.975)) %>% 
      as.data.frame() %>% 
      mutate(Distance = "50 km") 
    
    temp2 <- as.matrix(Dist.100km[ , 1:200])
    
    Catch.Quantiles.100 <-quantile(temp2, probs=c(0.025, 0.975)) %>% 
      as.data.frame() %>% 
      mutate(Distance = "100 km") 
    
    # Means.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>% 
    #   dplyr::select(Distance, Median_Pop)
    # SDs.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>% 
    #   dplyr::select(Distance, SD_Pop)
    Medians.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>% 
      dplyr::select(Distance, Median_Pop)
    Quantiles.Full <- rbind(Catch.Quantiles.10, Catch.Quantiles.50, Catch.Quantiles.100) %>%
      as.data.frame() 
    Quantiles.Full$Quantile <- rownames(Quantiles.Full)
    
    # Means[,YEAR] <- Means.Full$Mean_Pop
    # 
    # SDs[ ,YEAR] <- SDs.Full$SD_Pop
    
    Medians[ ,YEAR] <- Medians.Full$Median_Pop
    
    Quantiles[, YEAR] <- Quantiles.Full$. 
    Quantiles$Quantile <- Quantiles.Full$Quantile
    
    
  }
  
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
  
  # Pop.Catch.Mean[[S]] <- Means
  # Pop.Catch.SD[[S]] <- SDs
  Pop.Catch.Median[[S]] <- Medians
  Pop.Catch.Quant[[S]] <- Quantiles
  
}

# S00.means <- Pop.Catch.Mean[[1]]
S00.medians <- Pop.Catch.Median[[1]]
# S00.SD <- Pop.Catch.SD[[1]]
S00.Quant <- Pop.Catch.Quant[[1]]
S00.data <- cbind(S00.medians, S00.Quant$`2.5%`, S00.Quant$`97.5%`) %>% 
  rename(Q2.5 = "S00.Quant$`2.5%`",
         Q97.5 = "S00.Quant$`97.5%`")%>% 
  mutate(Q2.5 = ifelse(Year < 1988, Q2.5*0.7, Q2.5)) %>% 
  mutate(Q97.5 = ifelse(Year < 1988, Q97.5*1.5, Q97.5)) 


S00.CPUE.dist.plot <- ggplot()+
  geom_line(data=S00.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="#36753B")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S00.data, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#36753B", alpha=0.2)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1986, 2020)+
  ylim(0,2)+
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
  ggplot2::annotate("text", x=1988, y=2, label="(b)", size = 2.5, fontface=1, hjust=0)
S00.CPUE.dist.plot


# S01.means <- Pop.Catch.Mean[[2]]
S01.medians <- Pop.Catch.Median[[2]]
# S01.SD <- Pop.Catch.SD[[2]]
S01.Quant <- Pop.Catch.Quant[[2]]
S01.data <- cbind(S01.medians, S01.Quant$`2.5%`, S01.Quant$`97.5%`) %>% 
  rename(Q2.5 = "S01.Quant$`2.5%`",
         Q97.5 = "S01.Quant$`97.5%`") %>% 
  mutate(Q2.5 = ifelse(Year < 1988, Q2.5*0.7, Q2.5)) %>% 
  mutate(Q97.5 = ifelse(Year < 1988, Q97.5*1.5, Q97.5)) 

S01.CPUE.dist.plot <- ggplot()+
  geom_line(data=S01.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="#302383")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S01.data, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#302383", alpha=0.2)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1986,2020)+
  ylim(0,2)+
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
  ggplot2::annotate("text", x=1988, y=2, label="(d)", size = 2.5, fontface=1, hjust=0)
S01.CPUE.dist.plot

# S02.means <- Pop.Catch.Mean[[3]]
S02.medians <- Pop.Catch.Median[[3]]
# S02.SD <- Pop.Catch.SD[[3]]
S02.Quant <- Pop.Catch.Quant[[3]]
S02.data <- cbind(S02.medians, S02.Quant$`2.5%`, S02.Quant$`97.5%`) %>% 
  rename(Q2.5 = "S02.Quant$`2.5%`",
         Q97.5 = "S02.Quant$`97.5%`") %>% 
  mutate(Q2.5 = ifelse(Year < 1988, Q2.5*0.7, Q2.5)) %>% 
  mutate(Q97.5 = ifelse(Year < 1988, Q97.5*1.5, Q97.5)) 

S02.CPUE.dist.plot <- ggplot()+
  geom_line(data=S02.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="#66CCEE")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  geom_ribbon(data=S02.data, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#66CCEE", alpha=0.2)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1986, 2020)+
  ylim(0,2)+
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
  ggplot2::annotate("text", x=1988, y=2, label="(a)", size = 2.5, fontface=1, hjust=0)
S02.CPUE.dist.plot


# S03.means <- Pop.Catch.Mean[[4]]
S03.medians <- Pop.Catch.Median[[4]]
# S03.SD <- Pop.Catch.SD[[4]]
S03.Quant <- Pop.Catch.Quant[[4]]
S03.data <- cbind(S03.medians, S03.Quant$`2.5%`, S03.Quant$`97.5%`) %>% 
  rename(Q2.5 = "S03.Quant$`2.5%`",
         Q97.5 = "S03.Quant$`97.5%`") %>% 
  mutate(Q2.5 = ifelse(Year < 1988, Q2.5*0.7, Q2.5)) %>% 
  mutate(Q97.5 = ifelse(Year < 1988, Q97.5*1.5, Q97.5)) 

S03.CPUE.dist.plot <- ggplot()+
  geom_line(data=S03.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="#BBCC33")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"))+
  geom_ribbon(data=S03.data, aes(x=Year, y=Median_Catch, ymin=Q2.5, ymax=Q97.5, group=Distances), fill="#BBCC33", alpha=0.2)+
  theme_classic()+
  ylab(NULL)+
  xlab(NULL)+
  xlim(1986, 2020)+
  ylim(0,2)+
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
  ggplot2::annotate("text", x=1988, y=2, label="(c)", size = 2.5, fontface=1, hjust=0)
S03.CPUE.dist.plot

legend.plot <- ggplot()+
  geom_line(data=S03.data, aes(x=Year, y=Median_Catch, group=Distances, linetype=Distances),col="grey20")+
  scale_linetype_manual(values=c("0-10 km"="solid", "10-50 km"="dashed", "50-100 km" = "dotted"), name="Distance from\nboat ramp")+
  theme_classic()

## Put it all together
setwd(fig_dir)
x.label <- textGrob("Year", gp=gpar(fontsize=9))
y.label <- ylab(bquote(CPUE~(Fish~per~boat~days~per~km^2)))
legend <- gtable_filter(ggplotGrob(legend.plot), "guide-box")

CPUExDistance <-grid.arrange(arrangeGrob(S02.CPUE.dist.plot + theme(legend.position="none"),
                                          S00.CPUE.dist.plot + theme(legend.position="none"),
                                          S03.CPUE.dist.plot + theme(legend.position="none"),
                                          S01.CPUE.dist.plot + theme(legend.position="none"),
                                          left=y.label,
                                          bottom=x.label,
                                          right=legend))

ggsave(CPUExDistance, filename="CPUE_Distance.png",height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

#### CATCH PER UNIT EFFORT PLOTS DUPLICATE? ####


# #* Format catch data and add to effort ####
# setwd(pop_dir)
# 
# Pop.Catch <- list()
# 
# # Each layer is a simulation, rows are cells and columns are years
# Pop.Catch[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S00"))
# Pop.Catch[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S01"))
# Pop.Catch[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S02"))
# Pop.Catch[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S03"))
# 
# temp <- array(0, dim=c(NCELL,59,100))
# Pop.Catch.Mean <- NULL
# Pop.Catch.SD <- NULL
# 
# 
# for(S in 1:4){
#   
#   Scenario <- Pop.Catch[[S]]
#   
#   Means <- array(0, dim=c(1, 59)) %>% 
#     as.data.frame(.)
#   SDs <- array(0, dim=c(1, 59))%>% 
#     as.data.frame(.)
#   
#   for(SIM in 1:100){
#     temp[,,SIM] <- Scenario[[SIM]]
#   }
#   for(YEAR in 1:59){
#     
#     temp1 <- temp[,YEAR,]
#     
#     temp2 <- temp1 %>% 
#       as.data.frame() %>% 
#       summarise(across(where(is.numeric), sum)) %>% 
#       mutate(Median_Catch = rowMedians(.[1:100])) %>% 
#       mutate(SD_Catch = rowSds(as.matrix(.[,1:100]))) 
#     
#     Means[,YEAR] <- temp2$Mean_Catch
#     SDs[,YEAR] <- temp2$SD_Catch
#   }
#   Means <- as.data.frame(Means) %>% 
#     mutate(Scenario = Names[S]) %>% 
#     pivot_longer(cols=-c("Scenario"), values_to = "Mean_Catch", names_to = "Year") %>% 
#     mutate(Year = seq(1960,2018,1))
#   
#   SDs <- as.data.frame(SDs) %>% 
#     mutate(Scenario = Names[S]) %>% 
#     pivot_longer(cols=-c("Scenario"), values_to = "SD_Catch", names_to = "Year") %>% 
#     mutate(Year = seq(1960,2018,1))
#   
#   Pop.Catch.Mean <- rbind(Pop.Catch.Mean, Means)
#   Pop.Catch.SD <- rbind(Pop.Catch.SD, SDs)
#   
# }
# 
# CPUE <- Boat_Days_sum %>% 
#   mutate(Catch = Pop.Catch.Mean$Mean_Catch) %>% 
#   mutate(SD.Catch = Pop.Catch.SD$SD_Catch) %>% 
#   mutate(Year = rep(seq(1960,2018), 4)) %>% 
#   mutate(cpue = Catch/Effort) %>% 
#   mutate(cpue = ifelse(is.nan(cpue), 0, cpue))
# 
# #* Plot CPUE for each scenario ####
# 
# CPUE.Pre_1987 <- CPUE %>% 
#   filter(Year<=1986) %>%  
#   mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs") & Year>1985, "Historical and\ncurrent management", 
#                                                              ifelse(Scenario %in% c("Temporal Management Only") & Year>1985, "Temporal\nmanagement only", 
#                                                                     ifelse(Scenario %in% c("No NTZs or Temporal Management"), "No NTZs or\ntemporal management", "NTZs and\ntemporal management")))))
# 
# CPUE.plot <- CPUE %>% 
#   mutate(ColourGroup = ifelse(Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Historical and Current NTZs"), "Historical and\ncurrent management", 
#                                                              ifelse(Scenario %in% c("Temporal Management Only"), "Temporal\nmanagement only", 
#                                                                     ifelse(Scenario %in% c("No NTZs or Temporal Management"), "No NTZs or\ntemporal management", "NTZs and\ntemporal management")))))%>% 
#   filter(Year>1985) %>% 
#   ggplot(.)+
#   geom_line(aes(x=Year, y=cpue, group=Scenario, colour=ColourGroup), size=0.7)+
#   scale_fill_manual(values= c("Historical and\ncurrent management"="#36753B", "No NTZs or\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
#                               "Temporal\nmanagement only"="#BBCC33"),
#                     guide="none")+
#   scale_colour_manual(values = c( "Historical and\ncurrent management"="#36753B", "No NTZs or\ntemporal management"="#302383" ,"NTZs and\ntemporal management"="#66CCEE",
#                                   "Temporal\nmanagement only"="#BBCC33"), name= "Spatial and Temporal\nManagement Scenario")+ 
#   geom_line(data=CPUE.Pre_1987, aes(x=Year, y=cpue, group=Scenario, colour="grey20"), size=0.7)+
#   theme_classic()+
#   xlab(NULL)+
#   ylab("CPUE")+
#   xlim(1970,2020)+
#   ylim(0,0.3)+
#   #scale_linetype_manual(values = c("longdash", "solid" ), labels=c("Always Fished", "NTZ Area"), name="Model Area")+
#   theme(legend.title = element_text(size=9, face="bold"), #change legend title font size
#         legend.text = element_text(size=8), #change legend text font size
#         legend.spacing.y = unit(0.1, "cm"),
#         legend.key.size = unit(2,"line")) +
#   guides(color = guide_legend(byrow = TRUE))+
#   theme(axis.text=element_text(size=8),
#         axis.title = element_text(size=9))+
#   geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
#   geom_vline(xintercept=2005, colour="grey20")+
#   geom_vline(xintercept=2017, linetype="dotted", colour="grey20")
#   #ggplot2::annotate("text", x=1960, y=400, label="(b)", size = 2.5, fontface=2)
# CPUE.plot


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






















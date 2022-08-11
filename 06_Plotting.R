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

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
m_dir <- paste(working.dir, "Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")
pop_dir <-  paste(working.dir, "Output_Population", sep="/")


# Normal
# Sim 1 Nothing
# Sim 2 NTZs and Temporal Closure
# Sim 3 Just temporal closure, no sanctuary zones 


#### READ IN DATA ####
setwd(sg_dir)
NoTake <- readRDS("NoTakeList")

#* Scenario 1 Normal - IN PROGRESS ####
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
  
  year <- readRDS(paste0("YearlyTotal.", numYear[Y]))
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
  
  Population <-  readRDS(paste0("YearlyTotal.", numYears[YEAR]))
  
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

#* Scenario 2 - Nothing COMPLETED ####
#* Scenario 3 - NTZ and Temp Closure COMPLETED ####
#* Scenario 4 - Just Temp Closure COMPLETED ####

#### CALCULATE TOTAL AREA OF FISHED AND NO-TAKE ####
AreaFished <- water %>% 
  mutate(cell_area = as.numeric(cell_area)) %>% 
  mutate(Fished=as.factor(Fished)) %>% 
  filter(Fished=="Y") 

AreaFished <- (sum(AreaFished$cell_area))/100000

AreaNT <- water %>% 
  mutate(cell_area = as.numeric(cell_area)) %>% 
  mutate(Fished=as.factor(Fished)) %>% 
  filter(Fished=="N") 

AreaNT <- (sum(AreaNT$cell_area))/100000

#### TIME SERIES PLOT WITH ALL SCENARIOS ####

TotalPop <- TotalPop %>% 
  mutate(sTot.Pop = sTotalPop$sTot.Pop) %>% 
  rename(Normal = "Tot.Pop") %>% 
  rename(Simulation = "sTot.Pop") %>% 
  pivot_longer(cols=-c(Year), names_to="Scenario", values_to="Total.Population")

TimeSeries <- TotalPop %>% 
  mutate(ColourGroup = ifelse(Year<=1985, "Pre 1987", ifelse(Scenario %in% c("Normal") & Year>1985, "With NTZs", ifelse(Scenario %in% ("Simulation") & Year>1985, "Without NTZs", 0)))) %>% 
  mutate(ColourGroup = as.factor(ColourGroup)) %>% 
  mutate(ColourGroup = fct_relevel(ColourGroup, c("Pre 1987", "With NTZs", "Without NTZs"))) %>% 
  ggplot()+
  geom_line(aes(x=Year, y=Total.Population, group=Scenario, colour=ColourGroup)) +
  scale_x_continuous("Year", breaks = c(1960, 1970, 1980, 1990, 2000, 2010, 2020))+
  scale_colour_manual("Scenario",values=c( "gray20",  "#7ac0c2", "#487088"), labels=c("Pre 1987", 
                                                                                      "With NTZs in place",
                                                                                      "Without NTZs in place"))+
  xlab("Year")+
  ylab("Total Population")+
  geom_vline(xintercept=1987, linetype="dashed", color="grey20")+
  geom_vline(xintercept=2005, colour="grey20")+
  theme_classic()

#### LINE PLOTS BY AGE GROUP OF ALL SCENARIOS ####

AllNoTake <- rbind(NoTakeAges, sNoTakeAges) 
AllFished <- rbind(FishedAges, sFishedAges)

ScenarioWholePop <- rbind(AllNoTake, AllFished) %>%  
  pivot_longer(cols=-c(Age, Status, Scenario), names_to="Year", values_to="Number") %>% 
  mutate(Year = as.factor(Year)) %>% 
  mutate(NumKM2 = ifelse(Status %in% c("Fished"), Number/AreaFished, Number/AreaNT))

ScenarioNoRecruits <- ScenarioWholePop %>% 
  filter(Age!=1) 

line.recruits <- ScenarioWholePop %>% 
  filter(Age==1) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) %>% 
  ggplot(., aes(x=Year, y=NumKM2, group=interaction(Status,Scenario), colour=NTGroup, shape=interaction(Status, Scenario)))+
  geom_point(size=2.5)+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=9, colour="grey20")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="grey20")+
  scale_shape_manual(values= c(16,15,1,0), name="Area and scenario", 
                     labels= c("Fished area with NTZs in place", "NTZ area with NTZs in place", 
                               "Fished area with no NTZs in place", "NTZ area with no NTZs in place"))+
  scale_colour_manual(values = c("gray20", "#7ac0c2", "#47315e", "#487088"), guide=FALSE)+
  ylab(NULL)+
  xlab(NULL)+
  ggplot2::annotate("text", x=1.3, y=1.1, label="(a) Recruits", size = 3, fontface=2)

line.sublegal <- ScenarioWholePop %>% 
  filter(Age>1 & Age <=4) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) %>% 
  ggplot(., aes(x=Year, y=Total, group=interaction(Status,Scenario), colour=NTGroup, shape=interaction(Status, Scenario)))+
  geom_point(size=2.5)+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=9, colour="grey20")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="grey20")+
  scale_shape_manual(values= c(16,15,1,0), name="Area and scenario", 
                     labels= c("Fished area with NTZs in place", "NTZ area with NTZs in place", 
                               "Fished area with no NTZs in place", "NTZ area with no NTZs in place"))+
  scale_colour_manual(values = c("gray20", "#7ac0c2", "#47315e", "#487088"), guide=FALSE)+
  ylab(NULL)+
  xlab(NULL)+
  ggplot2::annotate("text", x=1.4, y=0.3, label="(b) Sub-legal", size = 3, fontface=2)

line.legal <- ScenarioWholePop %>% 
  filter(Age>4 & Age <=10) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) %>% 
  ggplot(., aes(x=Year, y=Total, group=interaction(Status, Scenario), colour=NTGroup, shape=interaction(Status, Scenario)))+
  geom_point(size=2.5)+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=9, colour="grey20")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="grey20")+
  scale_shape_manual(values= c(16,15,1,0), name="Area and scenario", 
                     labels= c("Fished area with NTZs in place", "NTZ area with NTZs in place", 
                               "Fished area with no NTZs in place", "NTZ area with no NTZs in place"))+
  scale_colour_manual(values = c("gray20", "#7ac0c2", "#47315e", "#487088"), guide=FALSE)+
  ylab(NULL)+
  xlab(NULL)+
  ggplot2::annotate("text", x=1.1, y=0.32, label="(c) Legal", size = 3, fontface=2)

line.biglegal <- ScenarioWholePop %>% 
  filter(Age>10) %>% 
  group_by(Scenario, Year, Status) %>% 
  mutate(Total = sum(NumKM2)) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) %>% 
  ggplot(., aes(x=Year, y=Total, group=interaction(Status,Scenario), colour=NTGroup, shape=interaction(Status,Scenario)))+
  geom_point(size=2.5)+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="grey20")+
  geom_vline(xintercept=9, colour="grey20")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="grey20")+
  scale_shape_manual(values= c(16,15,1,0), name="Area and scenario", 
                     labels= c("Fished area with NTZs in place", "NTZ area with NTZs in place", 
                               "Fished area with no NTZs in place", "NTZ area with no NTZs in place"))+
  scale_colour_manual(values = c("gray20", "#7ac0c2", "#47315e", "#487088"), guide=FALSE)+
  ylab(NULL)+
  xlab(NULL)+
  ggplot2::annotate("text", x=1.5, y=0.22, label="(d) Large Legal", size = 3, fontface=2)

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



#### PUT PLOTS TOGETHER FOR PUBLISHING ####
x.label <- textGrob("Year", gp=gpar(fontsize=14))
y.label <- textGrob("No. Fish per"~km^2, gp=gpar(fontsize=14), rot=90)
legend <- gtable_filter(ggplotGrob(line.biglegal), "guide-box")

LinePlotsxGroup <-grid.arrange(arrangeGrob(line.recruits + theme(legend.position="none"),
                                           line.sublegal + theme(legend.position="none"),
                                           line.legal + theme(legend.position="none"),
                                           line.biglegal + theme(legend.position="none"),
                                           nrow = 2,
                                           left = y.label,
                                           bottom = x.label), 
                               legend, 
                               widths=unit.c(unit(1, "npc") - legend$width, legend$width), 
                               nrow=1)




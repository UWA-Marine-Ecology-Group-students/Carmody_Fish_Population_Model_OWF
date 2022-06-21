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

#### READ IN DATA ####
setwd(pop_dir)

Year5 <- readRDS("YearlyTotal.5") 
Year5.NT <- Year5[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year5.F <- Year5[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

Year10 <- readRDS("YearlyTotal.10")
Year10.NT <- Year10[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year10.F <- Year10[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

Year15 <- readRDS("YearlyTotal.15") 
Year15.NT <- Year15[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year15.F <- Year15[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

Year20 <- readRDS("YearlyTotal.20")
Year20.NT <- Year20[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year20.F <- Year20[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

Year25 <- readRDS("YearlyTotal.25")
Year25.NT <- Year25[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year25.F <- Year25[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

# First sanctuaries at year 27
Year30 <- readRDS("YearlyTotal.30")
Year30.NT <- Year30[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year30.F <- Year30[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)
  
Year35 <- readRDS("YearlyTotal.35")
Year35.NT <- Year35[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year35.F <- Year35[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

Year40 <- readRDS("YearlyTotal.40")
Year40.NT <- Year40[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year40.F <- Year40[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

Year45 <- readRDS("YearlyTotal.45")
Year45.NT <- Year45[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year45.F <- Year45[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

## Sanctuaries expanded
Year50 <- readRDS("YearlyTotal.50")
Year50.NT <- Year50[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year50.F <- Year50[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

## Commonwealth sanctuary put in at year 53
Year55 <- readRDS("YearlyTotal.55")
Year55.NT <- Year55[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year55.F <- Year55[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

Year59 <- readRDS("YearlyTotal.59")
Year59.NT <- Year59[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year59.F <- Year59[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

#### CREATE POP TOTAL ####

NoTakeAges <- as.data.frame(cbind(Year5.NT, Year10.NT, Year15.NT, Year20.NT, Year25.NT, Year30.NT, Year35.NT, Year40.NT, Year45.NT, Year50.NT, Year55.NT, Year59.NT)) %>% 
  mutate(Status = "NTZ") %>% 
  mutate(Age = seq(1:30)) %>% 
  mutate(Scenario = "With NTZ") %>% 
  mutate(Scenario = as.factor(Scenario))
colnames(NoTakeAges) <- c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018", "Status", "Age", "Scenario")


FishedAges <- as.data.frame(cbind(Year5.F, Year10.F, Year15.F, Year20.F, Year25.F, Year30.F, Year35.F, Year40.F, Year45.F, Year50.F, Year55.F, Year59.F)) %>% 
  mutate(Status = "Fished") %>% 
  mutate(Age = seq(1:30)) %>% 
  mutate(Scenario = "With NTZ") %>% 
  mutate(Scenario = as.factor(Scenario))
colnames(FishedAges) <- c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018", "Status", "Age", "Scenario")

WholePop <- rbind(NoTakeAges, FishedAges) %>%  
  pivot_longer(cols=-c(Age, Status, Scenario), names_to="Year", values_to="Number") %>% 
  mutate(Year = as.factor(Year))

WholePop <- WholePop %>% 
  mutate(Year = fct_relevel(Year, c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018"))) %>% 
  mutate(Status = as.factor(Status)) 
  
NoRecruits <- WholePop %>% 
  filter(Age!=1) 

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

#### RIDGELINE PLOTS WITHOUT RECRUITS ####

ridge.norec.whole <- NoRecruits %>% 
  group_by(Age, Year) %>% 
  mutate(Tot.Num = sum(Number)) %>% 
  ungroup() %>% 
  ggplot(., aes(x=Age, y=Year, height=Tot.Num, fill=Tot.Num)) +
  geom_ridgeline_gradient(scale=0.0006) +
  scale_fill_distiller(palette="PuBu", breaks=seq(0,1600,200), direction=1, name="Total") +
  theme_classic()+
  geom_hline(yintercept=5.6, linetype="dashed", color="seagreen")+
  geom_hline(yintercept=9, colour="seagreen")+
  geom_hline(yintercept=11.5, linetype="dashed", colour="seagreen")

ridge.norec.F <- NoRecruits %>% 
  filter(Status=="Fished") %>% 
  ggplot(., aes(x=Age, y=Year, height=Number, fill=Number)) +
  geom_ridgeline_gradient(scale=0.0009) +
  scale_fill_distiller(palette="PuBu", breaks=seq(0,1000,100), direction=1) +
  theme_classic()+
  geom_hline(yintercept=5.6, linetype="dashed", color="seagreen")+
  geom_hline(yintercept=9, colour="seagreen")+
  geom_hline(yintercept=11.5, linetype="dashed", colour="seagreen")
  
ridge.norec.NTZ <- NoRecruits %>% 
  filter(Status=="NTZ") %>% 
  ggplot(., aes(x=Age, y=Year, height=Number, fill=Number)) +
  geom_ridgeline_gradient(scale=0.00125) +
  scale_fill_distiller(palette="PuBu", breaks=seq(0,800,100), direction=1) +
  theme_classic()+
  geom_hline(yintercept=5.6, linetype="dashed", color="seagreen")+
  geom_hline(yintercept=9, colour="seagreen")+
  geom_hline(yintercept=11.5, linetype="dashed", colour="seagreen") 

#### LINE PLOTS BY AGE GROUP ####

line.recruits <- WholePop %>% 
  filter(Age==1) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                  ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) %>% 
  ggplot(., aes(x=Year, y=Number, group=Status, colour=NTGroup, shape=Status))+
  geom_point(size=2.5)+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="seagreen")+
  geom_vline(xintercept=9, colour="seagreen")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="seagreen") 

line.sublegal <- WholePop %>% 
  filter(Age>1 & Age <=4) %>% 
  group_by(Year, Status) %>% 
  mutate(Total = sum(Number)) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) %>% 
  ggplot(., aes(x=Year, y=Total, group=interaction(Status,Age), colour=NTGroup, shape=Status))+
  geom_point(size=2.5)+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="seagreen")+
  geom_vline(xintercept=9, colour="seagreen")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="seagreen") 

line.legal <- WholePop %>% 
  filter(Age>4 & Age <=10) %>% 
  group_by(Year, Status) %>% 
  mutate(Total = sum(Number)) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) %>% 
  ggplot(., aes(x=Year, y=Total, group=Status, colour=NTGroup, shape=Status))+
  geom_point(size=2.5)+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="seagreen")+
  geom_vline(xintercept=9, colour="seagreen")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="seagreen") 

line.biglegal <- WholePop %>% 
  filter(Age>10) %>% 
  group_by(Year, Status) %>% 
  mutate(Total = sum(Number)) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) %>% 
  ggplot(., aes(x=Year, y=Total, group=interaction(Status,Age), colour=NTGroup, shape=Status))+
  geom_point(size=2.5)+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="seagreen")+
  geom_vline(xintercept=9, colour="seagreen")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="seagreen") 


#### Kobe Style Plot ####
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

#### PLOT SIMULATION ####

sYear5 <- readRDS("sim01_YearlyTotal.5") 
sYear5.NT <- sYear5[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
sYear5.F <- sYear5[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

sYear10 <- readRDS("sim01_YearlyTotal.10")
sYear10.NT <- sYear10[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
sYear10.F <- sYear10[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

sYear15 <- readRDS("sim01_YearlyTotal.15") 
sYear15.NT <- sYear15[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
sYear15.F <- sYear15[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

sYear20 <- readRDS("sim01_YearlyTotal.20")
sYear20.NT <- sYear20[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
sYear20.F <- sYear20[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

sYear25 <- readRDS("sim01_YearlyTotal.25")
sYear25.NT <- sYear25[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
sYear25.F <- sYear25[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

# First sanctuaries at year 27
sYear30 <- readRDS("sim01_YearlyTotal.30")
sYear30.NT <- sYear30[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
sYear30.F <- sYear30[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

sYear35 <- readRDS("sim01_YearlyTotal.35")
sYear35.NT <- sYear35[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
sYear35.F <- sYear35[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

sYear40 <- readRDS("sim01_YearlyTotal.40")
sYear40.NT <- sYear40[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
sYear40.F <- sYear40[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

sYear45 <- readRDS("sim01_YearlyTotal.45")
sYear45.NT <- sYear45[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
sYear45.F <- sYear45[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

## Sanctuaries expanded
sYear50 <- readRDS("sim01_YearlyTotal.50")
sYear50.NT <- sYear50[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
sYear50.F <- sYear50[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

## Commonwealth sanctuary put in at year 53
sYear55 <- readRDS("sim01_YearlyTotal.55")
sYear55.NT <- sYear55[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
sYear55.F <- sYear55[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

sYear59 <- readRDS("sim01_YearlyTotal.59")
sYear59.NT <- sYear59[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
sYear59.F <- sYear59[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)


#### CREATE POP TOTAL FOR SIMULATION ####

sNoTakeAges <- as.data.frame(cbind(sYear5.NT, sYear10.NT, sYear15.NT, sYear20.NT, sYear25.NT, sYear30.NT, sYear35.NT, sYear40.NT, sYear45.NT, sYear50.NT, sYear55.NT, sYear59.NT)) %>% 
  mutate(Status = "NTZ") %>% 
  mutate(Age = seq(1:30)) %>% 
  mutate(Scenario = "Without NTZ") %>% 
  mutate(Scenario = as.factor(Scenario))
colnames(sNoTakeAges) <- c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018", "Status", "Age", "Scenario")


sFishedAges <- as.data.frame(cbind(sYear5.F, sYear10.F, sYear15.F, sYear20.F, sYear25.F, sYear30.F, sYear35.F, sYear40.F, sYear45.F, sYear50.F, sYear55.F, sYear59.F)) %>% 
  mutate(Status = "Fished") %>% 
  mutate(Age = seq(1:30)) %>% 
  mutate(Scenario = "Without NTZ") %>% 
  mutate(Scenario = as.factor(Scenario))
colnames(sFishedAges) <- c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018", "Status", "Age", "Scenario")

sWholePop <- rbind(sNoTakeAges, sFishedAges) %>%  
  pivot_longer(cols=-c(Age, Status,Scenario), names_to="Year", values_to="Number") %>% 
  mutate(Year = as.factor(Year))

sWholePop <- sWholePop %>% 
  mutate(Year = fct_relevel(Year, c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2018"))) %>% 
  mutate(Status = as.factor(Status)) 

sNoRecruits <- sWholePop %>% 
  filter(Age!=1) #%>% 

#### RIDGELINE PLOTS WITHOUT RECRUITS ####

sridge.norec.whole <- sNoRecruits %>% 
  group_by(Age, Year) %>% 
  mutate(Tot.Num = sum(Number)) %>% 
  ungroup() %>% 
  ggplot(., aes(x=Age, y=Year, height=Tot.Num, fill=Tot.Num)) +
  geom_ridgeline_gradient(scale=0.0006) +
  scale_fill_distiller(palette="PuBu", breaks=seq(0,1600,200), direction=1, name="Total") +
  theme_classic()+
  geom_hline(yintercept=5.6, linetype="dashed", color="seagreen")+
  geom_hline(yintercept=9, colour="seagreen")+
  geom_hline(yintercept=11.5, linetype="dashed", colour="seagreen")

sridge.norec.F <- sNoRecruits %>% 
  filter(Status=="Fished") %>% 
  ggplot(., aes(x=Age, y=Year, height=Number, fill=Number)) +
  geom_ridgeline_gradient(scale=0.0009) +
  scale_fill_distiller(palette="PuBu", breaks=seq(0,1000,100), direction=1) +
  theme_classic()+
  geom_hline(yintercept=5.6, linetype="dashed", color="seagreen")+
  geom_hline(yintercept=9, colour="seagreen")+
  geom_hline(yintercept=11.5, linetype="dashed", colour="seagreen")

sridge.norec.NTZ <- sNoRecruits %>% 
  filter(Status=="NTZ") %>% 
  ggplot(., aes(x=Age, y=Year, height=Number, fill=Number)) +
  geom_ridgeline_gradient(scale=0.00125) +
  scale_fill_distiller(palette="PuBu", breaks=seq(0,800,100), direction=1) +
  theme_classic()+
  geom_hline(yintercept=5.6, linetype="dashed", color="seagreen")+
  geom_hline(yintercept=9, colour="seagreen")+
  geom_hline(yintercept=11.5, linetype="dashed", colour="seagreen") 

#### LINE PLOTS BY AGE GROUP ####

sline.recruits <- sWholePop %>% 
  filter(Age==1) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) %>% 
  ggplot(., aes(x=Year, y=Number, group=Status, colour=NTGroup, shape=Status))+
  geom_point(size=2.5)+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="seagreen")+
  geom_vline(xintercept=9, colour="seagreen")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="seagreen") 

sline.sublegal <- sWholePop %>% 
  filter(Age>1 & Age <=4) %>% 
  group_by(Year, Status) %>% 
  mutate(Total = sum(Number)) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) %>% 
  ggplot(., aes(x=Year, y=Total, group=interaction(Status,Age), colour=NTGroup, shape=Status))+
  geom_point(size=2.5)+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="seagreen")+
  geom_vline(xintercept=9, colour="seagreen")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="seagreen") 

sline.legal <- sWholePop %>% 
  filter(Age>4 & Age <=10) %>% 
  group_by(Year, Status) %>% 
  mutate(Total = sum(Number)) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) %>% 
  ggplot(., aes(x=Year, y=Total, group=Status, colour=NTGroup, shape=Status))+
  geom_point(size=2.5)+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="seagreen")+
  geom_vline(xintercept=9, colour="seagreen")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="seagreen") 

sline.biglegal <- sWholePop %>% 
  filter(Age>10) %>% 
  group_by(Year, Status) %>% 
  mutate(Total = sum(Number)) %>% 
  mutate(NTGroup = ifelse(Year %in% c("1960","1965","1970","1975","1980","1985"), 1, ifelse(Year %in% c("1990","1995","2000","2005"), 2, 
                                                                                            ifelse(Year %in% c("2010","2015"),3, ifelse(Year %in% c("2018"), 4, 0))))) %>% 
  mutate(NTGroup = as.factor(NTGroup)) %>% 
  ggplot(., aes(x=Year, y=Total, group=interaction(Status,Age), colour=NTGroup, shape=Status))+
  geom_point(size=2.5)+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="seagreen")+
  geom_vline(xintercept=9, colour="seagreen")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="seagreen") 


#### Kobe Style Plot ####
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

#### PUT PLOTS TOGETHER ####

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

## Put the plots together 
#Put together as one plot
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

## Line Plot of Whole Population
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

sTotalPop <- array(0, dim=c(59,2))
sTotalPop <- as.data.frame(sTotalPop)
sTotalPop <- sTotalPop %>% 
  rename(Year="V1") %>% 
  mutate(Year=seq(1960,2018,1)) %>% 
  rename(sTot.Pop="V2")

numYear <- seq(1,59,1)

for(Y in 1:59){
  
  year <- readRDS(paste0("sim01_YearlyTotal.", numYear[Y]))
  year <- rowSums(year[,,1:30], dim=2)
  year <- sum(year[,12])
  
  sTotalPop[Y,2] <- year
  
}

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


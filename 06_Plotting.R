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

Year5.F <- readRDS("YearlyTotal.1") %>% 
  colSums(., dim=2)
Year5.NT <- array(0, dim=c(30, 1))

Year10.F <- readRDS("YearlyTotal.2")%>% 
  colSums(., dim=2)
Year10.NT <- array(0, dim=c(30, 1))

Year15.F <- readRDS("YearlyTotal.3")%>% 
  colSums(., dim=2)
Year15.NT <- array(0, dim=c(30, 1))

Year20.F <- readRDS("YearlyTotal.4")%>% 
  colSums(., dim=2)
Year20.NT <- array(0, dim=c(30, 1))

Year25.F <- readRDS("YearlyTotal.5")%>% 
  colSums(., dim=2)
Year25.NT <- array(0, dim=c(30, 1))

# First sanctuaries at year 27
Year30 <- readRDS("YearlyTotal.6")
Year30.NT <- Year30[c(NoTake[[1]]),12, ] %>% 
  colSums(.)
Year30.F <- Year30[-c(NoTake[[1]]),12, ] %>% 
  colSums(.)
  
Year35 <- readRDS("YearlyTotal.7")
Year35.NT <- Year35[c(NoTake[[1]]),12, ] %>% 
  colSums(.)
Year35.F <- Year35[-c(NoTake[[1]]),12, ] %>% 
  colSums(.)

Year40 <- readRDS("YearlyTotal.8")
Year40.NT <- Year40[c(NoTake[[1]]),12, ] %>% 
  colSums(.)
Year40.F <- Year40[-c(NoTake[[1]]),12, ] %>% 
  colSums(.)

Year45 <- readRDS("YearlyTotal.9")
Year45.NT <- Year45[c(NoTake[[1]]),12, ] %>% 
  colSums(.)
Year45.F <- Year45[-c(NoTake[[1]]),12, ] %>% 
  colSums(.)

## Sanctuaries expanded
Year50 <- readRDS("YearlyTotal.10")
Year50.NT <- Year50[c(NoTake[[2]]),12, ] %>% 
  colSums(.)
Year50.F <- Year50[-c(NoTake[[2]]),12, ] %>% 
  colSums(.)

## Commonwealth sanctuary put in at year 53
Year55 <- readRDS("YearlyTotal.11")
Year55.NT <- Year55[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year55.F <- Year55[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)

Year59 <- readRDS("YearlyTotal.12")
Year59.NT <- Year59[c(NoTake[[3]]),12, ] %>% 
  colSums(.)
Year59.F <- Year59[-c(NoTake[[3]]),12, ] %>% 
  colSums(.)


#### CREATE POP TOTAL ####

NoTakeAges <- as.data.frame(cbind(Year5.NT, Year10.NT, Year15.NT, Year20.NT, Year25.NT, Year30.NT, Year35.NT, Year40.NT, Year45.NT, Year50.NT, Year55.NT, Year59.NT)) %>% 
  mutate(Status = "NTZ") %>% 
  mutate(Age = seq(1:30)) %>% 
colnames(NoTakeAges) <- c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2019", "Status", "Age")


FishedAges <- as.data.frame(cbind(Year5.F, Year10.F, Year15.F, Year20.F, Year25.F, Year30.F, Year35.F, Year40.F, Year45.F, Year50.F, Year55.F, Year59.F)) %>% 
  mutate(Status = "Fished") %>% 
  mutate(Age = seq(1:30))
colnames(FishedAges) <- c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2019", "Status", "Age")

WholePop <- rbind(NoTakeAges, FishedAges) %>%  
  pivot_longer(cols=-c(Age, Status), names_to="Year", values_to="Number") %>% 
  mutate(Year = as.factor(Year))

WholePop <- WholePop %>% 
  mutate(Year = fct_relevel(Year, c("1965", "1970", "1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2019"))) %>% 
  mutate(Status = as.factor(Status))
  
NoRecruits <- WholePop %>% 
  filter(Age!=1) %>% 
  mutate(Number = Number/1000) # Ridgeline doesn't like big numbers so we just make all the numbers smaller 

#### RIDGELINE PLOTS WITHOUT RECRUITS ####

ridge.norec.F <- NoRecruits %>% 
  filter(Status=="Fished") %>% 
  ggplot(., aes(x=Age, y=Year, height=Number)) +
  geom_ridgeline(scale=0.045) +
  theme_classic()+
  geom_hline(yintercept=5.6, linetype="dashed", color="seagreen")+
  geom_hline(yintercept=9, colour="seagreen")+
  geom_hline(yintercept=11.5, linetype="dashed", colour="seagreen")
  
ridge.norec.NTZ <- NoRecruits %>% 
  filter(Status=="NTZ") %>% 
  ggplot(., aes(x=Age, y=Year, height=Number)) +
  geom_ridgeline(scale=0.45) +
  theme_classic()+
  geom_hline(yintercept=5.6, linetype="dashed", color="seagreen")+
  geom_hline(yintercept=9, colour="seagreen")+
  geom_hline(yintercept=11.5, linetype="dashed", colour="seagreen") 

#### LINE PLOTS BY AGE GROUP ####

line.recruits <- WholePop %>% 
  filter(Age==1) %>% 
  ggplot(., aes(x=Year, y=Number, group=Status, colour=Status))+
  geom_point()+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="seagreen")+
  geom_vline(xintercept=9, colour="seagreen")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="seagreen") 

line.sublegal <- WholePop %>% 
  filter(Age>1 & Age <=4) %>% 
  group_by(Year, Status) %>% 
  mutate(Total = sum(Number)) %>% 
  ggplot(., aes(x=Year, y=Total, group=interaction(Status,Age), colour=Status))+
  geom_point()+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="seagreen")+
  geom_vline(xintercept=9, colour="seagreen")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="seagreen") 

line.legal <- WholePop %>% 
  filter(Age>4 & Age <=10) %>% 
  group_by(Year, Status) %>% 
  mutate(Total = sum(Number)) %>% 
  ggplot(., aes(x=Year, y=Total, group=Status, colour=Status))+
  geom_point()+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="seagreen")+
  geom_vline(xintercept=9, colour="seagreen")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="seagreen") 

line.biglegal <- WholePop %>% 
  filter(Age>10) %>% 
  group_by(Year, Status) %>% 
  mutate(Total = sum(Number)) %>% 
  ggplot(., aes(x=Year, y=Total, group=interaction(Status,Age), colour=Status))+
  geom_point()+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=5.6, linetype="dashed", color="seagreen")+
  geom_vline(xintercept=9, colour="seagreen")+
  geom_vline(xintercept=11.5, linetype="dashed", colour="seagreen") 















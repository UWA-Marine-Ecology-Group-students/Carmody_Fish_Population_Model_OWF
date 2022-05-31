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

Year5 <- readRDS("YearlyTotal.1") %>% 
  colSums(., dim=2)
Year10 <- readRDS("YearlyTotal.2")%>% 
  colSums(., dim=2)
Year15 <- readRDS("YearlyTotal.3")%>% 
  colSums(., dim=2)
Year20 <- readRDS("YearlyTotal.4")%>% 
  colSums(., dim=2)
Year25 <- readRDS("YearlyTotal.5")%>% 
  colSums(., dim=2)
Year30 <- readRDS("YearlyTotal.6")%>% 
  colSums(., dim=2)
Year35 <- readRDS("YearlyTotal.7")%>% 
  colSums(., dim=2)
Year40 <- readRDS("YearlyTotal.8")%>% 
  colSums(., dim=2)
Year45 <- readRDS("YearlyTotal.9")%>% 
  colSums(., dim=2)
Year50 <- readRDS("YearlyTotal.10")%>% 
  colSums(., dim=2)
Year55 <- readRDS("YearlyTotal.11")%>% 
  colSums(., dim=2)
Year59 <- readRDS("YearlyTotal.12")%>% 
  colSums(., dim=2)

#### CREATE POP TOTAL ####
WholePop <- as.data.frame(cbind(Year5, Year10, Year15, Year20, Year25, Year30, Year35, Year40, Year45, Year50, Year55, Year59)) %>% 
  mutate(Age = seq(1,30,)) %>% 
  pivot_longer(!Age, names_to="Year", values_to="Number") %>% 
  mutate(Year = as.factor(Year))

WholePop <- WholePop %>% 
  mutate(Year=fct_recode(Year, "1965"="Year5", "1970"="Year10", "1975"="Year15", "1980"="Year20", "1985"="Year25", "1990"="Year30", "1995"="Year35", "2000"="Year40", "2005"="Year45", "2010"="Year50", "2015"="Year55", "2019"="Year59"))
  
NoRecruits <- WholePop %>% 
  filter(Age!=1) 

# Need to standardise the values so they are between 1 and 100
NoRecruits.S <- NoRecruits %>% 
  group_by(Year) %>% 
  mutate(Standardised = (Number-min(Number))/(max(Number)-min(Number))) %>% 
  ungroup()

#### RIDGELINE WHOLE POPULATION ####

ridge.whole.pop <- ggplot(NoRecruits.S, aes(x=Age, y=Year, height=Standardised)) +
  geom_ridgeline() +
  theme_classic()
  
  


















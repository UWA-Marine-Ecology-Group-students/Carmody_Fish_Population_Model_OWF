## ---------------------------
##
## Script name: MDS Plot for NSW Fish Species 
##
## Purpose of script: To determine groups of co-occurring species which may influence site choice in recreational fishers 
##
## Author: Charlotte Aston
##
## ----------------------------
## Notes:
##  -Catch data is taken from 1 survey: NSW state-wide phone diary survey
## ---------------------------

#### Load libraries ####
library(tidyverse)
library(dplyr)
library(MASS)
library(vegan)

#### Set working directories ####
working.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

## Save these directory names to use later----
data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")

#### Read in Data ####
setwd(data_dir)
codes <- read.csv("NRFS_SpeciesCodes For Matt.csv")
catch_data <- read.csv("NSWACT_2013-14_Data for RUM_2.csv")

#### Data cleaning ####

catch <- catch_data%>%
  dplyr::select("SpeciesID", "NKept", "NReleased", "DiaryEventPerson_ComboID")%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  mutate(Total = NKept+NReleased)

codes <- codes%>%
  rename("SpeciesID"="NRFS.CODE")

NSWFishData <- merge(catch, codes)

### Create a MDS plot ###

MDS1Data <- NSWFishData%>%
  dplyr::select("DiaryEventPerson_ComboID", "Total", "Name")%>%
  pivot_wider(names_from=Name, values_from=Total)%>%
  dplyr::mutate(ID = row_number())%>%
  dplyr::select(-DiaryEventPerson_ComboID)%>%
  mutate_all(~replace(., is.na(.), 0))

MDS1Data <- MDS1Data[1:1000,]

MDS1<- metaMDS(comm = MDS1Data, distance = "bray", trace = FALSE, autotransform = FALSE)

plot(MDS1$points)

MDS1$stress












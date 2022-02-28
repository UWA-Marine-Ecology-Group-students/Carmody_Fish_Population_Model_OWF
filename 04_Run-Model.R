###################################################

# Script just for running the full model 
# Requires files that are made in the first script
# Will run the full model with fishing
# Make sure you reset things like NCELL if you've
# used the test script

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


#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
m_dir <- paste(working.dir, "Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")

#### PRE-SETS ####

## Create colours for the plot
cols <- brewer.pal(8, "RdBu")
levels_water <- data.frame(c("<1000", "1000-1100", "1100-1200", "1200-1300",
                             "1300-1400", "1400-1500", "1500-1600", ">1600"))
names(levels_water)[1] <- "Levels"
levels_water$Levels <- as.factor(levels_water$Levels)
names(cols) <- levels(levels_water$Levels)

## Read in functions
setwd(working.dir)
source("X_Functions.R")

#### LOAD FILES ####
setwd(sg_dir)
movement <- readRDS("movement")
juv_movement <- readRDS("juvmove")
fishing <- readRDS("Fishing")
NoTake <- readRDS("NoTake")
water <- readRDS("water")
  
#### PARAMETER VALUES ####

## Natural Mortality
# We have instantaneous mortality from Marriot et al 2011 and we need to convert that into monthly mortality
yearly_surv=exp(-0.146)
monthly_mort=1-(yearly_surv^(1/12))

M <- monthly_mort # Natural mortality rate per month

#Ricker recruitment model parameters (these are currently just made up values)
a <- 7
b <- 0.0017
M50 <- 2 # From Grandcourt et al. 2010
M95 <- 5 # From Grandcourt et al. 2010 (technically M100)

#Fish movement parameters
SwimSpeed <- 1.0 # Swim 1km in a day - this is completely made up 

## Fishing mortality parameters
A50 <- 4 # For L. miniatus from Williams et al. 2010 # L. miniatus becomes vulnerable to fishing at about age two
A95 <- 6 # For L. miniatus from Williams et al. 2012
q <- 0.5 # Apparently this is what lots of stock assessments set q to be...

#### RUN MODEL ####

# Fill the array with random numbers of fish for now
# Need to change the loop so that rather than filling in a column it adds a new column each time 

YearlyTotal <- array(0, dim = c(NCELL,12,100)) #This is our yearly population split by age category (every layer is an age group)
for(d in 1:100){
  YearlyTotal[,1,] <- matrix(floor(runif(NCELL, 1, 2000)))
}

PopTotal <- array(0, dim=c(NCELL, 12, 59)) # This is our total population, all ages are summed and each column is a month (each layer is a year)

Movement <- NULL

NCELL <- nrow(water)

# The time steps now represent 1 month (need to make sure you change the distance your fish can swim to match the timestep)

for(YEAR in 2:59){
  
  for(MONTH in 2:13){
    
      for(A in 2:dim(YearlyTotal)[3]){
        
        movement.func(Age=A, Month=MONTH, Population=YearlyTotal, Max.Cell=NCELL, Adult.Move= movement, Juv.Move=juv_movement)

      } #End bracket for movement 
    
      ## Fishing Mortality
      for(A in 1:dim(YearlyTotal)[3]){
        
        mortality.func(Age=A, mort.50=A50, mort.95=A95, Nat.Mort=M, NTZ=NoTake, Effort=Fishing, Cell=CELL, Max.Cell = NCELL,
                       Month=MONTH, Year=YEAR, Population=YearlyTotal)
      } # End Mortality
    
    # ## Recruitment
    # Recs <- matrix(0, nrow=1, ncol=1) #Blank matrix to add the recruits to
    # Recruitment <- as.matrix(colSums(Pop[,t,], dims=1)) 
    # for(A in 1:dim(Recruitment)[1]){
    #   Mature <- 1/(1+(exp(-log(19)*((A-M95)/(M95-M50))))) #Number of mature individuals in each age class
    #   S <- colSums(Recruitment)*Mature #Spawning stock
    #   Rec <- a*S*exp(-b*S) #Number of recruits from that age class
    #   Recs <- rbind(Recs, Rec) #Combine recruits from each age class into one dataframe
    # }
    
    PopTotal[ , , YEAR] <- rowSums(YearlyTotal) 
    # Collapses the 3D matrix and gives you the number of fish of all ages in each cell in each month
    
  }
  print(YEAR)
  water$pop <- rowSums(PopTotal[ ,12,YEAR]) # We just want the population at the end of the year
  
  water <- water%>%
    mutate(Population = ifelse(pop < 1000, "<1000",
                               ifelse (pop>1000 & pop<110, "1000-1100",
                                       ifelse (pop>1100 & pop<1200, "1100-1200",
                                               ifelse (pop>1200 & pop<1300, "1200-1300",
                                                       ifelse(pop>1300 & pop<1400, "1300-1400",
                                                              ifelse(pop>1400 & pop<1500, "1400-1500",
                                                                     ifelse(pop>1500 & pop<1600, "1500-1600", ">1600"
                                                                            
                                                                     ))))))))%>%
    mutate(Population = factor(Population))
  
  water$Population <- fct_relevel(water$Population, "<1000",  "1000-1100", "1100-1200", "1200-1300",
                                  "1300-1400", "1400-1500", "1500-1600", ">1600")
  
  print(map <- ggplot(water)+
          geom_sf(aes(fill=Population, color=Fished))+
          scale_fill_manual(name="Population", values= cols, drop=FALSE)+
          scale_color_manual(name="Fished", values=c("white", "black"))+
          theme_void())
  Sys.sleep(3)
}


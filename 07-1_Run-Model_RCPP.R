###################################################

# Running the full model with whole area and fishing
# Produces plots every five years
# Get spatial plots as well as age based histograms
# See script six for more complex plots
###################################################

library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(forcats)
library(RColorBrewer)
library(MQMF)

# rm(list = ls())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

# data_dir <- paste(working.dir, "Data", sep="/")
# fig_dir <- paste(working.dir, "Figures", sep="/")
# sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")
# pop_dir <-  paste(working.dir, "Output_Population", sep="/")
# sim_dir <- paste(working.dir, "Simulations", sep="/")

## Read in functions
setwd(working.dir)
library("Rcpp")
sourceCpp("Model_RccpArm.cpp")

#### PRE-SETS ####

## Create colours for the plot
pop.groups <- c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150)
my.colours <- "PuBu"

model.name <- "ningaloo"

#### LOAD FILES ####

## Normal Model Files
setwd(sg_dir)
movement <- readRDS(paste0(model.name, sep="_", "movement"))
# juv_movement <- readRDS(paste0(model.name, sep="_","juvmove")) 
recruitment <- readRDS(paste0(model.name, sep="_","recruitment"))
fishing <- readRDS(paste0(model.name, sep="_", "fishing"))
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))
water <- readRDS(paste0(model.name, sep="_","water"))
YearlyTotal <- readRDS(paste0(model.name, sep="_", "BurnInPop"))
selectivity <- readRDS("selret")
maturity <- readRDS("maturity")
weight <- readRDS("weight")

## Simulation Files
# setwd(sim_dir)
# fishing <- readRDS(paste0(model.name, sep="_", "S02_fishing"))

#### PARAMETER VALUES ####
## Natural Mortality
# We have instantaneous mortality from Marriott et al 2011 and we need to convert that into monthly mortality
M <- 0.146
step <- 1/12 # We're doing a monthly time step here

# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
alpha <- 0.4344209 #0.4344209
beta <-	0.01889882

NCELL <- nrow(water)
Ages <- seq(1,30) #These are the ages you want to plot 
Time <- seq(1,59) #This is how long you want the model to run for
# PlotTotal <- T #This is whether you want a line plot of the total or the map

Pop.Groups <- seq(1,12)

#### SET UP LISTS TO HOLD THE PLOTS ####
SpatialPlots <- list()
LengthPlots <- list()
AgePlots <- list()
TimesPlotted <- 0

#### SET UP INITIAL POPULATION ####
PopTotal <- array(0, dim=c(NCELL, 12, length(Time))) # This is our total population, all ages are summed and each column is a month (each layer is a year)
Total <- array(NA, dim=c(length(Time),1)) # For plotting
bio.catch <- array(0, dim=c(NCELL, length(Ages)))
monthly.catch <- array(0, dim=c(length(Ages), 12))
yearly.catch <- array(0, dim=(c(length(Time), 3)))

#### RUN MODEL ####
BurnIn = F #This is to swap the model between burn in and running the model properly
# setwd(sim_dir)





WL_a = 0.000028
WL_b = 2.8761
Linf = 664
vbK = 0.241
tzero = -0.375
Weight = as.matrix(data.frame(matrix(nrow=MaxAge,ncol=12)))
colnames(Weight)=1:12
DecAge=0
for (a in 1:MaxAge) {
  for (m in 1:12) {
    DecAge = DecAge + 1/12
    EstLen = Linf * (1 - exp(-vbK * (DecAge-tzero)))
    Weight[a,m] = WL_a * EstLen ^ WL_b
  }
}


YEAR=0
MONTH=0
Age=0
Max_Cell=NCELL
Adult_Move=movement
Population=YearlyTotal
MaxAge=30
Max_Year=2
BHa = as.numeric(alpha)
BHb= as.numeric(beta)
PF=0.5
Mature=as.vector(unlist(maturity))
settlement=as.vector(recruitment[,2])

Start=Sys.time()


ModelOutput <- RunModelfunc_cpp(MaxAge, MaxYear, Max_Cell, Nat_Mort, BHa, BHb, PF, Adult_Move, Mature, Weight, settlement, 
                                       Population, Select, Effort)


End=Sys.time()
Runtime = End - Start
Runtime


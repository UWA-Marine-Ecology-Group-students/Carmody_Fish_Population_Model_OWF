###################################################

# Setting up initial population to establish an 
# equilibrium level of mortality and recruitment

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
library(abind)


#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
m_dir <- paste(working.dir, "Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")

model.name <- "ningaloo"

#### PARAMETERS FOR THE FISH POPULATION #####

## Set timestep 
step = 1/12 # We're doing a monthly timestep here

## Von Bertalanffy Parameters
Linf = 664 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
k = 0.241 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
t0 = -0.375 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr

## Weight-Length Relationship
WLa = 0.000028 # This is for mm and g
WLb = 2.8761

## Proportion expected to be females
prop.fem = 0.5

## Natural Mortality
M = 0.146  #Marriot et al. 2011

# Beverton-Holt Parameters
h = 0.76 # Steepness 0.76
R0 = 1 # Initial recruitment

# Selectivity
A50 = 3#0.7 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
A95 = 4#1.4 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
# Next runs this probably needs to come up to five
MaxAge = 30
nYears= 59

PRM = 0.25 # Post release mortality for our retention function

# 270 in 1913 to 1961 from Fisheries Act Amendment
# 202mm in 1961 - From the Fisheries Act Amendment
# 280mm in 1991 https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1054&context=fr_fmp
# Was 410mm in 1995 https://www.mediastatements.wa.gov.au/Pages/Court/1995/04/Revised-bag-and-size-limits-for-recreational-fishing-announced.aspx

# Maturity
M50 = 3.62 # From Marriott et al. 2011
M95 = 5.97 # From Marriott et al. 2011

# Fishing parameters
eq.init.fish = 0.025

#### SET UP LIFE HISTORY VALUES FOR L. NEBULOSUS ####

Life.History <- as.data.frame(array(1, dim=c(361, 3))) %>% 
  rename(Age = "V1") %>% 
  rename(Length = "V2") %>% 
  rename(Weight = "V3")

for (r in 1:360){                                    # This sets the ages to be every month, but in decimals
  Life.History[r+1, 1] <- Life.History[r, 1] + step 
}

Life.History <- Life.History %>% 
  mutate(Length = Linf*(1-exp(-k*(Age-t0)))) %>% 
  mutate(Weight = (WLa*(Length^WLb)/1000)) # Divide by 1000 to get in kg 

#### SET UP SURVIVAL AND MATURITY FOR UNFISHED POPULATION ####

Unfished.Pop.SetUp <- as.data.frame(array(1, dim=c(361, 3))) %>% 
  rename(Unfish.Surv = "V1") %>% 
  rename(Unfish.Mat = "V2") %>% 
  rename(Unfish.Bio = "V3")

## Unfished Population Survival in Each Age Group
Unfished.Pop.SetUp[1, 1] <- R0*prop.fem

for (r in 2:361){                                    # This calculates the survival in the next age based on the previous age
  Unfished.Pop.SetUp[r, 1] <- Unfished.Pop.SetUp[r-1, 1]*exp(-M/12) 
}

## Unfished Population Proportion Mature per Recruit in Each Age Group

Unfished.Pop.SetUp <- Unfished.Pop.SetUp %>% 
  mutate(Age = Life.History$Age) %>% 
  mutate(Unfish.Mat = 1/(1+(exp(-log(19)*((Age-M50)/(M95-M50))))))

## Unfished Mature Biomass per Recruit in Each Age Group

Unfished.Pop.SetUp <- Unfished.Pop.SetUp %>% 
  mutate(Unfish.Bio = Unfish.Mat * Unfish.Surv * Life.History$Weight)

## Calculate the Total Mature Unfished Biomass per Recruit

Unfish.Fem.SSB <- Unfished.Pop.SetUp %>%  # This is the total mature female spawning biomass per recruit for the unfished population
  slice(which(row_number() %% 12 == 1)) %>% 
  summarise(., sum=sum(Unfish.Bio))

#### CALCULATE BEVERTON-HOLT PARAMETERS FOR UNFISHED POPULATION ####

alpha <- (Unfish.Fem.SSB/R0)*((1-h)/(4*h))

beta <- ((h-0.2)/(0.8*h*R0))

#### SET UP SURVIVAL AND MATURITY FOR FISHED POPULATION ####

Fished.Pop.SetUp <- as.data.frame(array(1, dim=c(361, 5))) %>% 
  rename(Selectivity = "V1") %>% 
  rename(Fishing.Mort = "V2") %>% 
  rename(TotMort = "V3") %>% 
  rename(Fish.Surv = "V4") %>% 
  rename(Fish.Mat = "V5")

## Selectivity for each age group - This is for the equilibrium population

Fished.Pop.SetUp <- Fished.Pop.SetUp %>% 
  mutate(Age = Life.History$Age) %>% 
  mutate(Length = Life.History$Length) %>% 
  mutate(Selectivity = 1/(1+(exp(-log(19)*((Age-A50)/(A95-A50))))))

## Next calculate fishing mortality 

Fished.Pop.SetUp <- Fished.Pop.SetUp %>% 
  mutate(Fishing.Mort = Selectivity * eq.init.fish)

## Caluclate total mortality

Fished.Pop.SetUp <- Fished.Pop.SetUp %>% 
  mutate(TotMort = Fishing.Mort + M)

## Calculate Fished Survival
Fished.Pop.SetUp[1, 4] <- R0*prop.fem

for (r in 2:361){                                    # This calculates the survival in the next age based on the previous age
  Fished.Pop.SetUp[r, 4] <- Fished.Pop.SetUp[r-1, 4]*exp(-Fished.Pop.SetUp[r-1,3]*step) # Divide by the time step here
}

## Calculate Mature Fished Biomass per Recruit

Fished.Pop.SetUp <- Fished.Pop.SetUp %>% 
  mutate(Fish.Mat = Unfished.Pop.SetUp$Unfish.Mat) %>% 
  mutate(Fish.Bio = Fish.Mat * Fish.Surv * Life.History$Weight)

## Calculate the Total Mature Fished Biomass per Recruit

Fish.Fem.SSB <- Fished.Pop.SetUp %>%  # This is the total mature female spawning biomass per recruit for the unfished population
  slice(which(row_number() %% 12 == 1)) %>% 
  summarise(., sum=sum(Fish.Bio))

##### CALCULATING SPAWNER PER RECRUIT AND EQUILIRBIRUM RECRUITMENT ####

SPR <- Fish.Fem.SSB/Unfish.Fem.SSB
  
Equil.Rec <- (Fish.Fem.SSB-alpha)/(beta*Fish.Fem.SSB)

#### CREATING AN INITIAL EQUILIBRIUM FISHED POPULATION FOR THE MODEL ####

## Our initial population uses all the numbers from the previous step of this code 
## The maturity at age, weight at age, fished per recruit survival
## the proportion of females at age
## We also need to decide on an initial level of recruitment (in thousands)

init.recruit <- 4000 # in thousands - normally 5000 for big model, 5 for small model

## Calculate initial fished recruitment

Init.Fish.Rec <- (Fish.Fem.SSB-alpha)/(Fish.Fem.SSB*beta)*init.recruit 

## Calculate initial unfished spawning biomass

Init.Unfish.Fem.SpawnBio <- Init.Fish.Rec * Unfish.Fem.SSB

## Calculate initial fished spawning biomass

Init.Fem.SpawnBio <- Init.Fish.Rec * Fish.Fem.SSB

## Calculate new Beverton-Holt Parameters

alpha <- (Init.Unfish.Fem.SpawnBio/Init.Fish.Rec)*((1-h)/(4*h))

beta <- ((h-0.2)/(0.8*h*Init.Fish.Rec))

## Calculate number of female recruits in the next time step

N.Female.Rec <- (Init.Fem.SpawnBio/(alpha+beta*Init.Fem.SpawnBio))*prop.fem

## Calculate survival into the next time step to create a full age structured population

Starting.Pop.For.Model <- as.data.frame(array(1, dim=c(361, 1)))

Starting.Pop.For.Model <- Starting.Pop.For.Model %>% 
  mutate(Age = Life.History$Age) 

Starting.Pop.For.Model[1,1] <- N.Female.Rec

for (r in 2:361){                                    # This calculates the survival in the next age based on the previous age using total mortality from our fished population
  Starting.Pop.For.Model[r, 1] <- Starting.Pop.For.Model[r-1, 1]*exp(-Fished.Pop.SetUp[r-1,3]*step) # Divide by the time step here
}

Starting.Pop.For.Model <- Starting.Pop.For.Model %>% 
  rename(N = "V1")

## These fish form our starting population for the model
## At the end of the year the spawning stock biomass of all females will be calculated to generate recruitment for the next year
## Alpha and Beta need to be recorded for use in the next step of the model

#### SELECTIVITY-RETENTION FOR THE FISHING IN THE MODEL ####
## Retention for each age group
# This is for when you actually run the model you don't use this for setting up the initial population
# Trying to account for the fact that fish that are below the legal size limit are likely to be thrown back and so won't necessarily die
Fished.Pop.SetUp <- Fished.Pop.SetUp %>% 
  mutate(Retention6091 = ifelse(Length<=200, 0, ifelse(Length > 200 & Length < 202, 0.5, ifelse(Length>=202, 0.95, 0)))) %>% 
  mutate(Retention9195 = ifelse(Length<=275, 0, ifelse(Length > 275 & Length < 280, 0.5, ifelse(Length>=280, 0.95, 0)))) %>% 
  mutate(Retention95 = ifelse(Length<=405, 0, ifelse(Length > 405 & Length < 410, 0.5, ifelse(Length>=410, 0.95, 0))))


## Landings and discards
# This gives us the proportion of fish that are kept and the proportion that are thrown back
## Selectivity-Retention Values Including post-release mortality

Fished.Pop.SetUp <- Fished.Pop.SetUp %>% 
  mutate(Landings6091 = Selectivity*Retention6091) %>% 
  mutate(Discards6091 = Selectivity*(1-Landings6091)) %>% 
  mutate(Landings9195 = Selectivity*Retention9195) %>% 
  mutate(Discards9195 = Selectivity*(1-Landings9195)) %>% 
  mutate(Landings95 = Selectivity*Retention95) %>% 
  mutate(Discards95 = Selectivity*(1-Landings95))

## Selectivity-Retention Values Including post-release mortality
Fished.Pop.SetUp <- Fished.Pop.SetUp %>% 
  mutate(SelRet6091 = Landings6091+(PRM*Discards6091)) %>% 
  mutate(SelRet9195 = Landings9195+(PRM*Discards9195)) %>% 
  mutate(SelRet95 = Landings95+(PRM*Discards95))


SelRet6091 <- Fished.Pop.SetUp$SelRet6091
SelRet6091 <- array(SelRet6091, dim=c(12,30))
SelRet6091 <- t(SelRet6091)

SelRet9195 <- Fished.Pop.SetUp$SelRet9195
SelRet9195 <- array(SelRet9195, dim=c(12,30))
SelRet9195 <- t(SelRet9195)

SelRet95 <- Fished.Pop.SetUp$SelRet95
SelRet95 <- array(SelRet95, dim=c(12,30))
SelRet95 <- t(SelRet95)

SelRet <- NULL
Change1 <- 31
Change2 <- 35

for(i in 1:Change1){
  SelRet <- abind(SelRet, SelRet6091, along=3)
}
for(i in 1:4){
  SelRet <- abind(SelRet, SelRet9195, along=3)
}
for(i in 1:24){
  SelRet <- abind(SelRet, SelRet95, along=3)
}

#### SAVING FILES ####

## We need our initial population 
setwd(sg_dir)
saveRDS(Starting.Pop.For.Model, file=paste0(model.name, sep="_", "Starting_Pop"))

## Also want to save our selectivity so we don't have to keep calculating it in the model

Selectivity <- Fished.Pop.SetUp$Selectivity
saveRDS(Selectivity, file="selectivity")
saveRDS(SelRet, file="selret")

## Similarly we don't want to have to keep calculating maturity

Maturity <- Fished.Pop.SetUp$Fish.Mat 
Maturity <- array(Maturity, dim=c(12,30))
Maturity <- t(Maturity)
saveRDS(Maturity, file="maturity")

## We need the weights for each age group as well

Weight <- Life.History$Weight
Weight <- array(Weight, dim=c(12,30))
Weight <- t(Weight)
saveRDS(Weight, file="weight")









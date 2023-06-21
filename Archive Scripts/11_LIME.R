###################################################

# Script for running LIME fishing mortality 
# analysis on the data from my model to get a value 
# of F for the sanctuary zones to make Kobe plots

###################################################
library(tidyverse)
library(dplyr)
library(ggplot2)
library(sf)
library(forcats)
library(RColorBrewer)
library(MQMF)
library(Rcpp)
library(RcppArmadillo)
library(raster)
library(sfnetworks)
library(abind)
library(matrixStats)
library(LIME)

rm(list = ls())
#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")
pop_dir <-  paste(working.dir, "Output_Population", sep="/")
sim_dir <- paste(working.dir, "Simulations", sep="/")

model.name <- "ningaloo"

## Read in functions
setwd(working.dir)
sourceCpp("X_Model_RccpArm.cpp")
source("X_Functions.R")

### Parameters
## Von Bertalanffy Parameters
Linf = 664 #(66.4 cm) # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
k = 0.241 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
t0 = -0.375 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr

## Weight-Length Relationship
WLa = 0.000028 # This is not in cm!
WLb = 2.8761

## Natural Mortality
M = 0.146  #Marriot et al. 2011
step <- 1/12

# Beverton-Holt Parameters
h = 0.76 # Steepness 0.76
R0 = 1 # Initial recruitment

# Selectivity
A50 = 0.7 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
A95 = 1.4 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr

# Maturity
M50 = 3.62 # From Marriott et al. 2011
M95 = 5.97 # From Marriott et al. 2011

##### READ IN AGE DATA AND FORMAT ####
setwd(sg_dir)
water <- readRDS(paste0(model.name, sep="_","water"))

setwd(pop_dir)

Age.NTZ <- list()

# Each layer is a simulation, rows are cells and columns are years
Age.NTZ[[1]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S00"))
Age.NTZ[[2]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S01"))
Age.NTZ[[3]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S02"))
Age.NTZ[[4]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ", sep="_", "S03"))

Names <- c("Historical and current management", "No NTZs or temporal management", 
           "NTZs and temporal management","Temporal management only" )

Ages.NTZ <- NULL
temp3 <- NULL

for(S in 1:4){
  
  temp <- Age.NTZ[[S]]
  
  for(SIM in 1:length(temp)){
    
    temp2 <- as.data.frame(colSums(temp[[SIM]])) %>% 
      mutate(Age = seq(1:30)) %>% 
      pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
      mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
    
    temp3 <- cbind(temp3, temp2$Number)
    
  }
  
  NTZ_Ages <- as.data.frame(temp3) %>% 
    mutate(Age = rep(1:30, each=59)) %>% 
    mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
    mutate(Mean_Pop = rowMeans(.[,1:(ncol(.)-2)])) %>%
    mutate(SD_Pop = rowSds(as.matrix(.[,1:(ncol(.)-2)]))) %>%
    mutate(Scenario = Names[S]) %>% 
    dplyr::select(Mean_Pop, SD_Pop, Scenario, Mod_Year, Age)
  
  Ages.NTZ <- rbind(Ages.NTZ, NTZ_Ages)
}

## Fished area 
Fished_ID <- water %>% 
  filter(Fished_2017=="Y") %>% 
  dplyr::select(ID)

Fished_ID <- Fished_ID$ID

setwd(pop_dir)
Age.F <- list()

# Each layer is a simulation, rows are cells and columns are years
Age.F[[1]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S00"))
Age.F[[2]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S01"))
Age.F[[3]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S02"))
Age.F[[4]] <- readRDS(paste0(model.name, sep="_", "Sp_Population_F", sep="_", "S03"))

Ages.F <- NULL
temp3 <- NULL

for(S in 1:4){
  
  temp <- Age.F[[S]]
  
  for(SIM in 1:length(temp)){
    
    temp2 <- as.data.frame(colSums(temp[[SIM]])) %>% 
      mutate(Age = seq(1:30)) %>% 
      pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
      mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
    
    temp3 <- cbind(temp3, temp2$Number)
    
  }
  
  F_Ages <- as.data.frame(temp3) %>% 
    mutate(Age = rep(1:30, each=59)) %>% 
    mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
    mutate(Mean_Pop = rowMeans(.[,1:(ncol(.)-2)])) %>%
    mutate(SD_Pop = rowSds(as.matrix(.[,1:(ncol(.)-2)]))) %>%
    mutate(Scenario = Names[S]) %>% 
    dplyr::select(Mean_Pop, SD_Pop, Scenario, Mod_Year, Age)
  
  Ages.F <- rbind(Ages.F, F_Ages)
}

## The whole model
setwd(pop_dir)
Age.All <- list()

# Each layer is a simulation, rows are cells and columns are years
Age.All[[1]] <- readRDS(paste0(model.name, sep="_", "Age_Distribution", sep="_", "S00"))
Age.All[[2]] <- readRDS(paste0(model.name, sep="_", "Age_Distribution", sep="_", "S01"))
Age.All[[3]] <- readRDS(paste0(model.name, sep="_", "Age_Distribution", sep="_", "S02"))
Age.All[[4]] <- readRDS(paste0(model.name, sep="_", "Age_Distribution", sep="_", "S03"))

Ages.All <- NULL
temp3 <- NULL

for(S in 1:4){
  
  # Ages
  temp <- Age.All[[S]]
  
  temp2 <- temp %>% 
    rowMeans(., dim=2) %>% 
    as.data.frame(.) %>% 
    mutate(Scenario = Names[S]) %>% 
    mutate(Age = seq(1,30,1))
  
  Ages.All <- rbind(Ages.All, temp2)
  
}

#### LIME ####
#* Specify biological inputs and starting values ####
lh <- create_lh_list(vbk = 0.241,
                     linf = 66.4,
                     t0 = -0.375,
                     lwa = 0.02105029,
                     lwb = 2.8761,
                     S50 = 3.6, #c(41.24376),
                     S95 = 4.5, #c(45.89186),
                     selex_input = "age",
                     selex_type = c("logistic"),
                     M50 = 3.62, # Is actually a little bit less than this but the maturity grid doesn't give 0.5 exactly, you get 0.51
                     M95 = 5.97,
                     h = 1, # matches what I have in my population set up
                     AgeMax = 30,
                     maturity_input = "age",
                     M = 0.146,
                     binwidth = 10,
                     #CVlen = 0.1, # Variability around the growth curve
                     SigmaR = 0.6, # SD around the recruitment variation
                     #SigmaF = 0.2, # SD around fishing mortality
                     #SigmaC = 0.2, # SD around the catch although I won't put catch data in
                     R0 = 1, # Initial unfished recruitment
                     #Frate = 0.2, # exponent determining the strength of coupling between effort and biomass
                     #Fequil = 0.5, # Prop.  SB to Unfished SB at bioeconomic equilib.
                     qcoef = 1e-06, # Catchability - changed it to be the same as the model
                     start_ages = 1,
                     rho = 0,
                     #theta = 1,
                     nseasons = 1, # Only change if you have shorter than annual data e.g. specues that only live 6 months
                     nfleets = 1 , # Will consider all recreational fishers to be one fleet
                     )

## Plot the inputs to make sure they look right
par(mfrow=c(2,2), mar=c(4,4,3,1))
plot(lh$L_a, type="l", lwd=4, xlab="Age", ylab="Length (cm)")
plot(lh$W_a, type="l", lwd=4, xlab="Age", ylab = "Weight (g)")
plot(lh$Mat_l, type="l", lwd=4, xlab="Length (cm)", ylab = "Proportion Mature")
plot(lh$S_fl[1,], type="l", lwd=4, xlab="Length (cm)", ylab="Proportion selected to gear")
#* Data inputs to lime ####

## Need to be in a matrix with year as the rows, length bins as the columns and the number in each length bin as the values 
labels <- seq(10,70,10)

S00 <- Ages.NTZ %>% 
  filter(Scenario %in% c("Historical and current management")) %>% 
  group_by(Mod_Year) %>% 
  mutate(Length = (Linf*(1-exp(-k*(Age-t0))))) %>% 
  ungroup() %>% 
  mutate(Bins = cut(.$Length, breaks=seq(0,70,10), labels=labels))

S00.S <- NULL

for(A in 1:30){
  temp <- S00 %>% 
    filter(Age == A)
  
  temp2 <- temp %>% 
    mutate(Mean_Pop = Mean_Pop * Selectivity[A,12])
  
  S00.S <- rbind(S00.S, temp2)
  
}


S00_Lime <- array(0, dim=c(59, 7))
years <- seq(1960,2018,1)
Lengths <- labels

for(Y in 1:59){
  for(L in 1:length(Lengths)){
    
    length <- labels[L]
    year <- years[Y]
    
    temp <- S00.S %>% 
      filter(Mod_Year==year) %>% 
      filter(Bins == length) %>% 
      summarise(sum(Mean_Pop))
    
    S00_Lime[Y,L] <- temp$`sum(Mean_Pop)`
    
  }
}


S00_Lime <- S00_Lime[11:59, ] %>% 
  replace(is.na(.), 0) %>% 
  as.matrix(.)

rownames(S00_Lime) <- seq(11,59,1)
colnames(S00_Lime) <- labels

## This is the example data
# true <- generate_data(modpath = NULL, itervec = 1, Fdynamics = c("Endogenous"),
#                       Rdynamics = "Constant", lh = lh, Nyears = 20, Nyears_comp = c(20), comp_sample = rep(200,20), init_depl = 0.7, seed = 123, fleet_proportions = 1)
# 

# Set up list with all the inputs 
data_all <- list(years=11:59, LF=S00_Lime)
inputs_all <- create_inputs(lh=lh, input_data = data_all) # Names of years in LF has to match years in data_all i.e. 1:59 not 1960-2018

#* Run LIME ####

lc_only <- run_LIME(modpath = NULL, input = inputs_all, data_avail = "LC", derive_quants=T)

gradient <- lc_only$opt$max_gradient <= 0.001 
hessian <- lc_only$Sdreport$pdHess # Needs to be positive and definite
hessian == TRUE & gradient == TRUE

plot_LCfits(LF_df = LF_df, Inputs = lc_only$Inputs, Report = lc_only$Report)

plot_output(Inputs = lc_only$Inputs, Report = lc_only$Report, Sdreport = lc_only$Sdreport,True =S00.1,  lh = lh, plot = c("Fish", "Rec", "SPR", "ML", "SB", "Selex"),
            set_ylim = list(Fish = c(0, 0.9), SPR = c(0, 1)))

lc_only$Derived
lc_only$Report$F_fy

#### RUN LIME FOR FISHED AREA ####
#* Data inputs to lime ####

## Need to be in a matrix with year as the rows, length bins as the columns and the number in each length bin as the values 
labels <- seq(10,70,10)

S00 <- Ages.F %>% 
  filter(Scenario %in% c("Historical and current management")) %>% 
  group_by(Mod_Year) %>% 
  mutate(Length = (Linf*(1-exp(-k*(Age-t0))))) %>% 
  ungroup() %>% 
  mutate(Bins = cut(.$Length, breaks=seq(0,70,10), labels=labels))

S00.S <- NULL

for(A in 1:30){
  temp <- S00 %>% 
    filter(Age == A)
  
  temp2 <- temp %>% 
    mutate(Mean_Pop = Mean_Pop * Selectivity[A,12])
  
  S00.S <- rbind(S00.S, temp2)
  
}

S00_Lime <- array(0, dim=c(59, length(labels)))
years <- seq(1960,2018,1)
Lengths <- labels

for(Y in 1:59){
  for(L in 1:length(Lengths)){
    
    length <- Lengths[L]
    year <- years[Y]
    
    temp <- S00.S %>% 
      filter(Mod_Year==year) %>% 
      filter(Bins == length) %>% 
      summarise(sum(Mean_Pop))
    
    S00_Lime[Y,L] <- temp$`sum(Mean_Pop)`
    
  }
}
S00_Lime <- S00_Lime[11:59, ] %>% 
  replace(is.na(.), 0) %>% 
  as.matrix(.)

rownames(S00_Lime) <- seq(11,59,1)
colnames(S00_Lime) <- Lengths

## This is the example data
# true <- generate_data(modpath = NULL, itervec = 1, Fdynamics = c("Endogenous"),
#                       Rdynamics = "Constant", lh = lh, Nyears = 20, Nyears_comp = c(20), comp_sample = rep(200,20), init_depl = 0.7, seed = 123, fleet_proportions = 1)
# 
# Set up list with all the inputs 
data_all <- list(years=11:59, LF=S00_Lime)
inputs_all <- create_inputs(lh=lh, input_data = data_all) # Names of years in LF has to match years in data_all i.e. 1:59 not 1960-2018
inputs_all$theta <- 50 # fix this because it's being estimated way too high
#* Run LIME ####

lc_only_F <- run_LIME(modpath = NULL, input = inputs_all, data_avail = "LC", derive_quants=T, est_selex_f = F)

gradient <- lc_only_F$opt$max_gradient <= 0.001 
hessian <- lc_only_F$Sdreport$pdHess # Needs to be positive and definite
hessian == TRUE & gradient == TRUE

plot_LCfits(LF_df = LF_df, Inputs = lc_only_F$Inputs, Report = lc_only_F$Report)

plot_output(Inputs = lc_only_F$Inputs, Report = lc_only_F$Report, Sdreport = lc_only_F$Sdreport,True =S00.1,  lh = lh, plot = c("Fish", "Rec", "SPR", "ML", "SB", "Selex"),
            set_ylim = list(Fish = c(0, 0.9), SPR = c(0, 1)))

lc_only_F$Derived
lc_only_F$Report$sigma_R

#### RUN LIME FOR WHOLE AREA ####

## Need to be in a matrix with year as the rows, length bins as the columns and the number in each length bin as the values 
labels <- seq(1,70,10)

S00 <- Ages.All %>%
  filter(Scenario %in% c("Historical and current management")) %>% 
  pivot_longer(cols = -c("Scenario", "Age"), values_to="Mean_Pop", names_to = "Mod_Year") %>% 
  mutate(Mod_Year = rep(seq(1960,2018,1), 30)) %>% 
  group_by(Mod_Year) %>% 
  mutate(Length = (Linf*(1-exp(-k*(Age-t0))))) %>% 
  ungroup() %>% 
  mutate(Bins = cut(.$Length, breaks=seq(0,70,10), label=labels))

S00.S <- NULL
# S00.S <- S00

for(A in 1:30){
  temp <- S00 %>% 
    filter(Age == A)
  
  temp2 <- temp %>% 
    mutate(Mean_Pop = Mean_Pop * Selectivity[A,12])
  
  S00.S <- rbind(S00.S, temp2)
  
}

S00_Lime <- array(0, dim=c(59, length(labels)))
years <- seq(1960,2018,1)
Lengths <- labels

for(Y in 1:59){
  for(L in 1:length(Lengths)){
    
    length <- Lengths[L]
    year <- years[Y]
    
    temp <- S00.S %>% 
      filter(Mod_Year==year) %>% 
      filter(Bins == length) %>% 
      summarise(sum(Mean_Pop))
    
    S00_Lime[Y,L] <- temp$`sum(Mean_Pop)`
    
  }
}

S00_Lime <- S00_Lime[1:58, ] %>% 
  replace(is.na(.), 0) %>% 
  as.matrix(.)

rownames(S00_Lime) <- seq(1,58,1)
colnames(S00_Lime) <- Lengths
## This is the example data
# true <- generate_data(modpath = NULL, itervec = 1, Fdynamics = c("Endogenous"),
#                       Rdynamics = "Constant", lh = lh, Nyears = 20, Nyears_comp = c(20), comp_sample = rep(200,20), init_depl = 0.7, seed = 123, fleet_proportions = 1)
# 

# Set up list with all the inputs 
data_all <- list(years=1:58, LF=S00_Lime)
inputs_all <- create_inputs(lh=lh, input_data = data_all) # Names of years in LF has to match years in data_all i.e. 1:59 not 1960-2018

#* Run LIME ####

lc_only_All <- run_LIME(modpath = NULL, input = inputs_all, data_avail = "LC", derive_quants=T)

gradient <- lc_only$opt$max_gradient <= 0.001 
hessian <- lc_only$Sdreport$pdHess # Needs to be positive and definite
hessian == TRUE & gradient == TRUE

plot_LCfits(LF_df = LF_df, Inputs = lc_only_All$Inputs, Report = lc_only_All$Report)

plot_output(Inputs = lc_only_All$Inputs, Report = lc_only_All$Report, Sdreport = lc_only_All$Sdreport, True =S00.1,  lh = lh, plot = c("Fish", "Rec", "SPR", "ML", "SB", "Selex"),
            set_ylim = list(Fish = c(0, 0.95), SPR = c(0, 1)))

lc_only_All$Derived$Fmsy
lc_only_All$Report$TB_t

#### CREATE KOBE PLOT #####

#* MSY plot ####
FMSY <- lc_only_All$Derived$Fmsy
BMSY <- lc_only_All$Derived$Bmsy

FMSY_All <- as.numeric(lc_only_All$Report$F_y)
# FMSY_All <- as.numeric(FMSY_All/FMSY) 

BBmsy_All <- as.numeric(lc_only_All$Report$TB_t)
# BBmsy_All <- as.numeric(BBmsy_All/BMSY) 

MSY_NTZ <- cbind(as.numeric(lc_only$Report$F_y), as.numeric(lc_only$Report$TB_t)) %>% 
  as.data.frame() %>%
  rename(FM = "V1",
         Bio = "V2") %>% 
  mutate(Mod_Year = seq(1970,2018,1)) %>%
  mutate(Rel.F = FM/lc_only$Derived$Fmsy) %>% 
  mutate(Rel.Bio = Bio/lc_only$Derived$Bmsy) %>%  
  mutate(Zone="NTZ")

MSY_F <- cbind(as.numeric(lc_only_F$Report$F_y), as.numeric(lc_only_F$Report$TB_t)) %>% 
  as.data.frame() %>% 
  rename(FM = "V1",
         Bio = "V2") %>% 
  mutate(Mod_Year = seq(1970,2018,1)) %>%
  mutate(Rel.F = FM/lc_only_F$Derived$Fmsy) %>% 
  mutate(Rel.Bio = Bio/lc_only_F$Derived$Bmsy) %>% 
  mutate(Zone="F") 

MSY_All <- rbind(as.data.frame(MSY_NTZ), as.data.frame(MSY_F)) %>% 
  mutate(Zone = ifelse(Mod_Year<=1987, "All", Zone)) %>% 
  mutate(FM = ifelse(Mod_Year <=1987, as.numeric(FMSY_All), FM)) %>% 
  mutate(Bio = ifelse(Mod_Year <=1987, as.numeric(BBmsy_All), Bio)) %>% 
  mutate(Rel.F = FM/FMSY) %>% 
  mutate(Rel.Bio = Bio/BMSY) %>%  
  mutate(Mod_Year = ifelse(Mod_Year %% 5 !=0, "", Mod_Year)) %>%  # %% means remainder of the division so here I only want years were dividing by 5 gives you a whole number
  as.data.frame()


MSY.Plot <- ggplot()+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 2.5, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 2.5,fill = "red")+
  annotate("rect", xmin = 1, xmax = 2.5, ymin = 1, ymax = 2.5,fill = "orange")+
  geom_point(data=MSY_All[c(1:49,67:98), ], aes(x=Rel.Bio, y=Rel.F, group=Zone))+
  geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_path(data=MSY_All[c(1:19),], aes(x=Rel.Bio, y=Rel.F, group=Zone))+
  geom_path(data=MSY_All[17:49,], aes(x=Rel.Bio, y=Rel.F))+
  geom_path(data=MSY_All[c(17,18, 67:98),], aes(x=Rel.Bio, y=Rel.F), linetype="dashed")+ # Fished
  geom_text(data=MSY_All[c(1:49,67:98), ], aes(x=Rel.Bio, y=Rel.F, group=Zone, label=Mod_Year), hjust=-0.25, vjust=-0.75, position=position_jitter(width=0.01,height=0.006))+
  scale_y_continuous(breaks=seq(0,2.5, 0.25), limits=c(0,2.5)) +
  scale_x_continuous(breaks=seq(0,2.5, 0.25), limits=c(0,2.5)) +
  theme_classic()+
  ylab("F/Fmsy")+
  xlab("B/Bmsy")
MSY.Plot

#* Conservation target plot ####

lime.weight <- lh$W_a
 
Bcon_All <- Ages.All %>% 
  filter(Scenario %in% c("Historical and current management")) %>% 
  dplyr::select(V1,Age) %>% 
  mutate(Weights = (V1*lime.weight)/1000) %>% 
  filter(!Age==1) %>% 
  filter(!Age==2) %>% 
  summarise(sum(Weights)) 
# Get F value from nearest biomass year in the LIME predictions
lc_only_All$Report$TB_t
lc_only_All$Report$N_t
Fcon_All <- 

Bcon_All <- Bcon_All$`sum(Weights)`*0.9

Fcon_NTZ
lc_only$Report$TB_t

Bcon_NTZ <- Ages.NTZ %>% 
  filter(Scenario %in% c("Historical and current management")) %>% 
  filter(Mod_Year == 1960) %>% 
  dplyr::select(Mean_Pop,Age) %>% 
  mutate(Weights = Mean_Pop*Weight[,12]) %>% 
  summarise(sum(Weights)) 
Bcon_NTZ <- Bcon_NTZ$`sum(Weights)`*0.9

Fcon_F <- 0.07770899
lc_only_F$Report$TB_t # Gives us year 21
Bcon_F <- Ages.F %>% 
  filter(Scenario %in% c("Historical and current management")) %>% 
  filter(Mod_Year == 1960) %>% 
  dplyr::select(Mean_Pop,Age) %>% 
  mutate(Weights = Mean_Pop*Weight[,12]) %>% 
  summarise(sum(Weights)) 
Bcon_F <- Bcon_F$`sum(Weights)`*0.9
  



#### Checking things ####
ages <- lh$ages
M <- lh$M
S_fa <- lh$S_fa
R0 <- lh$R0
W_a <- lh$W_a

#* calc_equil_abundance
calc_equil_abund <- function(ages, M, F, S_fa, R0){
  N_a <- rep(NA, length(ages))
  N_a[1] <- R0
  
  nfleets <- nrow(S_fa)
  Fmat <- matrix(NA, nrow=nfleets, ncol=length(ages))
  for(i in 1:nfleets){
    Fmat[i,] <- S_fa[i,]*F[i]
  }
  Ftotal <- colSums(Fmat)
  
  for(i in 2:length(ages)){
    if(i<length(ages)) N_a[i] <- N_a[i-1]*exp(-M-Ftotal[i-1])
    if(i==length(ages)) N_a[i] <- (N_a[i-1]*exp(-M-Ftotal[i-1]))/(1-exp(-M-Ftotal[i-1]))
  }
  # outputs <- list()
  # outputs[[1]] <- N_a
  # outputs[[2]] <- Fmat
  # return(outputs)
  return(N_a)
  
}

FM <- seq(0,0.5,0.1)

check <- NULL
num_age <- NULL
for (i in 1:6){
  fm <- FM[i]
  equil.abund <- calc_equil_abund(ages, M, F=fm, S_fa, R0)
  check <- rbind(check,equil.abund[[2]])
  num_age <- rbind(num_age, equil.abund[[1]])
}

#* calcMSY
calc_msy <- function(F, ages, M, R0, W_a, S_fa){
  
  Nage <- calc_equil_abund(ages=ages, M=M, F=F, R0=R0, S_fa=S_fa)
  YPR <- sum(Nage * W_a * (1 - exp(-M - F*S_fa[1,])) * (F*S_fa[1,]) / (M + F*S_fa[1,]))
  
  return(YPR)
}

MSY.calc <- calc_msy(F=0.3, ages,  M, R0, W_a, S_fa)

plot(x=ages, y=S_fa, type="line")
plot(x=ages, y=W_a, type="l")

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
Linf = 664 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
k = 0.241 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
t0 = -0.375 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr

## Weight-Length Relationship
WLa = 0.000028
WLb = 2.8761

## Natural Mortality
M = 0.146  #Marriot et al. 2011
step <- 1/12

# Beverton-Holt Parameters
h = 0.76 # Steepness 0.76
R0 = 1 # Initial recruitment

# Selectivity
A50 = 3 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
A95 = 4 # https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr

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


#* Convert age data to length data (mm)  and work out maturity/selectivity in length ####
Life.History <- as.data.frame(array(1, dim=c(361, 4))) %>% 
  rename(Age = "V1") %>% 
  rename(Length = "V2") %>% 
  rename(Maturity = "V3") %>% 
  rename(Selecivity = "V4")

for (r in 1:360){                                   
  Life.History[r+1, 1] <- Life.History[r, 1] + step 
}

Life.History <- Life.History %>% 
  mutate(Length = Linf*(1-exp(-k*(Age-t0)))) %>% 
  mutate(Selectivity = 1/(1+(exp(-log(19)*((Age-A50)/(A95-A50)))))) %>% 
  mutate(Maturity = 1/(1+(exp(-log(19)*((Age-M50)/(M95-M50))))))

## Length of every age class
Length <- Life.History$Length
Length <- array(Length, dim=c(12,30))
Length <- t(Length)

# Selectivity for every age class
Selectivity <- Life.History$Selectivity
Selectivity <- array(Selectivity, dim=c(12,30))
Selectivity <- t(Selectivity)
# S50 369.6095
# S95

## Proportion mature in every age class
Maturity <- Life.History$Maturity
Maturity <- array(Maturity, dim=c(12,30))
Maturity <- t(Maturity)
# M50
#M95

#### LIME ####
#* Specify biological inputs and starting values ####
lh <- create_lh_list(vbk = 0.241,
                     linf = 66.4,
                     t0 = -0.375,
                     lwa = 0.000028,
                     lwb = 2.8761,
                     S50 = c(36.96095),
                     S95 = c(43.26557),
                     selex_input = "length",
                     selex_type = c("logistic"),
                     M50 = 46.69928, # Is actually a little bit less than this but the maturity grid doesn't give 0.5 exactly, you get 0.51
                     M95 = 52.11340,
                     h = 0.76, # matches what I have in my population set up
                     AgeMax = 30,
                     maturity_input = "length",
                     M = 0.146,
                     binwidth = 1,
                     CVlen = 0.1, # Variability around the growth curve
                     SigmaR = 0.6, # SD aroud the recruitment variation
                     SigmaF = 0.2, # SD around fishing mortality
                     SigmaC = 0.1, # SD around the catch although I won't put catch data in
                     R0 = 1, # Initial unfished recruitment
                     Frate = 2, # exponent determining the strength of coupling between effort and biomass
                     Fequil = 0.25, # Prop.  SB to Unfished SB at bioeconomic equilib.
                     #qcoef = 1e-06, # Catchability - changed it to be the same as the model
                     start_ages = 0,
                     rho = 0,
                     theta = 1,
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
labels <- seq(5,70,5)

S00 <- Ages.NTZ %>% 
  filter(Scenario %in% c("Historical and current management")) %>% 
  group_by(Mod_Year) %>% 
  mutate(Length = (Length[,1])/10) %>% 
  ungroup() %>% 
  mutate(Bins = cut(.$Length, breaks=seq(0,70,5), labels=labels))

S00.S <- NULL

for(A in 1:30){
  temp <- S00 %>% 
    filter(Age == A)
  
  temp2 <- temp %>% 
    mutate(Mean_Pop = Mean_Pop * Selectivity[A,1])
  
  S00.S <- rbind(S00.S, temp2)
  
}

S00_Lime <- S00.S %>%
  pivot_wider(names_from = Bins, values_from = Mean_Pop) %>% 
  filter(Age>2) %>% 
  group_by(Mod_Year) %>% 
  mutate_all(funs(.[order(is.na(.))])) %>% 
  filter_at(vars(6:14), any_vars(!is.na(.))) %>% 
  ungroup()

S00_Lime <- S00_Lime[11:59, 6:14] %>% 
  replace(is.na(.), 0) %>% 
  as.matrix(.)

rownames(S00_Lime) <- seq(11,59,1)

## This is the example data
# true <- generate_data(modpath = NULL, itervec = 1, Fdynamics = c("Endogenous"),
#                       Rdynamics = "Constant", lh = lh, Nyears = 20, Nyears_comp = c(20), comp_sample = rep(200,20), init_depl = 0.7, seed = 123, fleet_proportions = 1)
# 

# Doesn't look great but I'm not sure there's much I can do about this just because of the way the lengths have ended up
LF_df <- LFreq_df(S00_Lime)
plot_LCfits(LF_df = LF_df)

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

#### RUN LIME FOR FISHED AREA ####
#* Data inputs to lime ####

## Need to be in a matrix with year as the rows, length bins as the columns and the number in each length bin as the values 
labels <- seq(5,70,5)

S00 <- Ages.F %>% 
  filter(Scenario %in% c("Historical and current management")) %>% 
  group_by(Mod_Year) %>% 
  mutate(Length = (Length[,1])/10) %>% 
  ungroup() %>% 
  mutate(Bins = cut(.$Length, breaks=seq(0,70,5), labels=labels))

S00.S <- NULL

for(A in 1:30){
  temp <- S00 %>% 
    filter(Age == A)
  
  temp2 <- temp %>% 
    mutate(Mean_Pop = Mean_Pop * Selectivity[A,1])
  
  S00.S <- rbind(S00.S, temp2)
  
}

S00_Lime <- S00.S %>%
  pivot_wider(names_from = Bins, values_from = Mean_Pop) %>% 
  filter(Age>2) %>% 
  group_by(Mod_Year) %>% 
  mutate_all(funs(.[order(is.na(.))])) %>% 
  filter_at(vars(6:14), any_vars(!is.na(.))) %>% 
  ungroup()

S00_Lime <- S00_Lime[11:59, 6:14] %>% 
  replace(is.na(.), 0) %>% 
  as.matrix(.)

rownames(S00_Lime) <- seq(11,59,1)

## This is the example data
# true <- generate_data(modpath = NULL, itervec = 1, Fdynamics = c("Endogenous"),
#                       Rdynamics = "Constant", lh = lh, Nyears = 20, Nyears_comp = c(20), comp_sample = rep(200,20), init_depl = 0.7, seed = 123, fleet_proportions = 1)
# 

# Doesn't look great but I'm not sure there's much I can do about this just because of the way the lengths have ended up
LF_df <- LFreq_df(S00_Lime)
plot_LCfits(LF_df = LF_df)

# Set up list with all the inputs 
data_all <- list(years=11:59, LF=S00_Lime)
inputs_all <- create_inputs(lh=lh, input_data = data_all) # Names of years in LF has to match years in data_all i.e. 1:59 not 1960-2018

#* Run LIME ####

lc_only_F <- run_LIME(modpath = NULL, input = inputs_all, data_avail = "LC", derive_quants=T)

gradient <- lc_only_F$opt$max_gradient <= 0.001 
hessian <- lc_only_F$Sdreport$pdHess # Needs to be positive and definite
hessian == TRUE & gradient == TRUE

plot_LCfits(LF_df = LF_df, Inputs = lc_only_F$Inputs, Report = lc_only_F$Report)

plot_output(Inputs = lc_only_F$Inputs, Report = lc_only_F$Report, Sdreport = lc_only_F$Sdreport,True =S00.1,  lh = lh, plot = c("Fish", "Rec", "SPR", "ML", "SB", "Selex"),
            set_ylim = list(Fish = c(0, 0.9), SPR = c(0, 1)))

lc_only_F$Derived

#### RUN LIME FOR WHOLE AREA ####

## Need to be in a matrix with year as the rows, length bins as the columns and the number in each length bin as the values 
labels <- seq(5,70,5)

S00 <- Ages.All %>%
  filter(Scenario %in% c("Historical and current management")) %>% 
  pivot_longer(cols = -c("Scenario", "Age"), values_to="Mean_Pop", names_to = "Mod_Year") %>% 
  group_by(Mod_Year) %>% 
  mutate(Length = (Length[,1])/10) %>% 
  ungroup() %>% 
  mutate(Bins = cut(.$Length, breaks=seq(0,70,5), labels=labels))

S00.S <- NULL

for(A in 1:30){
  temp <- S00 %>% 
    filter(Age == A)
  
  temp2 <- temp %>% 
    mutate(Mean_Pop = Mean_Pop * Selectivity[A,1])
  
  S00.S <- rbind(S00.S, temp2)
  
}


S00_Lime <- S00.S %>%
  pivot_wider(names_from = Bins, values_from = Mean_Pop) %>% 
  filter(Age>2) %>% 
  group_by(Mod_Year) %>% 
  mutate_all(funs(.[order(is.na(.))])) %>% 
  filter_at(vars(5:13), any_vars(!is.na(.))) %>% 
  ungroup()

S00_Lime <- S00_Lime[11:27, 5:13] %>% 
  replace(is.na(.), 0) %>% 
  as.matrix(.)

rownames(S00_Lime) <- seq(11,27,1)

## This is the example data
# true <- generate_data(modpath = NULL, itervec = 1, Fdynamics = c("Endogenous"),
#                       Rdynamics = "Constant", lh = lh, Nyears = 20, Nyears_comp = c(20), comp_sample = rep(200,20), init_depl = 0.7, seed = 123, fleet_proportions = 1)
# 

# Doesn't look great but I'm not sure there's much I can do about this just because of the way the lengths have ended up
LF_df <- LFreq_df(S00_Lime)
plot_LCfits(LF_df = LF_df)

# Set up list with all the inputs 
data_all <- list(years=11:27, LF=S00_Lime)
inputs_all <- create_inputs(lh=lh, input_data = data_all) # Names of years in LF has to match years in data_all i.e. 1:59 not 1960-2018

#* Run LIME ####

lc_only_All <- run_LIME(modpath = NULL, input = inputs_all, data_avail = "LC", derive_quants=T)

gradient <- lc_only$opt$max_gradient <= 0.001 
hessian <- lc_only$Sdreport$pdHess # Needs to be positive and definite
hessian == TRUE & gradient == TRUE

plot_LCfits(LF_df = LF_df, Inputs = lc_only$Inputs, Report = lc_only$Report)

plot_output(Inputs = lc_only$Inputs, Report = lc_only$Report, Sdreport = lc_only$Sdreport,True =S00.1,  lh = lh, plot = c("Fish", "Rec", "SPR", "ML", "SB", "Selex"),
            set_ylim = list(Fish = c(0, 0.95), SPR = c(0, 1)))

lc_only_All$Derived
lc_only_All$Report$TB_t

#### CREATE KOBE PLOT #####
FMSY <- 0.4165581
BMSY <- 3.71755

FMSY_All <- as.numeric(lc_only_All$Report$F_y)
# FMSY_All <- as.numeric(FMSY_All/FMSY) 

BBmsy_All <- as.numeric(lc_only_All$Report$TB_t)
# BBmsy_All <- as.numeric(BBmsy_All/BMSY) 

MSY_NTZ <- cbind(as.numeric(lc_only$Report$F_y), as.numeric(lc_only$Report$TB_t)) %>% 
  as.data.frame() %>%
  rename(FM = "V1",
         Bio = "V2") %>% 
  mutate(Mod_Year = seq(1970,2018,1)) %>%
  mutate(Rel.F = FM/FMSY) %>% 
  mutate(Rel.Bio = Bio/BMSY) %>%  
  mutate(Zone="NTZ")

MSY_F <- cbind(as.numeric(lc_only_F$Report$F_y), as.numeric(lc_only_F$Report$TB_t)) %>% 
  as.data.frame() %>% 
  rename(FM = "V1",
         Bio = "V2") %>% 
  mutate(Mod_Year = seq(1970,2018,1)) %>%
  mutate(Rel.F = FM/FMSY) %>% 
  mutate(Rel.Bio = Bio/BMSY) %>% 
  mutate(Zone="F") 

MSY_All <- rbind(as.data.frame(MSY_NTZ), as.data.frame(MSY_F)) %>% 
  mutate(Zone = ifelse(Mod_Year<1987, "All", Zone)) %>% 
  mutate(FM = ifelse(Mod_Year < 1987, as.numeric(FMSY_All), FM)) %>% 
  mutate(Bio = ifelse(Mod_Year < 1987, as.numeric(BBmsy_All), Bio)) %>% 
  mutate(Rel.F = FM/FMSY) %>% 
  mutate(Rel.Bio = Bio/BMSY) %>%  
  mutate(Mod_Year = ifelse(Mod_Year %% 5 !=0, "", Mod_Year)) # %% means remainder of the division so here I only want years were dividing by 5 gives you a whole number
  as.data.frame()


Bio.Plot <- ggplot()+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,fill = "yellow")+
  annotate("rect", xmin = 1, xmax = 2, ymin = 0, ymax = 1,fill = "green")+
  annotate("rect", xmin = 0, xmax = 1, ymin = 1, ymax = 2,fill = "red")+
  annotate("rect", xmin = 1, xmax = 2, ymin = 1, ymax = 2,fill = "orange")+
  geom_point(data=MSY_All[c(1:49,67:98), ], aes(x=Rel.Bio, y=Rel.F, group=Zone))+
  geom_text(hjust=0, vjust=0, position=position_jitter(width=0.01,height=0.006))+
  geom_path(data=MSY_All[c(1:19),], aes(x=Rel.Bio, y=Rel.F, group=Zone))+
  geom_path(data=MSY_All[17:49,], aes(x=Rel.Bio, y=Rel.F))+
  geom_path(data=MSY_All[c(17,18, 67:98),], aes(x=Rel.Bio, y=Rel.F), linetype="dashed")+ # Fished
  geom_text(data=MSY_All[c(1:49,67:98), ], aes(x=Rel.Bio, y=Rel.F, group=Zone, label=Mod_Year), hjust=-0.25, vjust=-0.75, position=position_jitter(width=0.01,height=0.006))+
  scale_y_continuous(breaks=seq(0,2, 0.25), limits=c(0,2)) +
  scale_x_continuous(breaks=seq(0,2, 0.25), limits=c(0,2)) +
  theme_classic()+
  ylab("F/Fmsy")+
  xlab("B/Bmsy")
Bio.Plot



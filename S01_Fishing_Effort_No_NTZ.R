###################################################

# Script for setting up the fishing surface for
# simulation 1, with no sanctuary zones
# Requires files that are made in the first script
# Also requires fishing effort data

###################################################


## NEED TO CHANGE THE INITIAL FISHING MORTALITY AT 1965 TO BE THE SAME AS THE EQLUILIBRIUM LEVEL ##
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

#### LOAD FILES ####

st_centroid_within_poly <- function (poly) { #returns true centroid if inside polygon otherwise makes a centroid inside the polygon
  
  # check if centroid is in polygon
  centroid <- poly %>% st_centroid() 
  in_poly <- st_within(centroid, poly, sparse = F)[[1]] 
  
  # if it is, return that centroid
  if (in_poly) return(centroid) 
  
  # if not, calculate a point on the surface and return that
  centroid_in_poly <- st_point_on_surface(poly) 
  return(centroid_in_poly)
}

## Data
setwd(data_dir)
boat_days <- read.csv("Boat_Days_Gascoyne.csv")

boat_days <- boat_days%>%
  mutate(NumMonth = as.numeric(NumMonth)) %>% 
  mutate(Month = as.factor(Month)) %>% 
  mutate(Gascoyne_Boat_Days = as.numeric(Gascoyne_Boat_Days)) %>% 
  mutate(Month = fct_relevel(Month, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))

## Spatial Data
setwd(sg_dir)

water <- readRDS("water")
NCELL <- nrow(water)

# Locations of the boat ramps
setwd(sp_dir)
BR <- st_read("Boat_Ramps.shp") %>% 
  st_transform(4283)%>%
  st_make_valid

# Parameter values
q <- 0.005

#### SET UP FISHING SURFACE ####

# Work out the number of boat days for the Exmouth region
# For this we're using data from the whole of the Gascoyne and then splitting up the number of boat days in Exmouth compared
# to Carnarvon/Shark Bay and then using visitor numbers to work out what percentage of the baot days were in Exmouth
# LGA statistics for visitors to the biggest shires where people would go fishing in the Gascoyne indicate that trips to
# Exmouth are about 33% of the visitation to the region so we'll allocate 33% of the Boat Days to Exmouth
# We then need to create a model for hindcasting fishing effort back to something like the 1960s to get an estimate of what 
# effort might have been like for the years where we don't have any data 

## Working out the fishing effort in the Gascoyne
# Plotting the data to see what it looks like
TotalYear <- boat_days %>% 
  group_by(Year) %>% 
  summarise(., Total_Boat_Days=sum(Gascoyne_Boat_Days))

YearPlot <- ggplot(TotalYear) + 
  geom_point(aes(x=Year, y=Total_Boat_Days))

YearModel <- lm(Total_Boat_Days~Year, data=TotalYear)
summary(YearModel)

#### HINDCASTING ####
Year2011_1990 <- as.data.frame(array(0, dim=c(21,1))) %>% 
  mutate(V1=seq(1990, 2010, by=1)) %>% 
  rename("Year"=V1)

predictions <- predict(YearModel, newdata=Year2011_1990)

Year2011_1990 <- Year2011_1990 %>% 
  mutate(Total_Boat_Days = predictions)

boat_days_hind <- rbind(TotalYear, Year2011_1990)

effort <- seq(0, 118466, length=30)
years <- seq(1960, 1989, by=1)

Years_1960_1989 <- as.data.frame(cbind(years, effort)) %>% 
  rename("Year" = years) %>% 
  rename("Total_Boat_Days" = effort)

#### JOIN FISHING EFFORT ####
boat_days_hind <- rbind(boat_days_hind, Years_1960_1989)

YearPlot <- ggplot(boat_days_hind) + 
  geom_line(aes(x=Year, y=Total_Boat_Days))+ 
  theme_classic()+
  ylab("Total Boat Days")

#### ALLOCATING MONLTHY EFFORT ####

# Work out proportion of fishing effort in each month for each year
boat_days <- boat_days %>% 
  group_by(Year) %>% 
  mutate(Year_Sum = sum(Gascoyne_Boat_Days)) %>%
  mutate(Month_Prop = Gascoyne_Boat_Days/Year_Sum) %>% 
  dplyr::select(-Year_Sum)

# Work out the average fishing effort for each month
boat_days <- boat_days %>% 
  group_by(Month) %>% 
  mutate(Ave_Month_Prop = mean(Month_Prop))

Month_Prop_Ave <- boat_days[1:12,c(2,6)]

check <- boat_days %>% 
  group_by(Year) %>% 
  summarise(Ave_Month_Prop = sum(Ave_Month_Prop))

check <- boat_days %>% 
  group_by(Year) %>% 
  summarise(Ave_Month = mean(Gascoyne_Boat_Days))

# Split up the hindcast data into monthly means 
boat_days_hind <- boat_days_hind %>% 
  bind_rows(replicate(11, boat_days_hind, simplify = FALSE)) %>% # Make a row for each month of each year in the hindcast
  mutate(Year = as.integer(Year)) %>% 
  arrange(-Year) %>% 
  mutate(NumMonth = (rep(1:12, 59))) %>% # Add numerical month
  mutate(Month = rep(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), 59)) %>% #Add word month
  filter(Year<2011) 


boat_days_hind <- boat_days_hind %>% 
  mutate(Gascoyne_Boat_Days = 0) %>% 
  left_join(., Month_Prop_Ave) %>% 
  group_by(Year) %>% 
  mutate(Gascoyne_Boat_Days = Total_Boat_Days*Ave_Month_Prop) %>% 
  ungroup() %>% 
  dplyr::select(-c(Total_Boat_Days))

check <- boat_days_hind %>% 
  group_by(Year) %>% 
  summarise(., sum(Gascoyne_Boat_Days))

# Put it all together into one big dataframe

Full_Boat_Days <- boat_days[,1:4] %>% 
  rbind(boat_days_hind) %>% 
  arrange(-Year)


# Plot and check that it looks right
MonthPlot <- Full_Boat_Days %>%
  mutate(Unique = paste(Year, NumMonth, sep="_")) %>%
  filter(Month %in% c("Jun")) %>%
  ggplot() +
  geom_point(aes(x=Unique, y=Gascoyne_Boat_Days))

#### ALLOCATION EFFORT TO EXMOUTH AND BOAT RAMPS ####
## We're using tourism numbers to allocate to Exmouth - saying about 65% of the trips to the region are to Ningaloo

# Create data frame that is just Exmouth Boat Days 
Exmouth_Boat_Days <- Full_Boat_Days %>% 
  mutate(Boat_Days = Gascoyne_Boat_Days*0.65) %>% 
  dplyr::select(-Gascoyne_Boat_Days)

# Split up the effort to Boat Ramps according to the information we have collected from Exmouth about how often people
# Use the boat ramps - the Exmouth Marina wasn't constructed until 1997 but there was a ramp at town beah from the 60s 
# onwards, this was built at the same time that the Tantabiddi ramp was built 
# Bundegi boat ramp wasn't built until 2008 and before that was just a beach with some concrete covered in sand
# So might be a good idea to reduce the number of people we assume launched from there as it was effectively a beach launch
# You also need to make sure that you standardise by number of hours spent at each boat ramp to account for sampling effort

# Tantabidi had 149.48 hours of sampling effort with 224 trips
# Bundegi had 133.9 hours of sampling effort with 157
# Exmouth Marina had 166.15 hours of sampling effort with 198 trips
# Coral Bay had 146.07 hours of sampling effort with 151 trips

BR_Trips <- data.frame(BoatRamp = c("Tantabiddi", "Bundegi", "ExmouthMar", "CoralBay"),
                       Effort = c(194.48, 133.9, 166.15, 146.07),
                       Trips = c(224, 157, 198, 151))

BR_Trips <- BR_Trips %>% 
  mutate(Trip_per_Hr = as.numeric(unlist((Trips/Effort)))) %>% # Standardise the no. trips based on how much time you spent sampling
  mutate(BR_Prop = Trip_per_Hr/sum(Trip_per_Hr)) %>% #Then work out the proportion of trips each hour that leave from each boat ramp
  mutate(BR_Prop_08 = c(0.3031544, 0.1577101, 0.3119251, 0.2272104)) # Create separate proportions for Bundegi before 2008, have 
# allocated 10% of its boats to the other Exmouth Boat Ramps

# Create a loop that allocates the correct proportion of boat effort to each of the BRs accounting for reduced effort at Bundegi before ther ramp was built
for(Y in 1:708){
  if(Exmouth_Boat_Days[Y,1]<2008){ # This is for when Bundegi wasn't a proper ramp so probably would have had less effort
    Exmouth_Boat_Days$Tb_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[1,6]
    Exmouth_Boat_Days$Bd_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[2,6]
    Exmouth_Boat_Days$ExM_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[3,6] 
    Exmouth_Boat_Days$CrB_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[4,6]
  } else {
    Exmouth_Boat_Days$Tb_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[1,5] 
    Exmouth_Boat_Days$Bd_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[2,5] 
    Exmouth_Boat_Days$ExM_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[3,5] 
    Exmouth_Boat_Days$CrB_BR = Exmouth_Boat_Days$Boat_Days*BR_Trips[4,5]
  }
}

check <- Exmouth_Boat_Days %>% 
  mutate(Total = Tb_BR+Bd_BR+ExM_BR+CrB_BR)

#### CONVERT NO-TAKE INTO NUMERIC VECTORS OF CELL IDS ####

# Creating new columns for each of our year groups when new NTZs were put in place

NoTake <- st_sf(water) %>% 
  st_drop_geometry() %>% 
  dplyr::select(Fished, ID, COMMENTS)

NoTake <- NoTake %>% 
  rename(Fished60_87 = "Fished") %>% 
  mutate(Fished60_87 = "Y") %>% 
  mutate(Fished87_05 = ifelse(str_detect(COMMENTS, c("Old")), "N", "Y")) %>% 
  mutate(Fished05_18 = ifelse(str_detect(COMMENTS, "[:alpha:]") & !str_detect(COMMENTS, ("Comm Cloates")), "N", "Y")) %>% 
  mutate(Fished18_21 = ifelse(str_detect(COMMENTS, "[:alpha:]"), "N", "Y")) %>% 
  dplyr::select(ID, Fished60_87, Fished87_05, Fished05_18, Fished18_21) %>% 
  mutate_at(vars(contains("Fished")), ~replace(., is.na(.), "Y")) %>% 
  mutate(ID = as.numeric(ID))

rownames(NoTake) <- seq(from = 1, to = 1901)

NoTake87_05 <- NoTake %>% 
  filter(Fished87_05=="N") %>% 
  mutate(ID = ID-1) %>%  # All of the cells are greater than 299 which is the one that I removed in an earlier script so they need to move back one position
  dplyr::select(ID)

NoTake05_18 <- NoTake %>% 
  filter(Fished05_18=="N") %>% 
  mutate(ID = ID-1) %>%
  dplyr::select(ID)

NoTake18_21 <- NoTake %>% 
  filter(Fished18_21=="N") %>% 
  mutate(ID = ID-1) %>%
  dplyr::select(ID)

NoTake <- list()

NoTake[[1]] <- as.numeric(NoTake87_05$ID)
NoTake[[2]] <- as.numeric(NoTake05_18$ID)
NoTake[[3]] <- as.numeric(NoTake18_21$ID)

#### ALLOCATING EFFORT TO CELLS ####

## Work out the probability of visiting a cell from each boat ramp based on distance and size
BR <- st_sf(BR)

water <- water %>% 
  mutate(DistBR = 0)

centroids <- st_centroid_within_poly(water)

# Working out distance from each BR to each cell
DistBR <- as.data.frame(array(0, dim = c(NCELL,3))) # This will contain the distance from each boat ramp to every cell

for(CELL in 1:NCELL){
  
  for(RAMP in 1:4){
    x <- st_distance(centroids[CELL,1], BR[RAMP,3])
    DistBR[CELL,RAMP] <- (x/1000) #to get the distance in km
    
  }
}

# Create a data frame with both the distances and the areas of the cells
Cell_Vars <- DistBR %>% 
  mutate(Area = as.vector((water$cell_area)/100000)) %>% #Cells are now in km^2 but with no units
  rename("Bd_BR"=V1) %>% 
  rename("ExM_BR" = V2) %>% 
  rename("Tb_BR" = V3) %>% 
  rename("CrB_BR"=V4)

## Whole dataframe with no NTZs 
BR_U <- as.data.frame(matrix(0, nrow=NCELL, ncol=4)) #Set up data frame to hold utilities of cells

Tot <- array(0, dim=c(NCELL, 4))
for(RAMP in 1:4){
  
  Tot[ ,RAMP] <- (Cell_Vars$Area/Cell_Vars[,RAMP])
}
Tot <- colSums(Tot)

for(RAMP in 1:4){
  for(CELL in 1:NCELL){
    BR_U[CELL,RAMP] <- ((Cell_Vars[CELL,5]/Cell_Vars[CELL,RAMP])/Tot[RAMP])
  }
} 
colSums(BR_U)

# Now we have BR_U which has the "utilities" for each cells based on it's size and distance from BR  

BR_Trips <- Exmouth_Boat_Days %>% # This is just the trips from each boat ramp
  ungroup() %>% 
  mutate(NumYear = rep(59:1, each=12)) %>% #This is to turn the years into a count for the loop
  dplyr::select(NumYear, NumMonth, Bd_BR, ExM_BR, Tb_BR, CrB_BR)  

Fishing <- array(0, dim=c(NCELL, 12, 59)) #This array has a row for every cell, a column for every month, and a layer for every year
Months <- array(0, dim=c(NCELL, 12))
Ramps <- array(0, dim=c(NCELL, 4))
layer <- 1

for(YEAR in 1:59){
  
  for(MONTH in 1:12){
    
    for(RAMP in 1:4){
      
      temp <- BR_Trips %>% 
        filter(NumYear==YEAR) %>% 
        dplyr::select(-c(NumYear, NumMonth))
      
      temp <- as.matrix(temp) 
      
      for(CELL in 1:NCELL){
        Ramps[CELL,RAMP] <- BR_U[CELL,RAMP] * temp[MONTH,RAMP]
      }
    }
    
    Months[,MONTH] <-  rowSums(Ramps)
  }
  Fishing[ , ,layer] <- Months 
  layer <- layer+1
}

Fishing <- Fishing*q # multiply by catchability

#### SAVE DATA ####
setwd(sim_dir)
saveRDS(Fishing, file="sim01_fishing")



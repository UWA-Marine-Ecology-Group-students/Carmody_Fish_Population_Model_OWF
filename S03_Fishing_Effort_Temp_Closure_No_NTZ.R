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

###################################################

# Script for setting up the fishing surface 
# Requires files that are made in the first script
# Also requires fishing effort data

rm(list = ls())

#### SET DIRECTORIES ####
working.dir <- dirname(rstudioapi::getActiveDocumentContext()$path) # to directory of current file - or type your own

data_dir <- paste(working.dir, "Data", sep="/")
fig_dir <- paste(working.dir, "Figures", sep="/")
m_dir <- paste(working.dir, "Matrices", sep="/")
sp_dir <- paste(working.dir, "Spatial_Data", sep="/")
sg_dir <- paste(working.dir, "Staging", sep="/")
sim_dir <- paste(working.dir, "Simulations", sep="/")

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

model.name <- "ningaloo"

## Data
setwd(data_dir)
boat_days <- read.csv("Boat_Days_Ningaloo.csv")

boat_days <- boat_days%>%
  mutate(NumMonth = as.numeric(NumMonth)) %>%
  mutate(Month = as.factor(Month)) %>%
  mutate(Survey_Year = as.character(Survey_Year)) %>% 
  mutate(Boat_Days_State = as.numeric(Boat_Days_State)) %>%
  mutate(Boat_Days_Commonwealth = as.numeric(Boat_Days_Commonwealth)) %>%
  mutate(Month = fct_relevel(Month, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))

## Spatial Data
setwd(sp_dir)

water <- readRDS(paste0(model.name, sep="_","water"))
NCELL <- nrow(water)

# Locations of the boat ramps
BR <- st_read("Boat_Ramps.shp") %>% 
  st_transform(4283)%>%
  st_make_valid

# No Take Zones
setwd(sg_dir)
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))

#### FILL IN MISSING VALUES FOR NINGALOO BOAT DAYS ####

# Work out the number of boat days for the Ningaloo region split by state and commonwealth
# We have estimates for just this area from Claire Smallwood but not for all years 
# We have calculated proportions for each month based on the Gascoyne data 
# So we are going to interpolate the months we are missing 

# Yes I'm aware that this code is horrendous and repetitive, again I had intended to put it into a loop and never got around to it... :')

## Interpolate the values for state and commonwealth separately
# January
boat_days_jan <- boat_days %>% 
  filter(NumMonth==1)

inter_jan_state <- approx(boat_days_jan$Year, boat_days_jan$Boat_Days_State, xout=c(2013,2015,2017)) %>% 
  unlist()
inter_jan_comm <- approx(boat_days_jan$Year, boat_days_jan$Boat_Days_Commonwealth, xout=c(2013,2015,2017)) %>% 
  unlist()

# February
boat_days_feb <- boat_days %>% 
  filter(NumMonth==2)

inter_feb_state <- approx(boat_days_feb$Year, boat_days_feb$Boat_Days_State, xout=c(2013,2015,2017)) %>% 
  unlist()
inter_feb_comm <- approx(boat_days_feb$Year, boat_days_feb$Boat_Days_Commonwealth, xout=c(2013,2015,2017)) %>% 
  unlist()

# March
boat_days_mar <- boat_days %>% 
  filter(NumMonth==3)

inter_mar_state <- approx(boat_days_mar$Year, boat_days_mar$Boat_Days_State, xout=c(2012,2013,2015,2017)) %>% 
  unlist()
inter_mar_comm <- approx(boat_days_mar$Year, boat_days_mar$Boat_Days_Commonwealth, xout=c(2012,2013,2015,2017)) %>% 
  unlist()

# April
boat_days_apr <- boat_days %>% 
  filter(NumMonth==4)

inter_apr_state <- approx(boat_days_apr$Year, boat_days_apr$Boat_Days_State, xout=c(2012,2013,2015,2017)) %>% 
  unlist()
inter_apr_comm <- approx(boat_days_apr$Year, boat_days_apr$Boat_Days_Commonwealth, xout=c(2012,2013,2015,2017)) %>% 
  unlist()

# May
boat_days_may <- boat_days %>% 
  filter(NumMonth==5)

inter_may_state <- approx(boat_days_may$Year, boat_days_may$Boat_Days_State, xout=c(2012,2014,2015,2017)) %>% 
  unlist()
inter_may_comm <- approx(boat_days_may$Year, boat_days_may$Boat_Days_Commonwealth, xout=c(2012,2014,2015,2017)) %>% 
  unlist()

# June
boat_days_jun <- boat_days %>% 
  filter(NumMonth==6)

inter_jun_state <- approx(boat_days_jun$Year, boat_days_jun$Boat_Days_State, xout=c(2012,2014,2015,2017)) %>% 
  unlist()
inter_jun_comm <- approx(boat_days_jun$Year, boat_days_jun$Boat_Days_Commonwealth, xout=c(2012,2014,2015,2017)) %>% 
  unlist()

# July
boat_days_jul <- boat_days %>% 
  filter(NumMonth==7)

inter_jul_state <- approx(boat_days_jul$Year, boat_days_jul$Boat_Days_State, xout=c(2012,2014,2015,2017)) %>% 
  unlist()
inter_jul_comm <- approx(boat_days_jul$Year, boat_days_jul$Boat_Days_Commonwealth, xout=c(2012,2014,2015,2017)) %>% 
  unlist()

# August
boat_days_aug <- boat_days %>% 
  filter(NumMonth==8)

inter_aug_state <- approx(boat_days_aug$Year, boat_days_aug$Boat_Days_State, xout=c(2012,2014,2015,2017)) %>% 
  unlist()
inter_aug_comm <- approx(boat_days_aug$Year, boat_days_aug$Boat_Days_Commonwealth, xout=c(2012,2014,2015,2017)) %>% 
  unlist()

# September
boat_days_sep <- boat_days %>% 
  filter(NumMonth==9)

inter_sep_state <- approx(boat_days_sep$Year, boat_days_sep$Boat_Days_State, xout=c(2012,2014,2016)) %>% 
  unlist()
inter_sep_comm <- approx(boat_days_sep$Year, boat_days_sep$Boat_Days_Commonwealth, xout=c(2012,2014,2016)) %>% 
  unlist()

# October
boat_days_oct <- boat_days %>% 
  filter(NumMonth==10)

inter_oct_state <- approx(boat_days_oct$Year, boat_days_oct$Boat_Days_State, xout=c(2012,2014,2016)) %>% 
  unlist()
inter_oct_comm <- approx(boat_days_oct$Year, boat_days_oct$Boat_Days_Commonwealth, xout=c(2012,2014,2016)) %>% 
  unlist()

# November
boat_days_nov <- boat_days %>% 
  filter(NumMonth==11)

inter_nov_state <- approx(boat_days_nov$Year, boat_days_nov$Boat_Days_State, xout=c(2012,2014,2016)) %>% 
  unlist()
inter_nov_comm <- approx(boat_days_nov$Year, boat_days_nov$Boat_Days_Commonwealth, xout=c(2012,2014,2016)) %>% 
  unlist()

# December
boat_days_dec <- boat_days %>% 
  filter(NumMonth==12)

inter_dec_state <- approx(boat_days_dec$Year, boat_days_dec$Boat_Days_State, xout=c(2012,2014,2016)) %>% 
  unlist()
inter_dec_comm <- approx(boat_days_dec$Year, boat_days_dec$Boat_Days_Commonwealth, xout=c(2012,2014,2016)) %>% 
  unlist()

## Add these values in to the data frame where we are missing the values (oh boy this code is grim)
boat_days[c(25,49,73), 4] <- inter_jan_state[4:6]
boat_days[c(25,49,73), 5] <- inter_jan_comm[4:6]

boat_days[c(26,50,74),4] <- inter_feb_state[4:6]
boat_days[c(26,50,74),5] <- inter_feb_comm[4:6]

boat_days[c(15,27,51,75),4] <- inter_mar_state[5:8]
boat_days[c(15,27,51,75),5] <- inter_mar_comm[5:8]

boat_days[c(16,28,52,76),4] <- inter_apr_state[5:8]
boat_days[c(16,28,52,76),5] <- inter_apr_comm[5:8]

boat_days[c(17,41,53,77),4] <- inter_may_state[5:8]
boat_days[c(17,41,53,77),5] <- inter_may_comm[5:8]

boat_days[c(18,42,54,78),4] <- inter_jun_state[5:8]
boat_days[c(18,42,54,78),5] <- inter_jun_comm[5:8]

boat_days[c(19,43,55,79),4] <- inter_jul_state[5:8]
boat_days[c(19,43,55,79),5] <- inter_jul_comm[5:8]

boat_days[c(20,44,56,80),4] <- inter_aug_state[5:8]
boat_days[c(20,44,56,80),5] <- inter_aug_comm[5:8]

boat_days[c(21,45,69),4] <- inter_sep_state[4:6]
boat_days[c(21,45,69),5] <- inter_sep_comm[4:6]

boat_days[c(22,46,70),4] <- inter_oct_state[4:6]
boat_days[c(22,46,70),5] <- inter_oct_comm[4:6]

boat_days[c(23,47,71),4] <- inter_nov_state[4:6]
boat_days[c(23,47,71),5] <- inter_nov_comm[4:6]

boat_days[c(24,48,72),4] <- inter_dec_state[4:6]
boat_days[c(24,48,72),5] <- inter_dec_comm[4:6]

## Add boat days together to get total boat days for both state and commonwealth waters
boat_days <- boat_days %>% 
  mutate(Total_Boat_Days = Boat_Days_State + Boat_Days_Commonwealth)

## There are few situations in early 2011 and late 2018 where we don't have any data so we're going to create 
## a linear model to try and estimate what those values would be

# (oh god more repetitve stuff)

# 2011 Jan and Feb
January_Boat <- boat_days %>% 
  filter(Month %in% c("Jan")) %>% 
  drop_na()

JanModel <- lm(Total_Boat_Days~Year, data=January_Boat)
summary(JanModel)

Jan_2011 <- as.data.frame(array(0, dim=c(1,2))) %>% 
  mutate(V1=seq(2011, 2011, by=1)) %>% 
  rename("Year"=V1)

predictions <- predict(JanModel, newdata=Jan_2011)
boat_days[1,7] <- predictions

Feb_Boat <- boat_days %>% 
  filter(Month %in% c("Feb")) %>% 
  drop_na()

FebModel <- lm(Total_Boat_Days~Year, data=Feb_Boat)
summary(FebModel)

Feb_2011 <- as.data.frame(array(0, dim=c(1,2))) %>% 
  mutate(V1=seq(2011, 2011, by=1)) %>% 
  rename("Year"=V1)

predictions <- predict(FebModel, newdata=Feb_2011)
boat_days[2,7] <- predictions

# September October November and December of 2018
Sep_Boat <- boat_days %>% 
  filter(Month %in% c("Sep")) %>% 
  drop_na()

SepModel <- lm(Total_Boat_Days~Year, data=Sep_Boat)
summary(SepModel)

Sep_2011 <- as.data.frame(array(0, dim=c(1,2))) %>% 
  mutate(V1=seq(2018, 2018, by=1)) %>% 
  rename("Year"=V1)

predictions <- predict(SepModel, newdata=Sep_2011)
boat_days[93,7] <- predictions

Oct_Boat <- boat_days %>% 
  filter(Month %in% c("Oct")) %>% 
  drop_na()

OctModel <- lm(Total_Boat_Days~Year, data=Oct_Boat)
summary(OctModel)

Oct_2011 <- as.data.frame(array(0, dim=c(1,2))) %>% 
  mutate(V1=seq(2018, 2018, by=1)) %>% 
  rename("Year"=V1)

predictions <- predict(OctModel, newdata=Oct_2011)
boat_days[94,7] <- predictions

Nov_Boat <- boat_days %>% 
  filter(Month %in% c("Nov")) %>% 
  drop_na()

NovModel <- lm(Total_Boat_Days~Year, data=Nov_Boat)
summary(NovModel)

Nov_2011 <- as.data.frame(array(0, dim=c(1,2))) %>% 
  mutate(V1=seq(2018, 2018, by=1)) %>% 
  rename("Year"=V1)

predictions <- predict(NovModel, newdata=Nov_2011)
boat_days[95,7] <- predictions

Dec_Boat <- boat_days %>% 
  filter(Month %in% c("Dec")) %>% 
  drop_na()

DecModel <- lm(Total_Boat_Days~Year, data=Dec_Boat)
summary(DecModel)

Dec_2011 <- as.data.frame(array(0, dim=c(1,2))) %>% 
  mutate(V1=seq(2018, 2018, by=1)) %>% 
  rename("Year"=V1)

predictions <- predict(DecModel, newdata=Dec_2011)
boat_days[96,7] <- predictions

total_plot <- boat_days %>% 
  unite("YearxMonth", c("Year", "NumMonth"), remove = T, sep = "-") %>% 
  mutate(YearxMonth = as.character(YearxMonth)) %>% 
  mutate(YearxMonth = factor(YearxMonth, levels=unique(YearxMonth))) %>% 
  ggplot(.) + 
  #geom_point(aes(x=YearxMonth, y=Total_Boat_Days))+
  geom_line(aes(x=YearxMonth, y=Total_Boat_Days, group=1))


#### HINDCASTING ####
## Going to do the hindcast based on the yearly total of boat days in the whole area
## Has now been changed to max point in 2000 but haven't changed the names of the data frames
TotalYear <- boat_days %>% 
  group_by(Year) %>% 
  summarise(Total=sum(Total_Boat_Days, na.rm = T)) %>% 
  ungroup()

YearModel <- lm(Total~Year, data=TotalYear)
summary(YearModel)

Year2011_1990 <- as.data.frame(array(0, dim=c(11,1))) %>% 
  mutate(V1=seq(2000, 2010, by=1)) %>% 
  rename("Year"=V1)

predictions <- predict(YearModel, newdata=Year2011_1990)

Year2011_1990 <- Year2011_1990 %>% 
  mutate(Total = predictions)

boat_days_hind <- rbind(TotalYear, Year2011_1990)

effort <- seq(0, 47485.01, length=41) # Use the max from the predictive model to then get a straight line back to 0 
years <- seq(1960, 2000, by=1)

Years_1960_1989 <- as.data.frame(cbind(years, effort)) %>% 
  rename("Year" = years) %>% 
  rename("Total" = effort) %>% 
  filter(Year < 2000)

#### JOIN FISHING EFFORT ####
boat_days_hind <- rbind(boat_days_hind, Years_1960_1989)

YearPlot <- ggplot(boat_days_hind) + # Check it looks right
  geom_line(aes(x=Year, y=Total))+ 
  theme_classic()+
  ylab("Total Boat Days")

#### ALLOCATING MONLTHY EFFORT ####

# Work out proportion of fishing effort in each month for each year
boat_days <- boat_days %>% 
  group_by(Year) %>% 
  mutate(Year_Sum = sum(Total_Boat_Days)) %>%
  mutate(Month_Prop = Total_Boat_Days/Year_Sum) %>% 
  dplyr::select(-Year_Sum)

# Work out the average fishing effort for each month
boat_days <- boat_days %>% 
  group_by(Month) %>% 
  mutate(Ave_Month_Prop = mean(Month_Prop))

Month_Prop_Ave <- boat_days[1:12,c(2,9)]

check <- boat_days %>% 
  group_by(Year) %>% 
  summarise(Ave_Month_Prop = sum(Ave_Month_Prop))

check <- boat_days %>% 
  group_by(Year) %>% 
  summarise(Ave_Month = mean(Total_Boat_Days))

# Split up the hindcast data into monthly means 
boat_days_hind <- boat_days_hind %>% 
  bind_rows(replicate(11, boat_days_hind, simplify = FALSE)) %>% # Make a row for each month of each year in the hindcast
  mutate(Year = as.integer(Year)) %>% 
  arrange(-Year) %>% 
  mutate(NumMonth = (rep(1:12, 59))) %>% # Add numerical month
  mutate(Month = rep(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), 59)) %>% #Add word month
  filter(Year<2011) 


boat_days_hind <- boat_days_hind %>% 
  mutate(Total_Boat_Days = 0) %>% 
  left_join(., Month_Prop_Ave) %>% 
  group_by(Year) %>% 
  mutate(Total_Boat_Days = Total*Ave_Month_Prop) %>% 
  ungroup() %>% 
  dplyr::select(-c(Total, Ave_Month_Prop))

check <- boat_days_hind %>% 
  group_by(Year) %>% 
  summarise(., sum(Total_Boat_Days)) # check that the values match

# Put it all together into one big dataframe

Full_Boat_Days <- boat_days[,c(1,2,3,7)] %>% 
  rbind(boat_days_hind) %>% 
  arrange(-Year)

# Plot and check that it looks right
MonthPlot <- Full_Boat_Days %>%
  mutate(Unique = paste(Year, NumMonth, sep="_")) %>%
  filter(Month %in% c("Jan")) %>%
  ggplot() +
  geom_point(aes(x=Unique, y=Total_Boat_Days))

#### ALLOCATION EFFORT TO EXMOUTH AND BOAT RAMPS ####

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
  if(Full_Boat_Days[Y,1]<2008){ # This is for when Bundegi wasn't a proper ramp so probably would have had less effort
    Full_Boat_Days$Tb_BR = Full_Boat_Days$Total_Boat_Days*BR_Trips[1,6]
    Full_Boat_Days$Bd_BR = Full_Boat_Days$Total_Boat_Days*BR_Trips[2,6]
    Full_Boat_Days$ExM_BR = Full_Boat_Days$Total_Boat_Days*BR_Trips[3,6] 
    Full_Boat_Days$CrB_BR = Full_Boat_Days$Total_Boat_Days*BR_Trips[4,6]
  } else {
    Full_Boat_Days$Tb_BR = Full_Boat_Days$Total_Boat_Days*BR_Trips[1,5] 
    Full_Boat_Days$Bd_BR = Full_Boat_Days$Total_Boat_Days*BR_Trips[2,5] 
    Full_Boat_Days$ExM_BR = Full_Boat_Days$Total_Boat_Days*BR_Trips[3,5] 
    Full_Boat_Days$CrB_BR = Full_Boat_Days$Total_Boat_Days*BR_Trips[4,5]
  }
}

check <- Full_Boat_Days %>% 
  mutate(Total = Tb_BR+Bd_BR+ExM_BR+CrB_BR)
#### ALLOCATING EFFORT TO CELLS ####

## Work out the probability of visiting a cell from each boat ramp based on distance and size
BR <- st_sf(BR)

water <- water %>% 
  mutate(DistBR = 0) %>% 
  mutate(cell_area = st_area(Spatial))

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

BR_Trips <- Full_Boat_Days %>% # This is just the trips from each boat ramp
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

#### Determining catchability ####
## Have estimated an increase in recreational efficiency of 30% from 1990 when colour sounders become more commonplace (GPS was 2002) according to Marriott et al. 2011
q <- as.data.frame(array(0, dim=c(59,1))) %>% 
  rename(Q = "V1")

q[1:30, 1] <- 0.005

# q has to now increase by 30% in equal steps over 29 years 

for (y in 31:59){
  q[y,1] <- q[y-1,1] + 0.00022413793
}

Fishing2 <- Fishing

# Now calculate F by multiplying our effort by q
for (y in 1:59){
  Fishing[,,y] <- Fishing[,,y]*q[y,1]
}

# Remove fishing effort for six months out of the year to simulate a temporal closure from 2008 which covers the spawning season
Fishing[ ,c(1,2,3,10,11,12), 28:59] <- 0

#### SAVE DATA ####
setwd(sim_dir)
saveRDS(Fishing, file=paste0(model.name, sep="_", "S03_fishing"))



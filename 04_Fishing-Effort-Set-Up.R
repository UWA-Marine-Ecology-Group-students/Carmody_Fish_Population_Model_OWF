###################################################

# Script for setting up the fishing surface 
# Requires files that are made in the first script
# Also requires fishing effort data

###################################################


## NEED TO CHANGE THE INITIAL FISHING MORTALITY AT 1960 TO BE THE SAME AS THE EQLUILIBRIUM LEVEL ##
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
library(sfnetworks)

rm(list = ls())

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
setwd(sg_dir)

water <- readRDS(paste0(model.name, sep="_","water"))
NCELL <- nrow(water)

setwd(sp_dir)
network <- st_read(paste0(model.name, sep="_","network.shapefile.shp"))

# Locations of the boat ramps
setwd(sp_dir)
BR <- st_read("Boat_Ramps.shp") %>% 
  st_transform(4283)%>%
  st_make_valid() 
BR <- BR[1:4,]

# No Take Zones
setwd(sg_dir)
NoTake <- readRDS(paste0(model.name, sep="_","NoTakeList"))

#### FILL IN MISSING VALUES FOR NINGALOO BOAT DAYS ####

# Work out the number of boat days for the Ningaloo region split by state and commonwealth
# We have estimates for just this area from Claire Smallwood but not for all years 
# We have calculated proportions for each month based on the Gascoyne data 
# So we are going to interpolate the months we are missing 

## Interpolate the values for state and commonwealth separately
state_values <- list()
comm_values <- list()

for(M in 1:12){
  if(M<3){
    boat_days_month <- boat_days %>% 
      filter(NumMonth==M)
    
    state_values[[M]] <- approx(boat_days_month$Year, boat_days_month$Boat_Days_State, xout=c(2013,2015,2017)) %>% 
      unlist()
    comm_values[[M]] <- approx(boat_days_month$Year, boat_days_month$Boat_Days_Commonwealth, xout=c(2013,2015,2017)) %>% 
      unlist()
  } else if (M==3|M==4){
    boat_days_month <- boat_days %>% 
      filter(NumMonth==M)
    
    state_values[[M]] <- approx(boat_days_month$Year, boat_days_month$Boat_Days_State, xout=c(2012,2013,2015,2017)) %>% 
      unlist()
    comm_values[[M]] <- approx(boat_days_month$Year, boat_days_month$Boat_Days_Commonwealth, xout=c(2012,2013,2015,2017)) %>% 
      unlist()
  } else if(M>4 & M<9) {
    boat_days_month <- boat_days %>% 
      filter(NumMonth==M)
    
    state_values[[M]] <- approx(boat_days_month$Year, boat_days_month$Boat_Days_State, xout=c(2012,2014,2015,2017)) %>% 
      unlist()
    comm_values[[M]] <- approx(boat_days_month$Year, boat_days_month$Boat_Days_Commonwealth, xout=c(2012,2014,2015,2017)) %>% 
      unlist()
  } else if(M>=9){
    boat_days_month <- boat_days %>% 
      filter(NumMonth==M)
    state_values[[M]] <- approx(boat_days_month$Year, boat_days_month$Boat_Days_State, xout=c(2012,2014,2016)) %>% 
      unlist()
    comm_values[[M]] <- approx(boat_days_month$Year, boat_days_month$Boat_Days_Commonwealth, xout=c(2012,2014,2016)) %>% 
      unlist()
    }
}


## Add these values in to the data frame where we are missing the values
shift <- c(0,1,0,1,0,1,2,3,0,1,2,3) # To account for the fact that different months and years are in different places in the big data frame

for(M in 1:12){
  if(M<3){
    
    boat_days[c(25+shift[M],49+shift[M],73+shift[M]), 4] <- state_values[[M]][4:6]
    boat_days[c(25+shift[M],49+shift[M],73+shift[M]), 5] <- comm_values[[M]][4:6]

  } else if (M==3|M==4){
    boat_days[c(15+shift[M],27+shift[M],51+shift[M],75+shift[M]),4] <- state_values[[M]][5:8]
    boat_days[c(15+shift[M],27+shift[M],51+shift[M],75+shift[M]),5] <- comm_values[[M]][5:8]
    
  } else if(M>4 & M<9) {

    boat_days[c(17+shift[M],41+shift[M],53+shift[M],77+shift[M]),4] <- state_values[[M]][5:8]
    boat_days[c(17+shift[M],41+shift[M],53+shift[M],77+shift[M]),5] <- comm_values[[M]][5:8]
    
  } else if(M>=9){
    boat_days[c(21+shift[M],45+shift[M],69+shift[M]),4] <- state_values[[M]][4:6]
    boat_days[c(21+shift[M],45+shift[M],69+shift[M]),5] <- comm_values[[M]][4:6]
  }
}

## Add boat days together to get total boat days for both state and commonwealth waters
boat_days <- boat_days %>% 
  mutate(Total_Boat_Days = Boat_Days_State + Boat_Days_Commonwealth)

## There are few situations in early 2011 and late 2018 where we don't have any data so we're going to create 
## a linear model to try and estimate what those values would be
Months <- c("Jan", "Feb", "Sep" ,"Oct", "Nov", "Dec")
rows <- c(1,2,93,94,95,96)

for(M in 1:length(Months)){
  if(M<3){
    Boat <- boat_days %>% 
      filter(Month %in% c(Months[M])) %>% 
      drop_na()
    
    MonthModel <- lm(Total_Boat_Days~Year, data=Boat)
    
    temp <- as.data.frame(array(0, dim=c(1,2))) %>% 
      mutate(V1=seq(2011, 2011, by=1)) %>% 
      rename("Year"=V1)
    
    predictions <- predict(MonthModel, newdata = temp)
    
    boat_days[rows[M], 7] <- predictions
    
  } else if(M>=3){
    Boat <- boat_days %>% 
      filter(Month %in% c(Months[M])) %>% 
      drop_na()
    
    MonthModel <- lm(Total_Boat_Days~Year, data=Boat)
    
    temp <- as.data.frame(array(0, dim=c(1,2))) %>% 
      mutate(V1=seq(2018, 2018, by=1)) %>% 
      rename("Year"=V1)
    
    predictions <- predict(MonthModel, newdata = temp)
    
    boat_days[rows[M], 7] <- predictions
  }
}



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

EffortPlot <- Full_Boat_Days %>%
  group_by(Year) %>% 
  summarise(., sum(Total_Boat_Days)) %>% 
  ggplot() +
  geom_line(aes(x=Year, y=`sum(Total_Boat_Days)`)) + 
  theme_classic()+
  ylab("Effort (Boat Days)")+
  geom_vline(xintercept=1991, linetype="dotted", color="#302383")+
  geom_vline(xintercept=1995, linetype="dashed", colour="#66CCEE")
EffortPlot

setwd(fig_dir)
a4.width <- 160
ggsave(EffortPlot, filename="Effort_Plot.png", height = a4.width*1, width = a4.width, units  ="mm", dpi = 300 )

temp <- Full_Boat_Days %>% 
  group_by(Year) %>% 
  summarise(., sum(Total_Boat_Days))

#### CALCULATE CATCHABILITY VARYING BY CELL SIZE ####
water <- water %>% 
  mutate(Area = as.vector((water$cell_area)/1000000))

water_area <- water %>% 
  dplyr::select(Fished_1960, Fished_1987, Fished_2005, Fished_2017, Area) %>% 
  mutate(Area_60 = ifelse(Fished_1960 %in% c("N"), 0, Area)) %>% 
  mutate(Area_87 = ifelse(Fished_1987 %in% c("N"), 0, Area)) %>% 
  mutate(Area_05 = ifelse(Fished_2005 %in% c("N"), 0, Area)) %>% 
  mutate(Area_17 = ifelse(Fished_2017 %in% c("N"), 0, Area)) %>% 
  dplyr::select(Area_60, Area_87, Area_05, Area_17) %>% 
  mutate(Sum_60 = sum(Area_60)) %>% 
  mutate(Sum_87 = sum(Area_87)) %>% 
  mutate(Sum_05 = sum(Area_05)) %>% 
  mutate(Sum_17 = sum(Area_17)) %>% 
  st_drop_geometry() %>% 
  mutate(q_60 = Area_60/Sum_60) %>% 
  mutate(q_87 = Area_87/Sum_87) %>% 
  mutate(q_05 = Area_05/Sum_05) %>% 
  mutate(q_17 = Area_17/Sum_17) %>% 
  mutate(ID = row_number())

spatial_q <- array(0.000005, dim=c(NCELL, 59))

for (y in 31:51){
  spatial_q[ ,y] <- spatial_q[ ,y-1] * 1.02
}

spatial_q[,52:59] <- spatial_q[,51]

for(COL in 1:27){
  for (ROW in 1:NCELL){
    
    spatial_q[ROW,COL] <- 0.3*(spatial_q[ROW,COL]/water_area[ROW,9])
    
  }
}

for(COL in 28:45){
  for (ROW in 1:NCELL){
    
    spatial_q[ROW,COL] <- 0.3*(spatial_q[ROW,COL]/water_area[ROW,10])
    
  }
}

for(COL in 46:57){
  for (ROW in 1:NCELL){
    
    spatial_q[ROW,COL] <- 0.3*(spatial_q[ROW,COL]/water_area[ROW,11])
    
  }
}

for(COL in 58:59){
  for (ROW in 1:NCELL){
    
    spatial_q[ROW,COL] <- 0.3*(spatial_q[ROW,COL]/water_area[ROW,12])
    
  }
}

spatial_q[spatial_q == Inf] <- 0

max(spatial_q)
## Save the monthly allocations for use in later things 
setwd(sg_dir)

saveRDS(Month_Prop_Ave, file="Average_Monthly_Effort")

#### FOR SENSITIVITY ANALYSIS - REMEMBER TO TURN OFF NINGALOO SAVE ####
## Changing F 

# Full_Boat_Days <- Full_Boat_Days %>%
#   mutate(Total_Boat_Days = Total_Boat_Days*2)

#### ALLOCATION EFFORT TO BOAT RAMPS ####

## NEED TO CHANGE TOTAL BOAT DAYS TO FISHING MORTALITY

# Full_Boat_Days <- Full_Boat_Days %>% 
#   mutate(Total_Boat_Days = Total_Boat_Days*(0.005/12))

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
  mutate(BR_Prop_08 = c(0.3031544, 0.1577101, 0.3119251, 0.2272104)) # Create separate proportions for Bundegi before 2008 as the boat ramp pretty much didn't exist, have allocated 10% of its boats to the other Exmouth Boat Ramps


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

#### CREATE SEPARATE VECTORS OF ROW IDS FROM THE NO TAKE LIST ####

# Get the cell IDs/rows for the no take cells in each portion of the model

NoTake87_05 <- NoTake[[1]] 

NoTake05_18 <- NoTake[[2]] 

NoTake18_21 <- NoTake[[3]] 


#### DISTANCE FROM EACH BOAT RAMP TO EACH CELL ####

## Work out the probability of visiting a cell from each boat ramp based on distance and size
BR <- st_as_sf(BR)
st_crs(BR) <- NA 

BR <- BR %>% 
  mutate(name = c("Bundegi","Exmouth","Tantabiddi","CoralBay"))

centroids <- st_centroid_within_poly(water)
points <- as.data.frame(st_coordinates(centroids))%>% #The points start at the bottom left and then work their way their way right
  mutate(ID=row_number()) 
points_sf <- st_as_sf(points, coords = c("X", "Y")) 
st_crs(points_sf) <- NA


network <- as_sfnetwork(network, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

net <- activate(network, "nodes")
network_matrix <- st_network_cost(net, from=BR, to=points_sf)
network_matrix <- network_matrix*111
dim(network_matrix)

DistBR <- as.data.frame(t(network_matrix)) %>% 
  rename("Bd_BR"=V1) %>% 
  rename("ExM_BR" = V2) %>% 
  rename("Tb_BR" = V3) %>% 
  rename("CrB_BR"=V4)


#### SETTING UP UTILITY FUNCTION ####
# Create a data frame with both the distances and the areas of the cells
Cell_Vars <- DistBR %>% 
  mutate(Area = as.vector((water$cell_area)/1000000))#Cells are now in km^2 but with no units

## Now need to create a separate fishing surface for each month of each year based on distance to boat ramp, size of each
## cell and multiply that by the effort in the cell to spatially allocate the effort across the area
## But we need to account for the fact that there will be sanctuary zones going in and the effort that would have gone in there will get allocated somewhere else
## Will then need to put the rows/columns back in as 0s 
NCELL_6086 <- NCELL
NCELL_8705 <- NCELL-length(NoTake[[1]])
NCELL_0517 <- NCELL-length(NoTake[[2]])
NCELL_1718 <- NCELL-length(NoTake[[3]])

# Add coefficients to each variable - all the BRs are the same but make them negative because the further away they are the less likely people are to visit
a = 0.0005
b = -0.001
c = -0.001
d = -0.001
e = -0.001

Vj <- Cell_Vars %>% 
  mutate(Bd_BR = Bd_BR*b,
         ExM_BR = ExM_BR*c,
         Tb_BR = Tb_BR*d,
         CrB_BR = CrB_BR*e,
         Area= Area*a) %>% 
  mutate(vj = Bd_BR+ExM_BR+Tb_BR+CrB_BR)

Vj_6086 <- Vj
Vj_8705 <- Vj[-c(NoTake[[1]]), ]
Vj_0517 <- Vj[-c(NoTake[[2]]), ]
Vj_1718 <- Vj[-c(NoTake[[3]]), ]



## 1960-1987 
BR_U_6086 <- as.data.frame(matrix(0, nrow=NCELL_6086, ncol=4)) #Set up data frame to hold utilities of cells

BR_U_6086 <- BR_U_6086 %>% #make sure you give the columns good names so that you know what they are
  rename("Bd_U"=V1) %>% 
  rename("ExM_U" = V2) %>% 
  rename("Tb_U" = V3) %>% 
  rename("CrB_U"=V4)

cellU <- matrix(NA, ncol=4, nrow=NCELL_6086)

for(RAMP in 1:4){
  for(cell in 1:NCELL_6086){
    U <- exp(Vj_6086[cell,5]/Vj_6086[cell,RAMP])
    cellU[cell, RAMP] <- U
  }
}

rowU <- as.data.frame(colSums(cellU))

for (RAMP in 1:4){
  for (cell in 1:NCELL_6086){
    BR_U_6086[cell,RAMP] <- (exp(Vj_6086[cell,5]/Vj_6086[cell,RAMP]))/rowU[RAMP,1]
  }
}
colSums(BR_U_6086)


## 1987-2005

BR_U_8705 <- as.data.frame(matrix(0, nrow=NCELL_8705, ncol=4)) #Set up data frame to hold utilities of cells

BR_U_8705 <- BR_U_8705%>% #make sure you give the columns good names so that you know what they are
  rename("Bd_U"=V1) %>% 
  rename("ExM_U" = V2) %>% 
  rename("Tb_U" = V3) %>% 
  rename("CrB_U"=V4)

cellU <- matrix(NA, ncol=4, nrow=NCELL_8705)

for(RAMP in 1:4){
  for(cell in 1:NCELL_8705){
    U <- exp(Vj_8705[cell,5]/Vj_8705[cell,RAMP])
    cellU[cell, RAMP] <- U
  }
}

rowU <- as.data.frame(colSums(cellU))

for (RAMP in 1:4){
  for (cell in 1:NCELL_8705){
    BR_U_8705[cell,RAMP] <- (exp(Vj_8705[cell,5]/Vj_8705[cell,RAMP]))/rowU[RAMP,1]
  }
}
colSums(BR_U_8705)

## 2005-2017

BR_U_0517 <- as.data.frame(matrix(0, nrow=NCELL_0517, ncol=4)) #Set up data frame to hold utilities of cells

BR_U_0517 <- BR_U_0517%>% #make sure you give the columns good names so that you know what they are
  rename("Bd_U"=V1) %>% 
  rename("ExM_U" = V2) %>% 
  rename("Tb_U" = V3) %>% 
  rename("CrB_U"=V4)

cellU <- matrix(NA, ncol=4, nrow=NCELL_0517)

for(RAMP in 1:4){
  for(cell in 1:NCELL_0517){
    U <- exp(Vj_0517[cell,5]/Vj_0517[cell,RAMP])
    cellU[cell, RAMP] <- U
  }
}

rowU <- as.data.frame(colSums(cellU))

for (RAMP in 1:4){
  for (cell in 1:NCELL_0517){
    BR_U_0517[cell,RAMP] <- (exp(Vj_0517[cell,5]/Vj_0517[cell,RAMP]))/rowU[RAMP,1]
  }
}
colSums(BR_U_0517)

## 2017-2018

BR_U_1718 <- as.data.frame(matrix(0, nrow=NCELL_1718, ncol=4)) #Set up data frame to hold utilities of cells

BR_U_1718 <- BR_U_1718 %>% #make sure you give the columns good names so that you know what they are
  rename("Bd_U"=V1) %>% 
  rename("ExM_U" = V2) %>% 
  rename("Tb_U" = V3) %>% 
  rename("CrB_U"=V4)

cellU <- matrix(NA, ncol=4, nrow=NCELL_1718)

for(RAMP in 1:4){
  for(cell in 1:NCELL_1718){
    U <- exp(Vj_1718[cell,5]/Vj_1718[cell,RAMP])
    cellU[cell, RAMP] <- U
  }
}

rowU <- as.data.frame(colSums(cellU))

for (RAMP in 1:4){
  for (cell in 1:NCELL_1718){
    BR_U_1718[cell,RAMP] <- (exp(Vj_1718[cell,5]/Vj_1718[cell,RAMP]))/rowU[RAMP,1]
  }
}
colSums(BR_U_1718)

#### ALLOCATING EFFORT TO CELLS ####

BR_Trips <- Full_Boat_Days %>% # This is just the trips from each boat ramp
  ungroup() %>% 
  mutate(NumYear = rep(59:1, each=12)) %>% #This is to turn the years into a count for the loop
  dplyr::select(NumYear, NumMonth, Bd_BR, ExM_BR, Tb_BR, CrB_BR)  

# 1960-1986
Fishing_6086 <- array(0, dim=c(NCELL_6086, 12, 27)) #This array has a row for every cell, a column for every month, and a layer for every year
Months <- array(0, dim=c(NCELL, 12))
Ramps <- array(0, dim=c(NCELL, 4))
layer <- 1

for(YEAR in 1:27){
  
  for(MONTH in 1:12){
    
    for(RAMP in 1:4){
      
      temp <- BR_Trips %>% 
        filter(NumYear==YEAR) %>% 
        dplyr::select(-c(NumYear, NumMonth))
      
      temp <- as.matrix(temp) 
      
      for(CELL in 1:NCELL_6086){
        Ramps[CELL,RAMP] <- BR_U_6086[CELL,RAMP] * temp[MONTH,RAMP]
      }
    }
    Months[,MONTH] <- rowSums(Ramps)
  }
  Fishing_6086[ , ,layer] <- Months 
  layer <- layer+1
}

# 1987-2004
Fishing_8705 <- array(0, dim=c(NCELL_8705, 12,18)) #This array has a row for every cell, a column for every month, and a layer for every year
Months <- array(0, dim=c(NCELL_8705, 12))
Ramps <- array(0, dim=c(NCELL_8705, 4))
layer <- 1

for(YEAR in 28:45){
  
  for(MONTH in 1:12){
    
    for(RAMP in 1:4){
      
      temp <- BR_Trips %>% 
        filter(NumYear==YEAR) %>% 
        dplyr::select(-c(NumYear, NumMonth))
      
      temp <- as.matrix(temp) 
      
      for(CELL in 1:NCELL_8705){
        Ramps[CELL,RAMP] <- BR_U_8705[CELL,RAMP] * temp[MONTH,RAMP]
      }
    }
    
    Months[,MONTH] <-  rowSums(Ramps)
  }
  Fishing_8705[ , ,layer] <- Months 
  layer <- layer+1
}

# 2005-2017
Fishing_0517 <- array(0, dim=c(NCELL_0517, 12,12)) #This array has a row for every cell, a column for every month, and a layer for every year
Months <- array(0, dim=c(NCELL_0517, 12))
Ramps <- array(0, dim=c(NCELL_0517, 4))
layer <- 1

for(YEAR in 46:57){
  
  for(MONTH in 1:12){
    
    for(RAMP in 1:4){
      
      temp <- BR_Trips %>% 
        filter(NumYear==YEAR) %>% 
        dplyr::select(-c(NumYear, NumMonth))
      
      temp <- as.matrix(temp) 
      
      for(CELL in 1:NCELL_0517){
        Ramps[CELL,RAMP] <- BR_U_0517[CELL,RAMP] * temp[MONTH,RAMP]
      }
    }
    
    Months[,MONTH] <-  rowSums(Ramps)
  }
  Fishing_0517[ , ,layer] <- Months 
  layer <- layer+1
}

# 2017-2018
Fishing_1718 <- array(0, dim=c(NCELL_1718, 12, 2)) #This array has a row for every cell, a column for every month, and a layer for every year
Months <- array(0, dim=c(NCELL_1718, 12))
Ramps <- array(0, dim=c(NCELL_1718, 4))
layer <- 1

for(YEAR in 58:59){
  
  for(MONTH in 1:12){
    
    for(RAMP in 1:4){
      
      temp <- BR_Trips %>% 
        filter(NumYear==YEAR) %>% 
        dplyr::select(-c(NumYear, NumMonth))
      
      temp <- as.matrix(temp) 
      
      for(CELL in 1:NCELL_1718){
        Ramps[CELL,RAMP] <- BR_U_1718[CELL,RAMP] * temp[MONTH,RAMP]
      }
    }
    
    Months[,MONTH] <-  rowSums(Ramps)
  }
  Fishing_1718[ , ,layer] <- Months 
  layer <- layer+1
}

## Add cells in the NTZs back in 
# Have to do this by adding in rows in to the matrix in the correct places but just give them 0

# 1987-2004
NTCells <- NoTake[[1]]
Fishing_8705_2 <- array(0, dim=c(NCELL,12,18))

Fishing_8705_2[-NTCells,,] <- Fishing_8705

# 2005-2016
NTCells <- NoTake[[2]]
Fishing_0517_2 <- array(0, dim=c(NCELL,12,12))

Fishing_0517_2[-NTCells,,] <- Fishing_0517

# 2017-2018
NTCells <- NoTake[[3]]
Fishing_1718_2 <- array(0, dim=c(NCELL,12,2))

Fishing_1718_2[-NTCells,,] <- Fishing_1718

## Put it all back together 
Fishing <- abind(Fishing_6086, Fishing_8705_2, along=3)
Fishing <- abind(Fishing, Fishing_0517_2, along=3)
Fishing <- abind(Fishing, Fishing_1718_2, along=3)

### ADD CATCHABILITY ####
for (YEAR in 1:59){
  Fishing[,,YEAR] <- Fishing[,,YEAR] * spatial_q[,YEAR]
}

#### SAVE DATA ####
setwd(sg_dir)

saveRDS(Fishing, file=paste0("ningaloo", sep="_", "fishing"))
# saveRDS(NoTake, file="NoTake")
# saveRDS(NoTake, file="NoTakeList")

#### SETTING UP EFFOR FOR BURN IN ####
## Need to set effort to be a consistent low level for the burn in, we can use the same value we're using to set up the initial population 

BR_Trips <- data.frame(BoatRamp = c("Tantabiddi", "Bundegi", "ExmouthMar", "CoralBay"),
                       Effort = c(194.48, 133.9, 166.15, 146.07),
                       Trips = c(224, 157, 198, 151))

BR_Trips <- BR_Trips %>% 
  mutate(Trip_per_Hr = as.numeric(unlist((Trips/Effort)))) %>% # Standardise the no. trips based on how much time you spent sampling
  mutate(BR_Prop = Trip_per_Hr/sum(Trip_per_Hr)) %>% #Then work out the proportion of trips each hour that leave from each boat ramp
  mutate(BR_Prop_08 = c(0.3031544, 0.1577101, 0.3119251, 0.2272104)) # Create separate proportions for Bundegi before 2008 as the boat ramp pretty much didn't exist, have allocated 10% of its boats to the other Exmouth Boat Ramps

# Fishing parameters
eq.init.fish = 0.025 
q = 0.0005
Effort = (-log(1-eq.init.fish))/q # We assume the same level of nominal effort in each year

## Split up this effort by the same proprtions as before and allocate it to the different access points
# Months
burn_in_effort <- Month_Prop_Ave[,2]*Effort

burn_in_effort <- as.data.frame(burn_in_effort) %>% 
  mutate(Tb_BR = 0,
         Bd_BR = 0,
         ExM_BR = 0,
         CrB_BR = 0) %>% 
  rename(Effort = "Ave_Month_Prop")


## Split up by boat ramp
for(M in 1:12){
  burn_in_effort$Tb_BR = burn_in_effort$Effort*BR_Trips[1,6]
  burn_in_effort$Bd_BR = burn_in_effort$Effort*BR_Trips[2,6]
  burn_in_effort$ExM_BR = burn_in_effort$Effort*BR_Trips[3,6] 
  burn_in_effort$CrB_BR = burn_in_effort$Effort*BR_Trips[4,6]
}

## Allocate to the cells using the same utilities that we set up earlier
Burn_In_Fishing <- array(0, dim=c(NCELL, 12, 75)) #This array has a row for every cell, a column for every month, and a layer for every year
Months <- array(0, dim=c(NCELL, 12))
Ramps <- array(0, dim=c(NCELL, 4))
layer <- 1

for(YEAR in 1:50){
  
  for(MONTH in 1:12){
    
    for(RAMP in 1:4){
      
      temp <- burn_in_effort %>% 
        dplyr::select(-c(Effort))
      
      temp <- as.matrix(temp) 
      
      for(CELL in 1:NCELL){
        Ramps[CELL,RAMP] <- BR_U_6086[CELL,RAMP] * temp[MONTH,RAMP] # Use the same utility from before as nothing should have changed
      }
    }
    
    Months[,MONTH] <-  rowSums(Ramps)
  }
  Burn_In_Fishing[ , ,layer] <- Months 
  Burn_In_Fishing[ , ,layer] <- Burn_In_Fishing[ , ,layer] * spatial_q[,1]
  layer <- layer+1
}

setwd(sg_dir)
saveRDS(Burn_In_Fishing, file=paste0("ningaloo", sep="_", "burn_in_fishing"))

#### FOR SMALL MODEL ####
## Small model burn in
model.name <- "small"

setwd(sp_dir)
small_water <- readRDS(paste0(model.name, sep="_", "water"))
water <- readRDS("water")

setwd(sg_dir)
burn_in <- readRDS("ningaloo_burn_in_fishing")

small_water_id <- st_intersects(small_water, water) %>%
  as.data.frame(.) %>%
  distinct(row.id, .keep_all = T)

small_burn_in <- burn_in[as.numeric(small_water_id$col.id),,]

setwd(sg_dir)
saveRDS(small_burn_in, file=paste0(model.name, sep="_", "burn_in_fishing"))
## Don't want to redo fishing effort because then we'd have an insane amount of fishing effort in this small area, just want to take the appropriate cells and keep the fishing effort as if it was the larger model

setwd(sg_dir)
Fishing <- readRDS("ningaloo_fishing")

small_fishing <- Fishing[as.numeric(small_water_id$col.id),,]
small_fishing <- small_fishing/10

setwd(sg_dir)
saveRDS(small_fishing, file=paste0(model.name, sep="_", "fishing"))





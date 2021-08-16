## ---------------------------
##
## Script name: National RUM
##
## Purpose of script: To develop a national RUM of recreational fishing
##
## Author: Dr. Matt Navarro
##
## ----------------------------
##
## Notes:
##  -Site choice data is taken from 2 surveys: NSW/ACT state-wide phone diary survey, and WA i survey data
##  -Models should be fit for: Boat based line fishing in marine waters.
##  -All covariates need to be Nationally avilable to allow extrapoliation
##  
##  Covariates:
##  -Depth and rugosity: GA 2009 Bathymetry and Topograpgy grid https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/67703
##  -sst taken from IMOS (monthly averages)
##  -wind/swell taken from CSIRO hindcasts: https://data.csiro.au/collections/#collection/CIcsiro:39819
##
##  Projections:
##  4283 is GDA94 (Main unprojected crs)
##  3112 is Australian Lambert (Main projected csr)
##  4326 is WGS84 (Used in sst raster extraction)
## ---------------------------

NECTAR = 1

## LOAD PACKAGES----

library(raster)
library(rgeos)
library(rgdal)
library(tidyverse)
library(broom)
library(devtools)
library(readxl)
library(sf)
library(sp)
library(magrittr)
library(ncdf4)
library(lubridate)
library(beepr)
library(reshape2)
library(googledrive)
library(foreach)
library(data.table)
library(lubridate)
library(tidync)
library(spatialEco)
library(units)
library(doSNOW)
require(doParallel)
require(parallel)
require(tcltk)
library(MASS)
library(matrixStats)
library(gridExtra)
library(grid.text)
library(gtable)
library(cowplot)


## LOAD FUNCTIONS----
devtools::source_url("https://github.com/UWA-Marine-Ecology-Group-projects/Useful-spatial-functions/blob/master/Useful-spatial-functions.R?raw=TRUE")
devtools::source_url("https://github.com/aodn/imos-user-code-library/blob/master/R/commons/NetCDF/ncParse.R?raw=TRUE")

outerBuffer<-function(x, dist){ #Calculates outer buffer of a polygon
  buff<-buffer(x, width = dist - 1, dissolve = T)
  e<-erase(buff,x)
  return(e)
}

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

#DEFAULT FIGURE SIZES
ax.title<-14
ax.text<-9
ax.legend<-12
strip<-9
a4_width<- 160
my_theme<-theme_classic()+ theme( axis.text=element_text(size=ax.text),
                                  text = element_text(size=ax.text),
                                  axis.title=element_text(size=ax.title),
                                  line = element_line(size = 0.3),
                                  strip.background = element_rect(colour = "white", fill = "gray80"),
                                  strip.text = element_text(size=strip))

colours<-c("#EDBF42", "#DF3C34","#A05FB3", "#01A8DE" , "#495BC9", "#01AC8C")

## LOAD DATA ONTO NECTAR----Only needs doing once

#Load data onto NECTAR (Only need to do this once)
#AusMap
# drive_download(as_id("1B-RGZnuQcpS1IXQqLXGQk8RJyYg9q54O"), path = "~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/data.zip", 
#   overwrite = TRUE)
# unzip("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/data.zip")
# 
# #wa data
# drive_download(as_id("1xFZvQBBK_B4NviLgTaMnVMVQlo_SCIiV"), path = "~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/wa_data.csv", 
#   overwrite = TRUE)
# drive.download(as_id("1LrwbWGqL3f-yqNVbDZO07KRkInjKk6Rn"))


## LOAD DATA ------
if(NECTAR == 1){#Nectar data loading
  #wa_trips<-read.csv("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/wa.trips.csv")
  
  #nsw_trips<- read.csv( "~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/nsw_trip_expanded.csv")
  #nsw_trip<- read_excel("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/NSWACT_2013-14_Data for RUM_2.xlsx", sheet = 1)
  
  Australia_map <-st_read("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/AusOutline_HighRes.shp")%>%
    st_transform(4283) %>% filter(!is.na(st_dimension(.))) 
  Australia_map_Lambert <- Australia_map %>% st_transform(3112)
  # 
  # wa_blocks <- st_read("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/Grid_Blocks_10nm.shp") %>%
  #   st_transform(4283) %>% filter(!is.na(st_dimension(.))) #change to GDA94 and filter out empty geometries
  # 
  # NSW_sites <- st_read("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/SWS-13_14_FishingEvent_Locations.shp") %>%
  #   st_transform(4283) %>% filter(!is.na(st_dimension(.)))
  
  #National_grid_Lambert <- st_read("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/National_grid_trimmed_Lambert.shp") %>%  #use this grid because it is the basis of the boat ramp grids
  #   filter(!is.na(st_dimension(.)))
  
  #This is the newer data with 31646 objects
  National_grid_Lambert2 <- st_read("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/National_grid_Zones_Lambert2.shp") %>%  #use this grid because it is the basis of the boat ramp grids
    filter(!is.na(st_dimension(.)))
  
  #Ramps <- st_read("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/Jaq_Ramps.shp") %>%  
  #  st_transform(4283) %>% filter(!is.na(st_dimension(.)))
  
  Ramps <- st_read("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/NationalBoatRampsReduced5.shp") %>%
    st_transform(3112) %>% dplyr::filter(!is.na(st_dimension(.)))
  #Quick bit of cleaning on ramps
  Ramps%<>%dplyr::rename("Ramp_name" = "Ramp_nm", "RampRegion_name"="RmpRgn_" ) 
  
  
  
  Bathy <- raster("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/ausbath_09_v4/w001001.adf")
    <- raster("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/VRM.tif")
  
}

if(NECTAR == 0){#laptop data loading
  wa_trips<-read.csv(   "C:/GitHub/Poject-NESPD6-MBH-NationalRUM/Data/wa.trips.csv")
  wa_trips_national<-read.csv("C:/GitHub/Poject-NESPD6-MBH-NationalRUM/Data/wa_trips_national.csv")
  nsw_trip<- read_excel("C:/GitHub/Poject-NESPD6-MBH-NationalRUM/Data/NSWACT_2013-14_Data for RUM_2.xlsx", sheet = 1)
  Australia_map <-st_read("C:/GitHub/Poject-NESPD6-MBH-NationalRUM/Data/Shapefiles/AusOutline_HighRes.shp")%>%
    st_transform(4283) %>% filter(!is.na(st_dimension(.))) 
  Australia_map_Lambert <- Australia_map %>% st_transform(3112)
  National_grid <- st_read("C:/GitHub/Poject-NESPD6-MBH-NationalRUM/Data/Shapefiles/TheOneGrid_ver2.shp") %>%  #use this grid because it is the basis of the boat ramp grids
    st_transform(4283) %>% filter(!is.na(st_dimension(.)))%>%
    #Just for toy version
    st_crop(xmin = 150.4078, xmax = 153.561, ymin = -35.73653, ymax = -29.10684)
}

## 1.0 NATIONAL GRID (n=31646)------

#1.1 GET DEPTH AND VRM
National_grid_wgs84<-st_transform(National_grid_Lambert2, 4326)

# #depth
National_grid_wgs84<-extract_mean_from_raster(Bathy, National_grid_wgs84)
National_grid_wgs84$depth<-National_grid_wgs84$mean

# # #vrm
# National_grid_wgs84<-extract_mean_from_raster(vrm, National_grid_wgs84)
# National_grid_wgs84$vrm<-National_grid_wgs84$mean

any(is.na(National_grid_wgs84$depth))
#any(is.na(National_grid_wgs84$vrm))
nrow(National_grid_wgs84)==31646

#1.2 GET ISLANDS
#produce island data sets
National_grid_Lam<-National_grid_wgs84 %>% st_transform(3112)

Australia_map_Lambert$area <- st_area(Australia_map_Lambert)
Island_small<-Australia_map_Lambert %>% filter(as.numeric(area) < 6e+10 & as.numeric(area) > 80000)
Island_large<-Australia_map_Lambert %>% filter(as.numeric(area) < 6e+10 & as.numeric(area) > 600000) 
Island_Mainland<-Australia_map_Lambert %>% filter(as.numeric(area) > 60000000) 

National_grid_Lam$Island_large<-ifelse(lengths(st_intersects(National_grid_Lam , Island_large)) >0, 1, 0)
National_grid_Lam$Island_small<-ifelse(lengths(st_intersects(National_grid_Lam , Island_small)) >0, 1, 0)
nrow(National_grid_Lam)==31646

#1.3 GET AREA
National_grid_Lam$area <- st_area(National_grid_Lam)
National_grid_Lam$areakm2 <- National_grid_Lam$area*1e-6

#1.3 GET FADS
fad<-read_excel("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/National FAD.xlsx")
fad<-st_as_sf(fad, coords = c("Long_DD", "Lat_DD"), crs = 4283, agr = "constant")
fad %<>% filter(`Deployed in 2018/19?`=="Yes")
fad %<>% st_transform(3112)
National_grid_Lam$fad<-ifelse(lengths(st_intersects(National_grid_Lam , fad)) >0, 1, 0)
nrow(National_grid_Lam)==31646

#1.4 GET OFFSHORE DISTANCE
National_grid_Lam_centroid<-st_centroid_within_poly(National_grid_Lam)
Off_distances <- st_distance(National_grid_Lam_centroid, Island_Mainland)
colnames(Off_distances)<-Island_Mainland$FID

Off_distances<- as.data.frame(Off_distances)
National_grid_Lam_centroid$off_dist <- transform(Off_distances, min = do.call(pmin, Off_distances))$min/1000
nrow(National_grid_Lam_centroid)==31646

#1.5 GET LAT/LONG
National_grid_centroid <-National_grid_Lam_centroid %>% st_transform(4283)
National_grid_centroid$Long<-st_coordinates(National_grid_centroid)[,1]
National_grid_centroid$Lat <- st_coordinates(National_grid_centroid)[,2]
nrow(National_grid_centroid)==31646

## 2.0 EXPAND WITH DATES -------------
d1 <- as.Date(paste0("201801","01"), "%Y%m%d")
d2 <- as.Date(paste0("201812","01"), "%Y%m%d")
seq_dates<-seq(d1,d2,by="month")

National_grid_cent_exp<-expand_grid(National_grid_centroid, seq_dates)
nrow(National_grid_cent_exp)==31646*12 #Check that the expansion was correctly done

## 3.0 EXTRACT WEATHER AND SST -------------

#Weather
temp<-get_hs_wsDOWNLOAD(National_grid_cent_exp$seq_dates, National_grid_cent_exp$Long,  National_grid_cent_exp$Lat)
National_grid_cent_exp$wind<-temp$wind
National_grid_cent_exp$wave<-temp$wave

#SST
temp2<-get_sst_OneMonthAverage(National_grid_cent_exp$seq_dates, National_grid_cent_exp$Long,  National_grid_cent_exp$Lat)
National_grid_cent_exp$sst<-temp2$sst
National_grid_cent_exp$sst_slope<-temp$sst_slope

National_grid_cent_exp_w_sst <- National_grid_cent_exp
nrow(National_grid_cent_exp_w_sst)==12*31646

#pitstop - with all covariates except tc
#write.csv(National_grid_cent_exp_w_sst, "~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/National_grid_cent_exp_w_sst.csv")
National_grid_cent_exp_w_sst<-read.csv("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/National_grid_cent_exp_w_sst.csv")
## 4.0 COMBINE WITH RAMPS WITHIN 100km ------------- 
#distance between ramps and gridss
National_grid_cent_exp_w_sst <- st_as_sf(National_grid_cent_exp_w_sst, coords = c("Long", "Lat"), crs = 4283, agr = "constant")
nrow(National_grid_cent_exp_w_sst)==12*31646

#need to split this bit up as I'm running out of RAM. 
Ramps %<>% ungroup %>% mutate(split = round(runif(nrow(Ramps), 0.51,20.4)))
Ramps %<>% st_transform(3112)
#and only use unique values for grids
Grid_temp<-National_grid_cent_exp_w_sst %>% filter(seq_dates == "2018-01-01") %>% st_transform(3112)

for(i in 1: 20){
  Ramps_temp<-Ramps %>% filter(split == i)
  temp<-sf::st_distance(Ramps_temp, Grid_temp, tolerance = 10 )
  temp<-as.data.frame(drop_units(temp))
  colnames(temp)<-Grid_temp$ID
  temp$Ramp <-  Ramps_temp$RampID
  
  temp %<>% pivot_longer(cols = -Ramp, names_to = 'Grid', values_to = 'Distance')
  #only keep distances less than 100km one way distance
  temp %<>% filter(Distance < 100000)
  
  if(i==1){out<-temp}else{out<-rbind(out, temp)}
} #out now contains all ramp grid cobos within 100km of each other (perfect for our RUM)
length(unique(out$Ramp))

#out: for each ramp list the Grids within 100km. So varies in length
nrow(out) #400,464 Ramp Grid combos
length(unique(out$Ramp))==nrow(Ramps) #Every ramp is present
length(unique(out$Grid));length(unique(National_grid_cent_exp_w_sst$ID)) #but only 17877 grids out of the 31,646 are present
out %>% group_by(Ramp) %>% summarise(n=n()) %>% summarise(mean = mean(n)) #on average each ramp has 252 adjacent grids

#inner join og National_grid_cent_exp_w_sst (combos of each grid and month) and out (combos within 100km of each ramp and grid)
Nat_RUM_data<-National_grid_cent_exp_w_sst %>% mutate(ID = as.character(ID)) %>% 
  inner_join(., out, by = c("ID"= "Grid")) %>% st_drop_geometry()
length(unique(Nat_RUM_data$ID)) #note that we now are down to 17877 grids (just those in both data sets)
length(unique(Nat_RUM_data$Ramp)) #1246 ramps

#add tc
Nat_RUM_data %<>% mutate(tc = 2*0.54*Distance/1000)
#Nat_RUM_data contains the basis for our RUM. WE have rows for each rampxmonthxgrid within 100km

## 5.0 FIT RUM -------------
#Create a pitstop here in case it all goes to crap - which it did :)
#write.csv(Nat_RUM_data, "~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/Nat_RUM_data3.csv")
Nat_RUM_data<-read.csv("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/Nat_RUM_data3.csv")
#temp<-Nat_RUM_data
#Nat_RUM_data<-temp
#Model options
callibrated=1 #modify the tc coefficient to produce more offshore trips
area="WA" #NSW

#5.1 READ IN AND PREPARE MODELS
if(area=="NSW"){
  b<-read_excel("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/nsw_model_1.xlsx") %>% dplyr::select(Vars, Coef)
  v<-read_excel("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/nsw_model_1.xlsx", sheet = "v", col_names = rep("x", 6)) %>% 
    dplyr::select(-x...1 )
}else{
  b<-read_excel("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/wa_model.xlsx") %>% dplyr::select(Vars, Coef)
  b$Vars<-c("tc", "areakm2", "wave", "depth")
  v<-read_excel("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/wa_model.xlsx", sheet = "v", col_names = rep("x", 6)) %>% 
    dplyr::select(-x...1 )
  v$x...2<-c("tc", "areakm2", "wave", "depth")
  Nat_RUM_data$depth<-Nat_RUM_data$depth*-1
}

colnames(v)<- c("vars", v$x...2)
v %<>% dplyr::select(-vars)
vars<-b$Vars #use this to refer to the columns
setdiff(vars,  colnames(Nat_RUM_data)) #check that names line up

#5.2 GENERATE COEF DRAWS
n1<- 1000 #number of multivariate normal draws (all estimates)
mvn_b<-mvrnorm(n1,mu=b$Coef, Sigma = as.matrix(v)) #generate draws
if(callibrated==1){
  mvn_b[,1] <-mvn_b[,1]+0.05 
}
#5.3 CALCULATE OBSERVED UTILITY FOR EACH DRAW
Nat_RUM_data$tripId<- paste(Nat_RUM_data$Ramp, Nat_RUM_data$seq_dates, sep=".") #trip Id for later
rum_matrix<-as.matrix(Nat_RUM_data[, vars])
Vj<-matrix(NA, nrow = nrow(rum_matrix), ncol = n1) #a column for each draw

for(i in 1:n1) {
  #fixed utility
  Vj[ ,i] <-as.vector(rum_matrix%*%as.matrix(mvn_b[i,]))  #observed fixed utility
}

#5.4 EXTRACT PROBABILITIES FOR GRID RAMP TRIPS
#get exponential of utility
Vj<-exp(Vj)

#probabilities
Vj<-as.data.table(Vj)
cols<-colnames(Vj) #save the names of the draws for later
Vj$trip <- Nat_RUM_data$tripId #add tripId to group_by
Vj$ID <- Nat_RUM_data$ID #add tripId to group_by
Vj$Zone.type <- Nat_RUM_data$Zone.type #add zone type so NTRs can be removed
Vj$tc <- Nat_RUM_data$tc #add zone type so NTRs can be removed

#Update state NTRs to have zero utility
Vj[Vj$Zone.type=="State MP NTR" |Vj$Zone.type=="State NTR",cols]<-0

#calculate probability and of visit and logsums without NPZs
f1<- function(x) if(is.integer(x)){x}else{  x/sum(x)}
p_wo<- Vj[ , lapply(.SD, f1), by = trip, .SDcols = c(cols, "ID")]
p_wo<-p_wo[order(trip, ID)]
Nat_RUM_data %<>% arrange(tripId, ID)
table(Nat_RUM_data$ID==p_wo$ID)#they do now

logsum_wo<- Vj[ , lapply(.SD, sum), by = trip, .SDcols = cols] 
p_wo<- p_wo[ , trip:=NULL]
p_wo<- p_wo[ , ID:=NULL]
p_wo<-as.matrix(p_wo) 

#get probability of visit summary statistics by row
#p_wo_mat<-  as.matrix(p_wo[ , ..cols]) #much faster in matrix form
Nat_RUM_data$p.mean_wo<- rowQuantiles(p_wo, probs = 0.5)
Nat_RUM_data$p.low_wo<- rowQuantiles(p_wo, probs = sqrt(0.025)) 
Nat_RUM_data$p.upp_wo<- rowQuantiles(p_wo, probs = 1-sqrt(0.025))
#write.csv(Nat_RUM_data, "NRDpitstop1.csv")
#write.csv(logsum_wo, "LSpitstop1.csv")

#now remove National Park Zones and re-estimate probabilities
Vj[Vj$Zone.type=="National Park Zone" ,cols]<-0
p_wo<- Vj[ , lapply(.SD, f1), by = trip, .SDcols = c(cols, "ID")]
p_wo<-p_wo[order(trip, ID)]
Nat_RUM_data %<>% arrange(tripId, ID)
table(Nat_RUM_data$ID==p_wo$ID)#they do now
logsum_w<- Vj[ , lapply(.SD, sum), by = trip, .SDcols = cols] #18876 trips

#get probability of visit summary statistics by row
p_wo<- p_wo[ , trip:=NULL]
p_wo<- p_wo[ , ID:=NULL]
p_wo<-as.matrix(p_wo)
Nat_RUM_data$p.mean_w<- rowQuantiles(p_wo, probs = 0.5)
Nat_RUM_data$p.low_w<- rowQuantiles(p_wo, probs = sqrt(0.025)) #saving outputs back to p
Nat_RUM_data$p.upp_w<- rowQuantiles(p_wo, probs = 1-sqrt(0.025)) #Note below about weird probs
#Note that when i combine these probabilites with those from the trip numbers normal confidence intervals are innappropriate
#as the probability of getting below the lower end of both CIs (prob trips and trip numbers) is 0.025*0.025 = 0.000625
#by using the sqrt(0.025) we are correcting for this and ensuring the final CIs are a true 95% CI range. sqrt(0.025)*sqrt(0.025) = 0.025
#write.csv(Nat_RUM_data, "NRDpitstop2.csv")
#write.csv(logsum_w, "LSpitstop2.csv")

#5.5 EXTRACT WELFARE IMPACTS FOR TRIPS
W<-logsum_w[ ,-1] -  logsum_wo[ ,-1]
W<-t(W)/mvn_b[,1 ]
W<-t(W)

#and get summaries across draws for
W <-as.matrix(W)
logsum_w$w.mean<- rowQuantiles(W, probs = 0.5)
logsum_w$w.low<- rowQuantiles(W, probs = sqrt(0.025)) #saving outputs back to p
logsum_w$w.upp<- rowQuantiles(W, probs = 1-sqrt(0.025))

Nat_RUM_data %<>% left_join(logsum_w[ , c("trip", "w.mean", "w.upp", "w.low")], by = c("tripId" = "trip"))
#write.csv(Nat_RUM_data, "NRDpitstop3.csv")


## 6.0 JOIN TO RAMPS -------------
#this section merges estimates of trip numbers at each ramp each month (Trips_Ramp_out) with estimates of the 
#probability of visiting each grid (by month) given a particular boat ramp was selected (Nat_RUM_data). 

#THIS IS PRODUCED FROM THE "TripAllocation.R" SCRIPT
Trips_Ramp_out<-read.csv("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/Trips_Ramp_out" ) 

#Check that merge variables match
Nat_RUM_data$Month<-lubridate::month(Nat_RUM_data$seq_dates, label = T)
setdiff(Trips_Ramp_out$Month,Nat_RUM_data$Month)
setdiff(Nat_RUM_data$Month, Trips_Ramp_out$Month) #Month matches
Nat_RUM_data$Month<-factor(Nat_RUM_data$Month, ordered = F)

setdiff(Trips_Ramp_out$RampID,Nat_RUM_data$Ramp)
setdiff(Nat_RUM_data$Ramp, Trips_Ramp_out$RampID) 
Nat_RUM_data %<>% filter(!Ramp %in% c(1275,28))

National_rum<-Nat_RUM_data %>% left_join(.,Trips_Ramp_out, by = c("Month" = "Month","Ramp" = "RampID" ), 
                                         suffix = c("", "_Ramp") )

#multiply by trips for welfare impacts
National_rum_W<-National_rum %>% distinct(.,tripId, .keep_all = TRUE) %>%
  mutate(across(starts_with("TRIPS_ramp_prop"), ~.x*w.mean, .names = "welf_mean_{.col}" ),
         across(starts_with("upp_TRIPS"), ~.x*w.upp, .names = "welf_upp_{.col}" ),
         across(starts_with("low_TRIPS"), ~.x*w.low, .names = "welf_low_{.col}" ))
setdiff(National_rum_W$Ramp, Ramps$RampID)
setdiff(Ramps$RampID, National_rum_W$Ramp)
National_rum_W %<>% left_join(Ramps[,c("RampID", "Network")], by = c("Ramp" = "RampID"), suffix = c("","_ramp"))

#multiplied trip numbers by corresponding probabilities. Note that because we adjusted the probs above to 
#sqrt(0.025) we should now have accurate 95% CIs. 

National_rum %<>% mutate(across(starts_with("TRIPS_ramp_prop"), ~.x*p.mean_w, .names = "mean_w_{.col}" ),
                         across(starts_with("TRIPS_ramp_prop"), ~.x*p.mean_wo, .names = "mean_wo_{.col}" ),
                         across(starts_with("upp_TRIPS"), ~.x*p.upp_w, .names = "upp_w_{.col}" ),
                         across(starts_with("upp_TRIPS"), ~.x*p.upp_wo, .names = "upp_wo_{.col}" ),
                         across(starts_with("low_TRIPS"), ~.x*p.low_w, .names = "low_w_{.col}" ),
                         across(starts_with("low_TRIPS"), ~.x*p.low_wo, .names = "low_wo_{.col}" ))

#Sum by grid (over months and ramps) to get use at grid level
National_rum_grid<-National_rum %>% group_by(ID) %>%  #was jsut grouping by ID
  summarise(across(starts_with(c("mean_w_","upp_w_","low_w_","mean_wo_", "upp_wo_", "low_wo_")), sum))

## 7.0 GENERATE OUTPUTS -------------
BrRegions<-st_read("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/BRRegions.shp") #bringing this in for some plots

#merge in details of grids
National_grid<-National_grid_Lambert2 %>% st_transform(4283)
National_rum_grid<-right_join(National_grid, National_rum_grid, by = c("ID"))

#merge in BrRegions
National_rum_grid<-st_join(National_rum_grid, BrRegions)
National_rum_grid_df<-National_rum_grid %>% st_drop_geometry()

#save for later
#write.csv(National_rum_grid_df, "~/R/Poject-NESPD6-MBH-NationalRUM/Output Datasets/WA_base_National_rum_grid_df.csv")
#write.csv(National_rum_W, "~/R/Poject-NESPD6-MBH-NationalRUM/Output Datasets/WA_base_National_rum_W.csv")

National_rum_grid_df %>% summarise(mean_wo_TRIPS_ramp_prop_Gravity=sum(mean_w_TRIPS_ramp_prop_Gravity) )

National_rum_grid_df %>% filter(Type=="AMP") %>% 
  summarise(across(starts_with(c("mean_w_","upp_w_","low_w_", "mean_wo_", "upp_wo_", "low_wo_")), sum)) 

#TRIPS BY NETWORK
u_y0<-National_rum_grid_df %>% filter(Type=="AMP") %>% group_by(Network) %>%
  summarise(across(starts_with(c("mean_w_","upp_w_","low_w_", "mean_wo_", "upp_wo_", "low_wo_")), sum)) %>%
  ggplot(.) + geom_bar(aes(x = Network, y = mean_w_TRIPS_ramp_prop_Gravity, fill = Network), stat="identity", colour = "black") +
  geom_errorbar(aes(x = Network, ymin = low_w_low_TRIPS_Gravity, ymax = upp_w_upp_TRIPS_Gravity), width = 0.2)+ 
  my_theme +
  ylim(c(0,90000)) +
  scale_x_discrete(drop = T, labels = function(Network) str_wrap(Network, width = 5)) +
  ylab("Annual number of trips\nin the AMPs") +
  annotate("text", x=0.8, y=80000, label="(A) y = 0", size = 4, fontface=2)+
  scale_fill_manual(values = colours)

u_y0.5<-National_rum_grid_df %>% filter(Type=="AMP") %>% group_by(Network) %>%
  summarise(across(starts_with(c("mean_w_","upp_w_","low_w_", "mean_wo_", "upp_wo_", "low_wo_")), sum)) %>%
  ggplot(.) + geom_bar(aes(x = Network, y = mean_w_TRIPS_ramp_prop_Gravity, fill = Network), stat="identity", colour = "black") +
  geom_errorbar(aes(x = Network, ymin = low_w_low_TRIPS_Gravity, ymax = upp_w_upp_TRIPS_Gravity), width = 0.2)+ 
  my_theme +
  ylim(c(0,90000)) +
  scale_x_discrete(drop = T, labels = function(Network) str_wrap(Network, width = 5)) +
  ylab("Annual number of trips\nin the AMPs") +
  annotate("text", x=0.8, y=80000, label="(B) y = 0.5", size = 4, fontface=2)+
  scale_fill_manual(values = colours)

y.label <- grid::textGrob("Estimated number of trips\nin AMPs (annual)", gp=grid::gpar(fontsize=14), rot=90)
x.label <- grid::textGrob("Network", gp=grid::gpar(fontsize=14))
legend <- gtable_filter(ggplotGrob(u_y0.5), "guide-box")

u <-grid.arrange(arrangeGrob(u_y0 + theme(legend.position="none")+ xlab(NULL) + ylab(NULL), 
                             u_y0.5 + theme(legend.position="none")+ xlab(NULL) + ylab(NULL),
                             nrow = 2,
                             left = y.label,
                             bottom = x.label),
                 legend, 
                 widths=grid::unit.c(unit(1, "npc") - legend$width, legend$width), 
                 nrow=1)

setwd(fig_sample_dir)
ggsave(file="~/R/Poject-NESPD6-MBH-NationalRUM/u_nsw.jpeg",  u, height = a4_width*1, width = a4_width, units  ="mm", dpi = 300)

#PLOT BY AMP
MainMPs<-c("Dampier", "Two Rocks","Geographe","South-west Corner","Western Eyre", "Murray",  "Hunter" , "Solitary Islands")

mpy0<-National_rum_grid_df %>% filter(Type=="AMP") %>% group_by(Network, MP.name) %>%
  summarise(across(starts_with(c("mean_w_","upp_w_","low_w_", "mean_wo_", "upp_wo_", "low_wo_")), sum)) %>%
  filter(MP.name %in% MainMPs) %>%
  mutate(MP.name = factor(MP.name, levels = MainMPs)) %>%
  ggplot(.) + geom_bar(aes(x = MP.name, y = mean_w_TRIPS_ramp_prop_Gravity, fill = Network), stat="identity", colour = "black") +
  geom_errorbar(aes(x = MP.name, ymin = low_w_low_TRIPS_Gravity, ymax = upp_w_upp_TRIPS_Gravity), width = 0.2)+ 
  my_theme +
  xlab("Marine Park")+
  ylim(c(0,5000)) +
  scale_x_discrete(drop = T, labels = function(MP.name) str_wrap(MP.name, width = 5)) +
  annotate("text", x=0.9, y=4500, label="(A) y = 0", size = 4, fontface=2)+
  ylab("Annual number of trips\nin the AMPs") +
  scale_fill_manual(values = colours[3:6])

mpy0.5<-National_rum_grid_df %>% filter(Type=="AMP") %>% group_by(Network, MP.name) %>%
  summarise(across(starts_with(c("mean_w_","upp_w_","low_w_", "mean_wo_", "upp_wo_", "low_wo_")), sum)) %>%
  filter(MP.name %in% MainMPs) %>%
  mutate(MP.name = factor(MP.name, levels = MainMPs)) %>%
  ggplot(.) + geom_bar(aes(x = MP.name, y = mean_w_TRIPS_ramp_prop_Grvt__5, fill = Network), stat="identity", colour = "black") +
  geom_errorbar(aes(x = MP.name, ymin = low_w_low_TRIPS_Grvt__5, ymax = upp_w_upp_TRIPS_Grvt__5), width = 0.2)+ 
  my_theme +
  xlab("Marine Park")+
  ylim(c(0,5000)) +
  scale_x_discrete(drop = T, labels = function(MP.name) str_wrap(MP.name, width = 5)) +
  annotate("text", x=0.9, y=4500, label="(B) y = 0.5", size = 4, fontface=2)+
  ylab("Annual number of trips\nin the AMPs") +
  scale_fill_manual(values = colours[3:6])

y.label <- grid::textGrob("Estimated number of trips\n(annual)", gp=grid::gpar(fontsize=14), rot=90)
x.label <- grid::textGrob("Marine Park", gp=grid::gpar(fontsize=14))
legend <- gtable_filter(ggplotGrob(u_y0.5), "guide-box")

u <-grid.arrange(arrangeGrob(mpy0 + theme(legend.position="none")+ xlab(NULL) + ylab(NULL), 
                             mpy0.5 + theme(legend.position="none")+ xlab(NULL) + ylab(NULL),
                             nrow = 2,
                             left = y.label,
                             bottom = x.label),
                 legend, 
                 widths=grid::unit.c(unit(1, "npc") - legend$width, legend$width), 
                 nrow=1)

#GRIDDED MAPS BY LOCATION
Ningaloo<-National_rum_grid %>% filter(Location=="Ningaloo") %>% 
  ggplot(.)+
  xlab("Longitude")+ylab("Latitude")+
  geom_sf(aes(fill = mean_w_TRIPS_ramp_prop_Gravity))+my_theme+
  labs(fill="Trips\nannually")

Capes<-National_rum_grid %>% filter(Location=="Capes") %>% 
  ggplot(.)+
  geom_sf(aes(fill = mean_w_TRIPS_ramp_prop_Gravity))+my_theme+
  xlab("Longitude")+ylab("Latitude")+
  labs(fill="Trips\nannually")

Coffs<-National_rum_grid %>% filter(Location=="Coffs Harbour") %>% 
  ggplot(.)+
  geom_sf(aes(fill = mean_w_TRIPS_ramp_prop_Gravity))+my_theme+
  xlab("Longitude")+ylab("Latitude")+
  labs(fill="Trips\nannually")

TwoRocks<-National_rum_grid %>% filter(Location=="Two Rocks") %>% 
  ggplot(.)+
  geom_sf(aes(fill = mean_w_TRIPS_ramp_prop_Grvt__5))+my_theme+
  xlab("Longitude")+ylab("Latitude")+
  labs(fill="Trips\nannually")

VictorHarbour<-National_rum_grid %>% filter(Location=="Victor Harbour") %>% 
  ggplot(.)+
  geom_sf(aes(fill = mean_wo_TRIPS_ramp_prop_Gravity))+my_theme+
  xlab("Longitude")+ylab("Latitude")+
  labs(fill="Trips\nannually")

Bicheno<-National_rum_grid %>% filter(Location=="Bicheno") %>% 
  ggplot(.)+
  geom_sf(aes(fill = mean_w_TRIPS_ramp_prop_Gravity))+my_theme+
  xlab("Longitude")+ylab("Latitude")+
  labs(fill="Trips\nannually")

plot_grid(Ningaloo, Coffs, labels = c("A", "B") )
plot_grid(TwoRocks, Bicheno, labels = c("C", "D") )

MP_use<-National_rum_grid_df %>% filter(Type=="AMP") %>% group_by(Network, MP.name) %>%
  summarise(across(starts_with(c("mean_w_","upp_w_","low_w_", "mean_wo_", "upp_wo_", "low_wo_")), sum)) %>% arrange(Network, MP.name)

write.csv(MP_use, "MP_use.csv")

#Welfare impacts
National_rum_W %>% 
  summarise(across(starts_with(c("welf_mean","welf_upp","welf_low")), sum))

w_y0<-National_rum_W %>% group_by(Network_ramp) %>%
  summarise(across(starts_with(c("welf_mean","welf_upp","welf_low")), sum)) %>%
  ggplot(.) + geom_bar(aes(x = Network_ramp, y = welf_mean_TRIPS_ramp_prop_Gravity, fill = Network_ramp), stat="identity",  colour = "black") +
  geom_errorbar(aes(x = Network_ramp, ymin = welf_low_low_TRIPS_Gravity, ymax = welf_upp_upp_TRIPS_Gravity), width = 0.2)+ 
  my_theme +
  #ylim(c(0,1000)) +
  annotate("text", x=0.8, y=70000, label="(A) y = 0", size = 4, fontface=2)+
  scale_x_discrete(drop = T, labels = function(Network_ramp) str_wrap(Network_ramp, width = 5)) +
  ylab("Annual welfare impacts on\nrecreational fishers ($)") +
  scale_fill_manual(values = colours)

w_y0.5<-National_rum_W %>% group_by(Network_ramp) %>%
  summarise(across(starts_with(c("welf_mean","welf_upp","welf_low")), sum)) %>%
  ggplot(.) + geom_bar(aes(x = Network_ramp, y = welf_mean_TRIPS_ramp_prop_Grvt__5, fill = Network_ramp), stat="identity", colour = "black") +
  geom_errorbar(aes(x = Network_ramp, ymin = welf_low_low_TRIPS_Grvt__5, ymax = welf_upp_upp_TRIPS_Grvt__5), width = 0.2)+ 
  my_theme +
  annotate("text", x=0.8, y=70000, label="(B) y = 0.5", size = 4, fontface=2)+
  #ylim(c(0,1000)) +
  scale_x_discrete(drop = T, labels = function(Network_ramp) str_wrap(Network_ramp, width = 5)) +
  ylab("Annual welfare impacts on\nrecreational fishers ($)") +
  scale_fill_manual(values = colours)


y.label <- grid::textGrob("Estimated welfare impacts on\nrecreational fishers ($ annual)", gp=grid::gpar(fontsize=14), rot=90)
x.label <- grid::textGrob("Network", gp=grid::gpar(fontsize=14))
legend <- gtable_filter(ggplotGrob(u_y0.5), "guide-box")

w <-grid.arrange(arrangeGrob(w_y0 + theme(legend.position="none")+ xlab(NULL) + ylab(NULL), 
                             w_y0.5 + theme(legend.position="none")+ xlab(NULL) + ylab(NULL),
                             nrow = 2,
                             left = y.label,
                             bottom = x.label),
                 legend, 
                 widths=grid::unit.c(unit(1, "npc") - legend$width, legend$width), 
                 nrow=1)

setwd(fig_sample_dir)
ggsave(file="~/R/Poject-NESPD6-MBH-NationalRUM/u_nsw.jpeg",  u, height = a4_width*1, width = a4_width, units  ="mm", dpi = 300)

# Plot by location - just for calibration
National_rum_grid_df %>% filter(Type=="AMP" & !is.na(Location)) %>% group_by(Location) %>%
  summarise(across(starts_with(c("mean_w_","upp_w_","low_w_", "mean_wo_", "upp_wo_", "low_wo_")), sum)) %>%
  ggplot(.) + geom_bar(aes(x = Location, y = mean_w_TRIPS_ramp_prop_Gravity), stat="identity") +
  geom_errorbar(aes(x = Location, ymin = low_w_low_TRIPS_Gravity, ymax = upp_w_upp_TRIPS_Gravity), width = 0.2)+ 
  my_theme +
  #ylim(c(0,1000)) +
  #scale_x_discrete(drop = T, labels = function(Network) str_wrap(Network, width = 5)) +
  ylab("Annual number of trips\nin the AMPs") +
  scale_fill_manual(values = colours)

test<-National_rum_grid_df %>% filter(Location=="Bicheno")

# Plot of map - just for callibration
National_rum_grid %>% filter(Location=="Two Rocks") %>% 
  ggplot(.)+
  geom_sf(aes(fill = mean_wo_TRIPS_ramp_prop_Gravity))+my_theme+
  labs(fill="Predicted annual trip\nnumbers (y = 0.5)")


#Troubleshooting
test<- National_rum_grid_df %>% filter(Network=="North")


###
##ARCHIVE----------------

## PREPARATION----
wa_trips$Date<-as_date(wa_trips$eventdate)
wa_trips_national$Date<-as.Date(wa_trips_national$Date)
wa_test<-wa_trips[runif(100000,1,400000),]

nsw_trips$Date<-as_date(nsw_trips$StartDate)
nsw_test<-nsw_trips[1,]


## EXTRACT SST----
start.time <- Sys.time()
wa_trips_national$sst_month<-get_sst_OneMonthAverage(wa_trips_national$Date,wa_trips_national$x,  wa_trips_national$y)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

waread.csv("wa_trips_national.csv")

## EXTRACT WEATHER----
start.time <- Sys.time()
wa<-get_hs_ws(test$Date,test$x,  test$y)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


wa_trips$hs<-wa$hs ; wa_trips$wspeed<-wa$wspeed

nsw_trips_test<- nsw_trips[1:100,]

start.time <- Sys.time()
nsw_weather<-get_hs_ws_day(nsw_trips$Date, nsw_trips$Long,  nsw_trips$Lat)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
nsw_trips %<>% bind_cols(nsw_trips, as.data.frame(nsw_weather))

write.csv(nsw_trips, "nsw_trips_expanded_weather.csv")

National_grid_expanded<- National_grid_centroid_expanded %>% filter(year(seq_dates)=="2018")

start.time <- Sys.time()
national_weather<-get_hs_ws_day(National_grid_expanded$seq_dates, National_grid_expanded$Long,  National_grid_expanded$Lat)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

write.csv(National_grid_expanded, "National_grid_expanded.csv")
National_grid_expanded<-read.csv("National_grid_expanded.csv")

National_grid_expanded %<>% arrange(seq_dates)
months<-unique(National_grid_expanded$seq_dates)

N1<-National_grid_expanded[National_grid_expanded$seq_dates%in%months[1:3],]
N2<-National_grid_expanded[National_grid_expanded$seq_dates%in%months[4:6],]
N3<-National_grid_expanded[National_grid_expanded$seq_dates%in%months[7:9],]
N4<-National_grid_expanded[National_grid_expanded$seq_dates%in%months[10:12],]


W1<-get_hs_wsDOWNLOAD(N1$seq_dates, N1$Long,  N1$Lat)
W2<-get_hs_wsDOWNLOAD(N2$seq_dates, N2$Long,  N2$Lat)
W3<-get_hs_wsDOWNLOAD(N3$seq_dates, N3$Long,  N3$Lat)
W4<-get_hs_wsDOWNLOAD(N4$seq_dates, N4$Long,  N4$Lat)

N1$wind<-W1$wind
N1$wave<-W1$wave
N2$wind<-W2$wind
N2$wave<-W2$wave
N3$wind<-W3$wind
N3$wave<-W3$wave
N4$wind<-W4$wind
N4$wave<-W4$wave

National_grid_expanded_w<-rbind(N1, N2, N3, N4)

National_grid_expanded_w

write.csv(National_grid_expanded_w, "Weathered.csv")


weather<-rbind(as.data.frame(W1),as.data.frame(W2), as.data.frame(W3), as.data.frame(W4.1),
               as.data.frame(W4.2), as.data.frame(W4.3), as.data.frame(W4.4))
National_grid_centroid_expanded %<>% bind_cols(weather)


data<-st_read("~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/National_grid_Lambert_centroid.shp")

data %<>% st_drop_geometry()



National_grid_centroid_expanded %<>% left_join(data[ , c("id", "areakm2", "off_dist")], by = c("id" = "id") )
write.csv(National_grid_centroid_expanded, "~/R/Poject-NESPD6-MBH-NationalRUM/Data_nectar/National_grid_weather.csv" )

## EXPORT DATA FOR STATA---------




#not exactly clear about what this is doing
test<-read.csv("wa_filler.csv")

missing<-wa_trips_national %>% filter(is.na(wave))
wa_trips_national %<>% filter(!is.na(wave))
missing %<>% arrange(X)
test %<>% arrange(X)

missing$wave<-test$wave

wa_trips_national<-rbind(wa_trips_national, missing)


temp<-left_join(wa_trips_national, test, by = "X")

write.csv(wa_trips_national, "wa_trips_national.csv")







### Things I need to try and do
## Make recruitment equilibrate to natural mortality - need to figure out what numbers mean that the number of recruits 
# produced each year equals the number of fish that die every year
# Can also work out a and b using the equations in Eva's sea cucumber paper
## Also need to determine the size of the cells and the swimming speed so we can figure out how far the fish can move in
# one time step
# 50% kernal density has a mean of 2.3km and 95% is 8.6 
# In the lagoon they moved on average 2.92km in 30 days in the lagoon and 4.21km on the reef slope (Pillans 2014)
# Also show seasonal migration based on spawning
# Also move from the sargassum to the reefs as part of an ontogenic shift - need to have a different matrix for the one year old fish to make them move to the reef
## Will need to make sure when you put mortaility back in that you're subsetting the population between the ones that are in/outside a sanctuary zone and apply the correct mortality
# this will have to be based on the intersection of our grid with a shape file of the sanctuary zones

#### OLD THINGS THAT I'VE TRIED AND DON'T USE ANYMORE ####

##INITIAL POPULATION MODEL##
# Bring in the habitat data 
setwd(sp_dir)
habitat <- st_read("SeamapAus_WA_DPAW_marine_habitatsPolygon.shp")%>%
  st_transform(4283)%>%
  st_make_valid()

ningaloo_habitat <- st_crop(habitat, xmin=112.5, xmax=114.3, ymin=-24, ymax=-21)
plot(ningaloo_habitat$geometry, add=TRUE)

# Set the different habitat types to be a value between 0 and 4
# 0 is places fish are unlikely to be found e.g.  mudflats and saltmarsh
# 2 is unnatractive places for fish including mobile sand, deep water habitat, pelagic
# 5 is macroalgae and mangroves where you're likely to get juvenile fish not targeted by fishers
# 10 is coral reef and subtidal bare reef (subtidal and intertidal)

ningaloo_habitat_cls <- ningaloo_habitat%>%
  mutate(classification = ifelse(str_detect(SM_HAB_CLS, "Mudflat|Saltmarsh"), 0,
                                 ifelse(str_detect(SM_HAB_CLS, "Mobile sand|Pelagic"), 2,
                                        ifelse(str_detect(SM_HAB_CLS, "Macroalgae|Mangroves"), 5, 10))))

# Create an intersection of the habitat and the grid cells to see what overlaps where, then
# calculate the area of each polygon
intersection <- st_intersection(ningaloo_habitat_cls, water)
intersection$area <- st_area(intersection)

#calculate the area of the water grid cells
water$area <- st_area(water)

#calculate the percentage of each habitat type found in the different water grid cells
hab_perc <- intersection%>%
  st_drop_geometry()%>%
  as.data.frame()%>%
  mutate(classification=factor(classification))%>%
  mutate(ID=factor(ID))%>%
  dplyr::select(classification, ID, area)%>%
  group_by(ID, classification)%>%
  summarise(total_area = sum(area))

# Join the habitat percentage data to the grid cell data and calculate percent reef
per_reef <- water%>%
  st_drop_geometry()%>%
  full_join(hab_perc)%>%
  filter(classification==10)%>%
  mutate(percent_reef=(total_area/area)*100)%>%
  full_join(water, by="ID")%>%
  dplyr::select(ID, percent_reef)

# 
# # For cells where we are missing data about habitat take raster data and then calculate
# # the average rugosity, for now we'll just say anything over XXX is reef
# setwd(sp_dir)
# bathy <- raster("Carnarvon_Shelf_Bathymetry_3_2008_3m_cog.tiff")
# bathy <- flip(bathy, direction="y")
# proj4string(bathy) #check projection
# 
# 
# 
# water_ras <- st_transform(water, "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# 
# grd_transform <- st_transform(grd, "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# 
# plot(bathy)
# plot(water_ras$water, add=TRUE)
# 
# st_crs(grd_transform)

# This needs to be habitat matrix*coefficient + distance*coefficient in some form or another. Need to figure out how to decide
# what the coefficients should be so we can weight the importance of habitat relative to distance

## Settlement
Settlement <- matrix(runif(143, 1, 10),nrow=143) # Random for now

# Need to normalise so that it sums to 1 so we're not adding or subtracting any fish anywhere
Settlement_norm <- (Settlement[,1] = (Settlement[,1]/sum(Settlement[,1])))
Settlement_norm <- as.matrix(Settlement_norm)

# ## Movement 
# Movement <- array(1, dim= c(143,143,1))
# for(c in 1:143){
#   Movement_norm[,c,] <- (Movement[,c,] = (Movement[,c,]/sum(Movement[,c,])))
# }
# Movement_norm <- as.matrix(Movement_norm)

## Set up the model
Pop <- array(0, dim = c(143,100,10)) #This is our population

# Fill the array with random numbers of fish for now
for(d in 1:10){
  Pop[,1,d] <- matrix(runif(143, 1, 10000))
}
Pop <- round(Pop, digits=0)
Recruits <- 2000


##Population Model###
## Run the model
for(t in 2:5){
  for(A in 1:dim(Pop)[3]){
    if(A==1){Pop[,t-1,A] <- Recruits}
    else {Pop[,t,A] <-  Pop[,t-1,A-1]}
    #else if (d>=6) {Days[,t,d] <-  Movement[,,1] %*% Days[,t-1,d-1]}
    #else {Days[,t,d] <- Movement[,,1] %*% Days[,t-1,d-1]}
  }
  
  
  # Mortality
  for(A in 1:dim(Pop)[3]){
    sa <- 1/(1+(exp(-log(19)*((A-A95/A95-A50)))))
    Pop[,t,A] <- Pop[,t,A]*exp(-M)*exp(-sa*Fishing)
  }
  
  
  # Ricker recruitment model with maturation
  Recs <- matrix(0, nrow=1, ncol=1) #Blank matrix to add the recruits to
  Recruitment <- as.matrix(colSums(Pop[,t,], dims=1)) 
  for(A in 1:dim(Recruitment)[1]){
    Mature <- 1/(1+(exp(-log(19)*((A-M95)/(M95-M50))))) #Number of mature individuals in each age class
    S <- colSums(Recruitment)*Mature #Spawning stock
    Rec <- a*S*exp(-b*S) #Number of recruits from that age class
    Recs <- rbind(Recs, Rec) #Combine recruits from each age class into one dataframe
  }
  R <- colSums(Recs) #Add up the number of recruits produced from all age classes
  Recruits <- as.matrix(Settlement_norm*R) #Distribute the recruits among the 143 sites
  
  # Plotting
  water$pop <- rowSums(Pop[,t,])
  
  water <- water%>%
    mutate(Population = ifelse(pop < 1000, "<1000",
                               ifelse (pop>1000 & pop<5000, "1000-5000",
                                       ifelse (pop>5000 & pop<10000, "5000-10000",
                                               ifelse (pop>10000 & pop<15000, "10000-15000",
                                                       ifelse(pop>15000 & pop<20000, "15000-20000",
                                                              ifelse(pop>20000 & pop<25000, "20000-25000",
                                                                     ifelse(pop>25000 & pop<30000, "25000-30000", ">30000"
                                                                            
                                                                     ))))))))%>%
    mutate(Population = factor(Population))
  water$Population <- fct_relevel(water$Population, "<1000",  "1000-5000", "5000-10000", "10000-15000",
                                  "15000-20000", "20000-25000", "25000-30000", ">30000")
  
  print(map <- ggplot(water)+
          geom_sf(aes(fill=Population))+
          scale_fill_manual(name="Population", values= cols, drop=FALSE)+
          theme_void())
  Sys.sleep(3)
  
}

## Making the fish move south
# Matrix that makes more southerly grid cells more preferable
latlong <- points[,1:2]
bearings <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (i in 1:NCELL){
  b <- bearing(latlong[i,], latlong)
  bearings[i, 1:NCELL] <- b
}

# We want all the bearings between 90-180 and -90 - -180 to be positive and everything else to be negative
# These are the bearings that indicate a southerly movement
bearings <- as.data.frame(bearings)%>%
  mutate_all(~ ifelse(.<(-90), .*-1,
                      ifelse(.>-90 & .<0, .*1,
                             ifelse(.>0 & .<90, .*-1, .))))

bearings <- as.matrix(bearings)

#### ARCHIVED FUNCTIONS ####

#### Mortality ####

mortality.func <- function (Age, mort.50, mort.95, Nat.Mort, NTZ, Effort, Cell, Max.Cell, Month, Year, Population){
  
  if(BurnIn){
    
    survived <- Population*exp(-Nat.Mort)
    
  }
  sa <- 1/(1+(exp(-log(19)*((Age-mort.95/mort.95-mort.50))))) # This puts in size selectivity
  
  tot.survived <- NULL

  for(Cell in 1:Max.Cell){

  survived <- Population[Cell,Month-1,Age]*exp(-sa*Effort[Cell,Month-1,Year])*exp(-Nat.Mort)
  tot.survived <- rbind(tot.survived, survived)

  } # End fishing mortality for each cell
  
  return(tot.survived)
}

#### Movement ####

movement.func <- function (Age, Month, Population, Max.Cell, Adult.Move, Juv.Move){
  
  All.Movers <- NULL
  
  ## Juvenile Movement
  if(Age<=4){
    
    Juv.Pop <- matrix(Population[ , Month-1, Age-1]) # This gives you the fish in all the sites at time step Month-1 of age A-1
    # Juv.Movers <- sum(Juv.Pop)
    # Juv.Moving <- NULL
    # Juv.Moved <- NULL
    
    Juv.Pop2 <- sapply(seq(Max.Cell), function(Cell){
      Juv.Pop2 <- as.matrix(Juv.Move[Cell, ] * Juv.Pop[Cell,1]) #This should give you the number of fish that move from each cell to all the other sites
      Juv.Pop2 <- t(Juv.Pop2) # Gives you a row with each cell being the number of fish that moved from cell 1 to all other cells 
      # Juv.Moving <- rbind(Juv.Moving, Juv.Pop2) #This should give you an array with rows representing the fish that move from each site to all the other sites 
    }
    )
    
    Juv.Movement2 <- colSums(Juv.Pop2)
    Juv.Moved <- sum(Juv.Movement2)
    # print(isTRUE(all.equal(Juv.Movers, Juv.Moved)))
    # if((isTRUE(all.equal(Juv.Movers, Juv.Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
    # {print(Juv.Movers)
    #   print(Juv.Moved)} else{ } Can Put This Check Back in Periodically to Confirm model Functioning
    
    All.Movers <- cbind(All.Movers, Juv.Movement2)
    
    return(All.Movers)
    
    
  } else {  
    
    ## Adult Movement
    Pop <- matrix(Population[ , Month-1, Age-1]) #This gives you the fish in all the sites at timestep t-1 of age A-1
    # Movers <- sum(Pop)
    # Moving <- NULL
    # Moved <- NULL
    
    Pop2 <- sapply(seq(Max.Cell), function(Cell){
      Pop2 <- as.matrix(Adult.Move[Cell, ] * Pop[Cell,1]) #This should give you the number of fish that move from site s to all the other sites
      Pop2 <- t(Pop2)
      # Moving <- rbind(Moving, Pop2) #This should give you an array with 143 rows representing the fish that move from each site to all the other sites 
    }
    ) #End bracket for movement in each cell
    
    Movement2 <- colSums(Pop2)
    Moved <- sum(Movement2)
    # print(isTRUE(all.equal(Movers, Moved)))
    # if((isTRUE(all.equal(Movers, Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
    # {print(Movers)
    #   print(Moved)} else{ }
    
    All.Movers <- cbind(All.Movers, Movement2)
    
    return(All.Movers)
    
  } #End adult movement 
  
}

#### Recruitment Function ####

## Recruitment

recruitment.func <- function(Population, Age, mat.95, mat.50, settlement, Max.Cell, relationship){
  adults <- Population[ ,1,Age]
  tot.recs <- data.frame(matrix(0, nrow = Max.Cell, ncol = dim(Population)[3]))

  for(Age in 1:dim(Population)[3]){

    mature <- 1/(1+(exp(-log(19)*((Age-mat.95)/(mat.95-mat.50))))) #Proportion of mature individuals in each age class

    recs <- sapply(seq(Max.Cell), function(Cell){
      S <- adults[Cell]*mature #Spawning stock
      rec <- S*relationship #Number of recruits from that age class -  will need to put the SST thing in here
      recs
    })

    tot.recs[ ,Age] <- recs

  }
  
  settle.recs <- tot.recs
  
  ## This was a loop before
  settle.recs <- sapply(1:dim(Population)[3], function(Age){
    temp <- sum(settle.recs[,Age])
    settle.recs[ ,Age] <- temp*settlement[,2]
  })
  
  settled <- rowSums(settle.recs)
  return(settled)
}

# Then get added to January Population 














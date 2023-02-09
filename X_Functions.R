###################################################

# Script that has all the functions
# Functions are only needed when running the model
# Have functions for movement, mortality, recruitment
# and for plotting

###################################################

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


## Also need to make sure that month 12 becomes month 1 of the next age group

#### Mortality ####

mortality.func <- function (Age, Nat.Mort, Effort, Max.Cell, Month, Year, Population, Select){
  
  tot.survived <- sapply(seq(Max.Cell), function(Cell) {
    tot.survived <- Population[Cell,Month,Age]*exp(-Select[(((Age-1)*12)+1),Year]*Effort[Cell,Month,Year])*exp(-(Nat.Mort/12)) 
    
    #Calculate Catch
    sel.effort <- (Select[(((Age-1)*12)+1),1]*Effort[Cell,Month,Year])
    catch <- sel.effort/((Nat.Mort/12)+sel.effort)
    catch <- tot.survived[Cell]*catch
  })
  
  # Calculate catch
  
  list.surv.catch <- list(tot.survived, catch)
  return(list.surv.catch)
  
}


#### Movement ####

movement.func <- function (Age, Month, Population, Max.Cell, Adult.Move){

  All.Movers <- NULL

  ## Adult Movement

  Moving <- sapply(1:30, function(a){
    All.Movers <- sapply(seq(Max.Cell), function(Cell){
      temp <- matrix(Population[ ,Month, a])
      Pop2 <- matrix(Adult.Move[Cell, ] * temp[Cell,1]) %>%  #This should give you the number of fish that move from sites to all the other sites
        unlist()
      
    })
    Movement <- rowSums(All.Movers)
  })


#dim(Moving) <- c(Max.Cell,1,30)
  # #End bracket for movement in each cell
  # 
  # Movement2 <- colSums(Pop2)
  # Moved <- sum(Movement2)
  # # print(isTRUE(all.equal(Pop2, Moved)))
  # # if((isTRUE(all.equal(Pop2, Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
  # # {print(Pop2)
  # #   print(Moved)} else{ }
  # #
  # All.Movers <- cbind(All.Movers, Movement2)

  return(Moving)
}

#### Recruitment Function ####

## Recruitment

recruitment.func <- function(Population, settlement, Max.Cell, BHa, BHb, Mature, Weight, PF){
  adults <- Population[ ,10, ] %>% 
    colSums(.) # Gives us just females because they are the limiting factor for reproduction
  adults <- adults * PF
  
  tot.recs <- lapply(1:dim(Population)[3], function(Age){
    SB <- adults[Age] * Mature[Age,1] * Weight[(Age*12)+1] #Gives us spawning biomass in each age group at the end of the year, hence the x 12+1 as it starts at 1 not zero
    TotMatBio <- sum(SB) #Gives us total mature spawning biomass
    recs <- (SB/BHa+BHb*SB) # This gives us males and females to go into the next generation
  })
  tot.recs <- do.call(rbind, tot.recs)

  settle.recs <- sum(tot.recs)
  
  settled <- settle.recs*settlement[,2]
  return(settled)
}

#### MSY Functions ####
## Need to get yield per recruit, biomass and mature/spawning biomass 
# We need number survived in each age class, fishing mort. and total mort. and weight at age
ypr.func <- function(Survived, Fishing.Mort, Nat.Mort, Weight, Age, Select){
  
  Survived.age <- Survived[Age,12]
  
  fishing.mort <- Fishing.Mort[FM]*Select[(((Age-1)*12)+1),3]
  total.mort <- fishing.mort+Nat.Mort
  
  ypr <- Survived.age*(fishing.mort/total.mort)*(1-exp(-total.mort))*Weight[(Age*12)+1]
  
  return(ypr)
  
}

bio.func <- function(Survived, Weight, Population){
  
  for (A in 1:30){
    adults <- Survived[ ,12]
    TotBio<- lapply(1:dim(Population)[3], function(Age){
      SB <- adults[Age] * Weight[(Age*12)+1] #Gives us spawning biomass in each age group at the end of the year, hence the x 12+1 as it starts at 1 not zero
      Bio <- sum(SB) #Gives us total biomass
    })
    TotalBio <- do.call(rbind, TotBio)
  }
  return(TotalBio)
}
  
SB.func <- function(Survived, Weight, Population, Maturity){
  for (A in 1:30){
    adults <- Survived[ ,10]
    adults <- adults * 0.5
    
    MatBio<- lapply(1:dim(Population)[3], function(Age){
      SB <- adults[Age] * Maturity[Age,1] * weight[(Age*12)+1] #Gives us spawning biomass in each age group at the end of the year, hence the x 12+1 as it starts at 1 not zero
      TotMatBio <- sum(SB) #Gives us total mature spawning biomass
    })
    MatBio <- do.call(rbind, MatBio)
  }
  return(MatBio)
}

#### Plotting ####

total.plot.func <- function (pop) {
  TimeSeries <- ggplot(pop)+
    geom_line(aes(x=Year, y=Tot.Pop)) +
    scale_x_continuous("Year", breaks = c(1960, 1970, 1980, 1990, 2000, 2010, 2020))+
    xlab("Year")+
    ylab("Total Population")+
    theme_classic()+
    geom_vline(xintercept=1987, linetype="dashed", color="grey20")+
    geom_vline(xintercept=2005, colour="grey20")+
    geom_vline(xintercept=2017, linetype="dashed", colour="grey20")
    
    return(TimeSeries)
}

spatial.plot.func <- function (area, pop, pop.breaks, pop.labels, colours){
  
  water <- area %>%
    mutate(pop = round(pop, digits=0)) %>% 
    mutate(pop_level = cut(pop, pop.breaks, include.lowest=T)) 
  
  nb.cols <- length(pop.breaks)
  mycols <- colorRampPalette(rev(brewer.pal(8, colours)))(nb.cols)
  
  map <- ggplot(water)+
    geom_sf(aes(fill=pop_level), lwd=0)+
    scale_fill_manual(name="Population", values= mycols, drop=FALSE)+
    #scale_color_manual(name="Fished", values=c("black", "black"))+
    theme_void()
  
  return(map)
  
}

age.plot.func <- function (pop, NTZs){
  
  if(YEAR >=27 & YEAR <= 45){
    NoTakeAges <- pop[c(NTZs[[1]]),12, ]
    FishedAges <- pop[-c(NTZs[[1]]),12, ]
    
  } else if(YEAR>45 & YEAR<=57){
    NoTakeAges <- pop[c(NTZs[[2]]),12, ]
    FishedAges <- pop[-c(NTZs[[2]]),12, ]
    
  } else if (YEAR >57){
    NoTakeAges <- pop[c(NTZs[[3]]),12, ]
    FishedAges <- pop[-c(NTZs[[3]]),12, ]
    
  } else if (YEAR<27){
    NoTakeAges <- as.data.frame(array(0, dim=c(100, 30)))
    FishedAges <- pop[ ,12, ]
  }
  
  NoTakeAges <- as.data.frame(colSums(NoTakeAges)) %>% 
    rename(Number = "colSums(NoTakeAges)") %>% 
    mutate(Status = "NTZ") %>% 
    mutate(Age = seq(1:30))
    
    
  FishedAges <- as.data.frame(colSums(FishedAges)) %>% 
    rename(Number = "colSums(FishedAges)") %>% 
    mutate(Status = "Fished") %>% 
    mutate(Age = seq(1:30))
  
  AllAges <- rbind(NoTakeAges, FishedAges) %>% 
    mutate(Age = as.factor(Age)) %>%
    group_by(Status) %>% 
    mutate(TotalZone=sum(Number)) %>% 
    ungroup() %>% 
    mutate(Proportion = (Number/TotalZone))
  
barplot.ages <- ggplot(AllAges)+
  geom_bar(aes(y=Proportion, x=Age, fill=Status), position="dodge", stat="identity")+
  theme_classic()

return(barplot.ages)
   
}

msy.plot.func <- function(yield, biomass,spawning){
  
  yield.plot <- ggplot()+
    geom_line(aes(x=yield[,1], y=yield[,2]))+
    theme_classic()+
    xlab("Fishing Mortality")+
    ylab("Yeild per Recruit")
  
  biomass.plot <- ggplot()+
    geom_line(aes(x=biomass[,1], y=biomass[,2]))+
    theme_classic()+
    xlab("Fishing Mortality")+
    ylab("Total Biomass")
  
  spawning.plot <- ggplot()+
    geom_line(aes(x=spawning[,1], y=spawning[,2]))+
    theme_classic()+
    xlab("Fishing Mortality")+
    ylab("Total Spawning Biomass")
  
  msy.plots <- list(yield.plot, biomass.plot, spawning.plot)
  return(msy.plots)
  
}
  


  
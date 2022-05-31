###################################################

# Script that has all the functions
# Functions are only needed when running the model
# Have functions for movement, mortality, recruitment
# and for plotting

###################################################

## Also need to make sure that month 12 becomes month 1 of the next age group

#### Mortality ####

mortality.func <- function (Age, mort.50, mort.95, Nat.Mort, NTZ, Effort, Max.Cell, Month, Year, Population, Select){
  
  if(BurnIn){
    
    tot.survived <- Population[ , Month, Age]*exp(-(Nat.Mort/12))

  } else {
    tot.survived <- sapply(seq(Max.Cell), function(Cell) {
      tot.survived <- Population[Cell,Month,Age]*exp(-Select[(((Age-1)*12)+1)]*Effort[Cell,Month,Year])*exp(-(Nat.Mort/12)) 
    })
  }
  return(tot.survived)
}

#### Movement ####

movement.func <- function (Age, Month, Population, Max.Cell, Adult.Move, Juv.Move){
  
  All.Movers <- NULL
  
  ## Juvenile Movement
  if(Age<=4){
    
  Juv.Pop <- matrix(Population[ , Month, Age]) # This gives you the fish in all the sites at time step Month-1 of age A-1
  
  Juv.Pop2 <- sapply(seq(Max.Cell), function(Cell){
    Juv.Pop2 <- as.matrix(Juv.Move[Cell, ] * Juv.Pop[Cell,1]) #This should give you the number of fish that move from each cell to all the other sites
    
  })
  
  Juv.Movement2 <- rowSums(Juv.Pop2)
  Juv.Moved <- sum(Juv.Movement2)
  # print(isTRUE(all.equal(Juv.Pop2, Juv.Moved)))
  # if((isTRUE(all.equal(Juv.Pop2, Juv.Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
  # {print(Juv.Pop2)
  #   print(Juv.Moved)} else{ } # Can Put This Check Back in Periodically to Confirm model Functioning
  All.Movers <- cbind(All.Movers, Juv.Movement2)
  
  return(All.Movers)
  
  } else {  
    
    ## Adult Movement
    Pop <- matrix(Population[ , Month, Age]) #This gives you the fish in all the sites at timestep t-1 of age A-1
    
    Pop2 <- sapply(seq(Max.Cell), function(Cell){
      Pop2 <- as.matrix(Adult.Move[Cell, ] * Pop[Cell,1]) #This should give you the number of fish that move from sites to all the other sites
      
    }) #End bracket for movement in each cell
    
    Movement2 <- rowSums(Pop2)
    Moved <- sum(Movement2)
    # print(isTRUE(all.equal(Pop2, Moved)))
    # if((isTRUE(all.equal(Pop2, Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
    # {print(Pop2)
    #   print(Moved)} else{ }
    
    All.Movers <- cbind(All.Movers, Movement2)
    
    return(All.Movers)

  } #End adult movement 
  
}

#### Recruitment Function ####

## Recruitment

recruitment.func <- function(Population, mat.95, mat.50, settlement, Max.Cell, BHa, BHb, Mature, Weight, PF){
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

# Then get added to January Population 
#### Plotting ####

total.plot.func <- function (pop) {
    TimeSeries <- ggplot()+
      geom_line(aes(x=1:nrow(pop), y=pop))+
      xlab("Year")+
      ylab("Total Population")+
      theme_classic()
    
    return(TimeSeries)
}

spatial.plot.func <- function (area, pop, pop.breaks, pop.labels, colours){
  
  water <- water%>%
    mutate(pop = round(pop, digits=0)) %>% 
    mutate(pop_level = cut(pop, pop.breaks, include.lowest=T)) 
  
  nb.cols <- length(pop.breaks)
  mycols <- mycolors <- colorRampPalette(brewer.pal(8, "RdBu"))(nb.cols)
  
  map <- ggplot(water)+
    geom_sf(aes(fill=pop_level, color=Fished))+
    scale_fill_manual(name="Population", values= mycols, drop=FALSE)+
    scale_color_manual(name="Fished", values=c("white", "black"))+
    theme_void()
  
  return(map)
  
}

age.plot.func <- function (pop, NTZs){
  
  if(YEAR >=27 & YEAR <= 45){
    NoTakeAges <- pop[c(NTZs[[1]]),12, ]
    FishedAges <- pop[-c(NTZs[[1]]),12, ]
    
  } else if(YEAR>45 & YEAR<=53){
    NoTakeAges <- pop[c(NTZs[[2]]),12, ]
    FishedAges <- pop[-c(NTZs[[2]]),12, ]
    
  } else if (YEAR >53){
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

  
  


  
  
###################################################

# Script that has all the functions
# Functions are only needed when running the model
# Have functions for movement, mortality, recruitment
# and for plotting

###################################################

## Also need to make sure that month 12 becomes month 1 of the next age group

#### Mortality ####

mortality.func <- function (Age, mort.50, mort.95, Nat.Mort, NTZ, Effort, Cell, Max.Cell, Month, Year, Population, Select){
  
  if(BurnIn){
   
    tot.survived <- Population[, Month-1,Age]*exp(-Nat.Mort)
    
  }
  #sa <- 1/(1+(exp(-log(19)*((Age-mort.50/mort.95-mort.50))))) # This puts in size selectivity - GOT A50 IN WRONG PLACE, FIXED IT NOW
  
  tot.survived <- sapply(seq(Max.Cell), function(Cell) {
    survived <- Population[Cell,Month-1,Age]*exp(-Select[((Age*12)+1),1]*Effort[Cell,Month-1,Year])*exp(-Nat.Mort) # Can't do selectivity like this - has to be specific for each month as it changes! 
    survived
  })
  
  return(tot.survived)
}

#### Movement ####

movement.func <- function (Age, Month, Population, Max.Cell, Adult.Move, Juv.Move){
  
  All.Movers <- NULL
  
  ## Juvenile Movement
  if(Age<=4){
    
  Juv.Pop <- matrix(Population[ , Month-1, Age-1]) # This gives you the fish in all the sites at time step Month-1 of age A-1
  
  Juv.Pop2 <- sapply(seq(Max.Cell), function(Cell){
    Juv.Pop2 <- as.matrix(Juv.Move[Cell, ] * Juv.Pop[Cell,1]) #This should give you the number of fish that move from each cell to all the other sites
  })
  
  Juv.Movement2 <- rowSums(Juv.Pop2)
  # Juv.Moved <- sum(Juv.Movement2)
  # print(isTRUE(all.equal(Juv.Movers, Juv.Moved)))
  # if((isTRUE(all.equal(Juv.Movers, Juv.Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
  # {print(Juv.Movers)
  #   print(Juv.Moved)} else{ } Can Put This Check Back in Periodically to Confirm model Functioning
  
  All.Movers <- cbind(All.Movers, Juv.Movement2)
  
  return(All.Movers)
  
  
  } else {  
    
    ## Adult Movement
    Pop <- matrix(Population[ , Month-1, Age-1]) #This gives you the fish in all the sites at timestep t-1 of age A-1
    
    Pop2 <- sapply(seq(Max.Cell), function(Cell){
      Pop2 <- as.matrix(Adult.Move[Cell, ] * Pop[Cell,1]) #This should give you the number of fish that move from site s to all the other sites
    }
    ) #End bracket for movement in each cell
    
    Movement2 <- rowSums(Pop2)
    # Moved <- sum(Movement2)
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

recruitment.func <- function(Population, Age, mat.95, mat.50, settlement, Max.Cell, BHa, BHb, Mature, Weight, PF){
  adults <- Population[ ,10, ]*PF # Gives us just females because they are the limiting factor for reproduction
  tot.recs <- data.frame(matrix(0, nrow = Max.Cell, ncol = dim(Population)[3]))
  
  tot.recs <- lapply(1:dim(Population)[3], function(Age){
    SB <- adults[ , Age] * Mature[Age,1] * Weight[(Age*12)+1] #Gives us spawning biomass in each age group at the end of the year, hence the x 12+1 as it starts at 1 not zero
    TotMatBio <- sum(SB) #Gives us total mature spawning biomass
    recs <- (SB/BHa+BHb*SB) # This gives us males and females to go into the next generation
  })
  tot.recs <- do.call(rbind, tot.recs)

  settle.recs <- sum(tot.recs)
  
  settled <- tot.recs*settlement[,2]
  return(settled)
}

# Then get added to January Population 
#### Plotting ####

plotting.func <- function (area, pop ,pop.breaks, pop.labels, colours) {
  
  if(PlotTotal){
    
    TimeSeries <- ggplot()+
      geom_line(aes(x=1:nrow(Total), y=Total))+
      xlab("Year")+
      ylab("Total Population")+
      theme_classic()
    
    return(TimeSeries)

  } else {
    
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

}


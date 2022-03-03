###################################################

# Script that has all the functions
# Functions are only needed when running the model
# Have functions for movement and mortality 

###################################################

##NEED TO GET ALL OF THESE TO RETURN SOMETHING SO THAT IT ACUTALLY DOES WHAT YOU WANT IT TO
## Also need to make sure that month 12 becomes month 1 of the next age group

#### Mortality ####

mortality.func <- function (Age, mort.50, mort.95, Nat.Mort, NTZ, Effort, Cell, Max.Cell, Month, Year, Population){
  
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
  
  if(Age==2 & Month==12){Population[, 12, Age-1] <- matrix(floor(runif(Max.Cell, 1, 1000)))
  } else{ }
  
  ## Juvenile Movement
  if(Age<=4){
    
  Juv.Pop <- matrix(Population[ , Month-1, Age-1]) # This gives you the fish in all the sites at time step Month-1 of age A-1
  Juv.Movers <- sum(Juv.Pop)
  Juv.Moving <- NULL
  Juv.Moved <- NULL
  
  for(Cell in 1:Max.Cell){
    Juv.Pop2 <- as.matrix(Juv.Move[Cell, ] * Juv.Pop[Cell,1]) #This should give you the number of fish that move from each cell to all the other sites
    Juv.Pop2 <- t(Juv.Pop2) # Gives you a row with each cell being the number of fish that moved from cell 1 to all other cells 
    Juv.Moving <- rbind(Juv.Moving, Juv.Pop2) #This should give you an array with rows representing the fish that move from each site to all the other sites 
  } #End bracket for movement in each cell
  
  Juv.Movement2 <- as.matrix(colSums(Juv.Moving))
  Juv.Moved <- sum(Juv.Movement2)
  print(isTRUE(all.equal(Juv.Movers, Juv.Moved)))
  if((isTRUE(all.equal(Juv.Movers, Juv.Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
  {print(Juv.Movers)
    print(Juv.Moved)} else{ }
  
  All.Movers <- cbind(All.Movers, Juv.Movement2)
  
  return(All.Movers)
  
  
  } else {  
    
    ## Adult Movement
    Pop <- matrix(Population[ , Month-1, Age-1]) #This gives you the fish in all the sites at timestep t-1 of age A-1
    Movers <- sum(Pop)
    Moving <- NULL
    Moved <- NULL
    
    for(Cell in 1:Max.Cell){
      Pop2 <- as.matrix(Adult.Move[Cell, ] * Pop[Cell,1]) #This should give you the number of fish that move from site s to all the other sites
      Pop2 <- t(Pop2)
      Moving <- rbind(Moving, Pop2) #This should give you an array with 143 rows representing the fish that move from each site to all the other sites 
    } #End bracket for movement in each cell
    
    Movement2 <- as.matrix(colSums(Moving))
    Moved <- sum(Movement2)
    print(isTRUE(all.equal(Movers, Moved)))
    if((isTRUE(all.equal(Movers, Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
    {print(Movers)
      print(Moved)} else{ }
    
    All.Movers <- cbind(All.Movers, Movement2)
    
    return(All.Movers)

  } #End adult movement 
  
}

#### Plotting ####

plotting.func <- function (area, n.breaks, colours) {
  
  water <- water%>%
    mutate(pop = round(pop, digits=0)) %>% 
    mutate(pop_level = cut_interval(pop, n=n.breaks)) 
  
  cols <- brewer.pal(n.breaks, colours)
  
  map <- ggplot(water)+
    geom_sf(aes(fill=pop_level, color=Fished))+
    scale_fill_manual(name="Population", values= cols, drop=FALSE)+
    scale_color_manual(name="Fished", values=c("white", "black"))+
    theme_void()
  
  return(map)
}

## NEED TO SET THE GROUPS TO BE THE SAME EVERY TIME SO THAT COlOURS CHANGE ON THE MAP PROPERLY
  










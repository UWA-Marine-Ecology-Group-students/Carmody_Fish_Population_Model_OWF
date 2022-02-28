###################################################

# Script that has all the functions
# Functions are only needed when running the model
# Have functions for movement and mortality 

###################################################

#### Mortality ####

mortality.func <- function (Age, mort.50, mort.95, Nat.Mort, NTZ, Effort, Cell, Max.Cell, Month, Year, Population){
  
  sa <- 1/(1+(exp(-log(19)*((Age-mort.95/mort.95-mort.50))))) # This puts in size selectivity

for(Cell in 1:Max.Cell){
  
  # Up until 1987 when there were no sanctuaries
  if(Year>0&Year<=26){Population[Cell,Month,Age] <- Population[Cell,Month,Age]*exp(-sa*Effort[Cell,Month,Year])*(1-Nat.Mort)}
  
  # From 1987 to 2005 when the first sanctuaries were in place
  else if(Year>=27&Year<=44){
    # If it's a NTZ then we don't have any fishing mortality in that cell
    if(NTZ[Cell, 3]=="Y"){Population[Cell,Month,Age] <- YearlyTotal[Cell,Month,Age]*exp(-sa*Effort[Cell,Month,Year])*(1-Nat.Mort)}
    else { }  
  }
  
  # From 2005 to 2017 when we had all the state sanctuaries but no commonwealth one
  else if(Year>=45&Year<=57){
    # If it's a NTZ then we don't have any fishing mortality in that cell
    if(NTZ[Cell, 4]=="Y"){Population[Cell,Month,Age] <- Population[Cell,Month,Age]*exp(-sa*Effort[Cell,Month,Year])*(1-Nat.Mort)}
    else { }  
  }
  
  # From 2017 onwards when we had all the sanctuaries 
  else if(Year>57){
    # If it's a NTZ then we don't have any fishing mortality in that cell
    if(NTZ[Cell, 5]=="Y"){Population[Cell,Month,Age] <- Population[Cell,Month,Age]*exp(-sa*Effort[Cell,Month,Year])*(1-Nat.Mort)}
    else { }  
  }
} # End fishing mortality for each cell

}

#### Movement ####

movement.func <- function (Age, Month, Population, Max.Cell, Adult.Move, Juv.Move){
  
  if(Age==2 & Month==12){Population[, 12, Age-1] <- matrix(floor(runif(Max.Cell, 1, 2000)))
  } else{ }
  
  ## Juvenile Movement
  if(Age<=4){
    
  Juv.Pop <- matrix(Population[ , Month-1, Age-1]) # This gives you the fish in all the sites at time step Month-1 of age A-1
  Juv.Movers <- sum(Juv.Pop)
  Juv.Movement <- NULL
  Juv.Movement2 <- NULL
  
  for(Cell in 1:Max.Cell){
    Juv.Pop2 <- as.matrix(Juv.Move[Cell, ] * Juv.Pop[Cell,1]) #This should give you the number of fish that move from each cell to all the other sites
    Juv.Pop2 <- t(Juv.Pop2) # Gives you a row with each cell being the number of fish that moved from cell 1 to all other cells 
    Juv.Movement <- rbind(Juv.Movement, Juv.Pop2) #This should give you an array with rows representing the fish that move from each site to all the other sites 
  } #End bracket for movement in each cell
  
  Juv.Movement2 <- as.matrix(colSums(Juv.Movement))
  Juv.Moved <- sum(Juv.Movement2)
  print(isTRUE(all.equal(Juv.Movers, Juv.Moved)))
  if((isTRUE(all.equal(Juv.Movers, Juv.Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
  {print(Juv.Movers)
    print(Juv.Moved)} else{ }
  
  Population[ , Month, Age-1] <- Juv.Movement2 # End Juvenile Movement
  } else {  
    
    ## Adult Movement
    Pop <- matrix(Population[ , Month-1, Age-1]) #This gives you the fish in all the sites at timestep t-1 of age A-1
    Movers <- sum(Pop)
    Movement <- NULL
    Movement2 <- NULL
    
    for(Cell in 1:Max.Cell){
      Pop2 <- as.matrix(Adult.Move[Cell, ] * Pop[Cell,1]) #This should give you the number of fish that move from site s to all the other sites
      Pop2 <- t(Pop2)
      Movement <- rbind(Movement, Pop2) #This should give you an array with 143 rows representing the fish that move from each site to all the other sites 
    } #End bracket for movement in each cell
    
    Movement2 <- as.matrix(colSums(Movement))
    Moved <- sum(Movement2)
    print(isTRUE(all.equal(Movers, Moved)))
    if((isTRUE(all.equal(Movers, Moved))) == FALSE) #This just prints the values of the fish that moved if it's not the same as the fish that were meant to move
    {print(Movers)
      print(Moved)} else{ }
    
    Population[ , Month, Age-1] <- Movement2
    
  } #End adult movement 
  
}


  










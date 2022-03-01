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
  
  tot.survived <- array(0, dim=c(Max.Cell,12,30))

for(Cell in 1:Max.Cell){
  
  # Up until 1987 when there were no sanctuaries
  if(Year>0&Year<=26){survived <- Population[Cell,Month-1,Age]*exp(-sa*Effort[Cell,Month-1,Year])*(1-Nat.Mort)
  tot.survived[Cell, Month-1, Age] <- survived}
  
  # From 1987 to 2005 when the first sanctuaries were in place
  else if(Year>=27&Year<=44){
    # If it's a NTZ then we don't have any fishing mortality in that cell
    if(NTZ[Cell, 3]=="Y"){survived <- YearlyTotal[Cell,Month-1,Age]*exp(-sa*Effort[Cell,Month-1,Year])*(1-Nat.Mort)
    tot.survived[Cell, Month-1, Age] <- survived}
    else { }  
  }
  
  # From 2005 to 2017 when we had all the state sanctuaries but no commonwealth one
  else if(Year>=45&Year<=57){
    # If it's a NTZ then we don't have any fishing mortality in that cell
    if(NTZ[Cell, 4]=="Y"){survived <- Population[Cell,Month-1,Age]*exp(-sa*Effort[Cell,Month-1,Year])*(1-Nat.Mort)
    tot.survived[Cell, Month-1, Age] <- survived}
    else { }  
  }
  
  # From 2017 onwards when we had all the sanctuaries 
  else if(Year>57){
    # If it's a NTZ then we don't have any fishing mortality in that cell
    if(NTZ[Cell, 5]=="Y"){survived <- Population[Cell,Month-1,Age]*exp(-sa*Effort[Cell,Month-1,Year])*(1-Nat.Mort)
    tot.survived[Cell, Month-1, Age] <- survived}
    else { }  
  }
  

} # End fishing mortality for each cell
  
  return(tot.survived)

}

#### Movement ####

movement.func <- function (Age, Month, Population, Max.Cell, Adult.Move, Juv.Move){
  
  All.Movers <- NULL
  
  if(Age==2 & Month==12){Population[, 12, Age-1] <- matrix(floor(runif(Max.Cell, 1, 2000)))
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
  # 
  # if(Month==13){
  #   eval.parent(substitute(Population[ , 12, Age-1] <- Juv.Movement2)) # Deals with issue of there not being a 13th column
  # } else {
  #   eval.parent(substitute(Population[ , Month, Age-1] <- Juv.Movement2))} # End Juvenile Movement
  
  
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
    
    # if(Month==13){
    #   eval.parent(substitute(Population[ , 12, Age-1] <- Movement2)) # Deals with issue of there not being a 13th column
    # } else {
    #   eval.parent(substitute(Population[ , Month, Age-1] <- Movement2))} # End Juvenile Movement
    # 
  } #End adult movement 
  
}


  










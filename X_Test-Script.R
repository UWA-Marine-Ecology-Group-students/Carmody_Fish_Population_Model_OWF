###################################################

# Script just for fiddling with the model and trying 
# out new things
# Need to make a smaller grid at some point

###################################################
YearlyTotal <- array(0, dim = c(NCELL,12,8)) #This is our yearly population split by age category (every layer is an age group)
for(d in 1:8){
  YearlyTotal[,1,d] <- matrix(floor(runif(NCELL, 1, 2000)))
}

PopTotal <- array(0, dim=c(NCELL, 12, 59)) # This is our total population, all ages are summed and each column is a month (each layer is a year)

for(YEAR in 2:2){
  
  for(MONTH in 2:12){
    
    ## Movement - this is where you'd change things to suit the months
    for(A in 2:dim(YearlyTotal)[3]){
      
      movement.func(Age=A, Month=MONTH, Population=YearlyTotal, Max.Cell=NCELL, Adult.Move= movement, 
                    Juv.Move=juv_movement)

    }  # End bracket for age classes

  } #End bracket for months
  
    # ## Recruitment
    # Recs <- matrix(0, nrow=1, ncol=1) #Blank matrix to add the recruits to
    # Recruitment <- as.matrix(colSums(Pop[,t,], dims=1)) 
    # for(A in 1:dim(Recruitment)[1]){
    #   Mature <- 1/(1+(exp(-log(19)*((A-M95)/(M95-M50))))) #Number of mature individuals in each age class
    #   S <- colSums(Recruitment)*Mature #Spawning stock
    #   Rec <- a*S*exp(-b*S) #Number of recruits from that age class
    #   Recs <- rbind(Recs, Rec) #Combine recruits from each age class into one dataframe
    # }
    
    #PopTotal[ , , YEAR] <- rowSums(YearlyTotal) # This flattens the matrix to give you the number of fish present in the population each month, with layers representing the years
     PopTotal[ , , YEAR] <- YearlyTotal[ , , 1] 
    
  
  print(YEAR)
  water$pop <- PopTotal[ , 12, YEAR] # We just want the population at the end of the year
  
  water <- water%>%
    mutate(Population = ifelse(pop < 1000, "<1000",
                               ifelse (pop>1000 & pop<110, "1000-1100",
                                       ifelse (pop>1100 & pop<1200, "1100-1200",
                                               ifelse (pop>1200 & pop<1300, "1200-1300",
                                                       ifelse(pop>1300 & pop<1400, "1300-1400",
                                                              ifelse(pop>1400 & pop<1500, "1400-1500",
                                                                     ifelse(pop>1500 & pop<1600, "1500-1600", ">1600"
                                                                            
                                                                     ))))))))%>%
    mutate(Population = factor(Population))
  
  water$Population <- fct_relevel(water$Population, "<1000",  "1000-1100", "1100-1200", "1200-1300",
                                  "1300-1400", "1400-1500", "1500-1600", ">1600")
  
  print(map <- ggplot(water)+
          geom_sf(aes(fill=Population, color=Fished))+
          scale_fill_manual(name="Population", values= cols, drop=FALSE)+
          scale_color_manual(name="Fished", values=c("white", "black"))+
          theme_void())
  Sys.sleep(3)
}

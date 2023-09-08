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
  
  Survived.age <- Survived[Age]
  
  fishing.mort <- Fishing.Mort[FM]*Select[Age,12,1]
  total.mort <- fishing.mort+Nat.Mort
  
  ypr <- Survived.age*(fishing.mort/total.mort)*(1-exp(-total.mort))*Weight[Age,12]
  
  return(ypr)
  
}

bio.func <- function(Survived, Weight, Population){
    adults <- Survived
      SB <- adults * Weight[ ,12] 
      TotalBio <- sum(SB) #Gives us total biomass
  return(TotalBio)
}
  
SB.func <- function(Survived, Weight, Population, Maturity){
  for (A in 1:30){
    adults <- Survived
    adults <- adults * 0.5
    
    MatBio<- lapply(1:dim(Population)[3], function(Age){
      SB <- adults[Age] * Maturity[Age,12] * Weight[Age,12] #Gives us spawning biomass in each age group at the end of the year, hence the x 12+1 as it starts at 1 not zero
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
  
#### FORMATTING FUNCTIONS ####
#* Fish by age and by zone ####


zone.fish.age <- function(NTZ.Ages, F.Ages){
  
  NTZ_Ages <- array(0, dim=c(30,59,100))
  F_Ages <- array(0, dim=c(30,59,100))
  
  for(S in 1:4){
    
    SP_Pop_NTZ <- NTZ.Ages[[S]]
    SP_Pop_F <- F.Ages[[S]]
    
    for(SIM in 1:length(SP_Pop_NTZ)){
      
      temp <- as.data.frame(colSums(SP_Pop_NTZ[[SIM]])) %>% 
        mutate(Age = seq(1:30)) %>% 
        pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
        mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", ""))) %>% 
        mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
        group_by(Age, Mod_Year) %>% 
        summarise(across(where(is.numeric), sum)) %>% 
        ungroup() %>%
        dplyr::select(!Num_Year) %>%
        pivot_wider(names_from="Mod_Year",id_col="Age", values_from="Number",values_fn = list(count=list)) %>% 
        dplyr::select(!Age) %>% 
        unlist()
      
      temp2 <- array(temp, dim=c(30,59))
      
      NTZ_Ages[,,SIM] <- temp2
      
      temp <- as.data.frame(colSums(SP_Pop_F[[SIM]])) %>% 
        mutate(Age = seq(1:30)) %>% 
        pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
        mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", ""))) %>% 
        mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
        group_by(Age, Mod_Year) %>% 
        summarise(across(where(is.numeric), sum)) %>% 
        ungroup() %>% 
        pivot_wider(names_from="Mod_Year",id_col="Age", values_from="Number") %>% 
        dplyr::select(!Age) %>% 
        unlist()
      
      temp2 <- array(temp, dim=c(30,59))
      
      F_Ages[,,SIM] <- temp2
    }
    
    Age.Dist.NTZ[[S]] <- NTZ_Ages
    Age.Dist.F[[S]] <- F_Ages
    
  }
  
  Age.Dist[[1]] <- Age.Dist.NTZ
  Age.Dist[[2]] <- Age.Dist.F
  
  return(Age.Dist)
} 

#* Catch by age ####

age.catch <- function(Catches, NTZ.Cells, F.Cells){
  Age.Catch <- list()
  temp4 <- list()
  temp5 <- list()
  catch <- array(0, dim=c(NCELL,59,100))
  
  for (S in 1:4){
    
    temp <- Catches[[S]] # NTZ and F are the same files 
    
    for(SIM in 1:100){
      catch[,,SIM] <- temp[[SIM]]
    }
    
    temp2 <- catch[as.numeric(NTZ.Cells),,]
    temp3 <- catch[as.numeric(F.Cells),,]
    
    temp4[[S]] <- temp2
    
    temp5[[S]] <- temp3
    
  }
  
  Age.Catch[[1]] <- temp4
  Age.Catch[[2]] <- temp5
  
  return(Age.Catch)
  
}


#* Means  ####

means.func <- function(Age, Catch, Zones, Age.NTZ, Age.F, Catch.NTZ, Catch.F){
  Ages.All <- list()
  Ages.F <- list()
  Ages.NTZ <- list()
  
  Catches.All <- list()
  Catches.NTZ <- list()
  Catches.F <- list()
  
  catch <- array(0, dim=c(30,59,100))
  catch.NTZ <- array(0, dim=c(30,59,100))
  catch.F <- array(0, dim=c(30,59,100))
  
  for(S in 1:4){
    
    #*************AGES***************
    temp.All <- Age[[S]]
    
    temp2 <- temp.All %>% 
      rowMeans(., dim=2) %>% 
      as.data.frame(.) 
    
    Ages.All[[S]] <- temp2

    
    if(Zones == TRUE){
      # Ages NTZ
      temp.NTZ <- Age.NTZ[[S]]
      
      temp2 <- temp.NTZ %>% 
        rowMeans(., dim=2) %>% 
        as.data.frame(.) 
      
      Ages.NTZ[[S]] <- temp2
      
      # Ages F
      temp.F <- Age.F[[S]]
      
      temp2 <- temp.F %>% 
        rowMeans(., dim=2) %>% 
        as.data.frame(.) 
      
      Ages.F[[S]] <- temp2
    } else {}
    
    
    #*********************************
    #*************CATCHES*************
    ## All
    temp.c <- Catch[[S]]
    for(SIM in 1:100){
      
      catch[,,SIM] <- temp.c[[SIM]]
      
    }
    
    temp2.c <- catch %>% 
      rowMeans(., dim=2) %>% 
      as.data.frame(.) 
    
    Catches.All[[S]] <- temp2.c
    
    if(Zones == TRUE){
      ## NTZ
      temp.c.NTZ <- Catch.NTZ[[S]]
      for(SIM in 1:100){
        
        catch.NTZ[,,SIM] <- temp.c.NTZ[[SIM]]
        
      }
      
      temp2.c.NTZ <- catch.NTZ %>% 
        rowMeans(., dim=2) %>% 
        as.data.frame(.) 
      
      if(S==1|S==3){
        temp2.c.NTZ[,] <- 0
      } else { }
      
      Catches.NTZ[[S]] <- temp2.c.NTZ
      
      ## F
      temp.c.F <- Catch.F[[S]]
      for(SIM in 1:100){
        
        catch.F[,,SIM] <- temp.c.F[[SIM]]
        
      }
      
      temp2.c.F <- catch.F %>% 
        rowMeans(., dim=2) %>% 
        as.data.frame(.) 
      
      Catches.F[[S]] <- temp2.c.F
    } else { }
    
    
    #*********************************  
  }
  Age.and.Catch <- list()
  
  Age.and.Catch[[1]] <- Ages.All
  Age.and.Catch[[2]] <- Ages.F
  Age.and.Catch[[3]] <- Ages.NTZ
  
  Age.and.Catch[[4]] <- Catches.All
  Age.and.Catch[[5]] <- Catches.NTZ
  Age.and.Catch[[6]] <- Catches.F
  
  return(Age.and.Catch)
  
}

#* Biomass function #####

biomass.func <- function(Age, Zones, Age.NTZ, Age.F, Weights){
  Biomass.All <- list()
  Biomass.NTZ <- list()
  Biomass.F <- list()
  
  for(S in 1:4){
    ## All
    temp.All <- Age[[S]]
    
    temp3 <- temp.All %>% 
      rowMeans(., dim=2) %>% 
      as.data.frame(.) %>% 
      mutate_all(.,function(col){Weights$V12*col}) 
    
    Biomass.All[[S]] <- temp3
    
    if(Zones == TRUE){
      ## NTZ
      temp.NTZ <- Age.NTZ[[S]]
      
      temp3.NTZ <- temp.NTZ %>% 
        rowMeans(., dim=2) %>% 
        as.data.frame(.) %>% 
        mutate_all(.,function(col){Weights$V12*col}) 
      
      Biomass.NTZ[[S]] <- temp3.NTZ
      
      ## F
      temp.F <- Age.F[[S]]
      
      temp3.F<- temp.F %>% 
        rowMeans(., dim=2) %>% 
        as.data.frame(.) %>% 
        mutate_all(.,function(col){Weights$V12*col}) 
      
      Biomass.F[[S]] <- temp3.F
    } else {  }
  }
 
  
  Biomass <- list()
  Biomass[[1]] <- Biomass.All
  Biomass[[2]] <- Biomass.NTZ
  Biomass[[3]] <- Biomass.F
  
  return(Biomass)
  
}

#### Data formatting ####

## Whole population

total.pop.format <- function(pop.file.list, scenario.names, nsim, nyears, startyear, maxage, mat, kg){
  
  total_pop <- NULL
  
  for(i in 1:length(pop.file.list)){
    temp <- pop.file.list[[i]]
    
    median.ages <- NULL
    
    for(year in 26:nyears){
      
      age <- temp[,year,]
      
      age.median <- age %>% 
        as.data.frame()
      
      age.median <- age.median[ , ]*mat[,12]
      age.median <- age.median[ , ]*kg[,12]
        
      age.median <- colSums(age.median[,1:nsim]) %>%   
        as.data.frame() %>% 
        mutate(Year = year) %>% 
        rename(MatBio = ".")
      
      median.ages <- rbind(median.ages, age.median)
    }
    median.ages <- median.ages %>% 
      group_by(Year) %>% 
      summarise(Median_MatBio = median(MatBio), 
                P_0.025 = quantile(MatBio, probs=c(0.025)),
                P_0.975 = quantile(MatBio, probs=c(0.975))) %>% 
      mutate(Scenario = scenario.names[i])

    total_pop <- rbind(total_pop, median.ages)
  
  }
  
  return(total_pop)
  
}

# Separate NTZs and Fished areas 

zone.pop.format <- function(ntz.list, f.list, scenario.name, nsim){

NTZ_Ages_S00 <- NULL
F_Ages_S00 <- NULL

for(SIM in 1:nsim){
  
  temp <- as.data.frame(colSums(ntz.list[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  NTZ_Ages_S00 <- cbind(NTZ_Ages_S00, temp$Number)
  
  temp <- as.data.frame(colSums(f.list[[SIM]])) %>% 
    mutate(Age = seq(1:30)) %>% 
    pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
    mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
  
  F_Ages_S00 <- cbind(F_Ages_S00,temp$Number)
}

NTZ_Ages_S00 <- as.data.frame(NTZ_Ages_S00) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaNT)) %>%
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:nsim]))) %>% 
  mutate(P_0.025 = rowQuantiles(as.matrix(.[,3:nsim]), probs=c(0.025))) %>%
  mutate(P_0.975 = rowQuantiles(as.matrix(.[,3:nsim]), probs=c(0.975))) %>% 
  mutate(Scenario = scenario.name)  %>% 
  dplyr::select(Mod_Year, Stage, Median_Pop, P_0.025, P_0.975, Scenario)


F_Ages_S00 <- as.data.frame(F_Ages_S00) %>% 
  mutate(Age = rep(1:30, each=59)) %>% 
  mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
  mutate(Stage = ifelse(Age==1, "Recruit",
                        ifelse(Age>1 & Age<3, "Sublegal",
                               ifelse(Age>=3 & Age<=10, "Legal",
                                      ifelse(Age>10, "Large Legal",NA))))) %>% 
  group_by(Stage, Mod_Year) %>% 
  summarise(across(where(is.numeric) & !Age, sum)) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric) & !Mod_Year, ~./AreaFished)) %>% 
  mutate(Median_Pop = rowMedians(as.matrix(.[,3:nsim]))) %>% 
  mutate(P_0.025 = rowQuantiles(as.matrix(.[,3:nsim]), probs=c(0.025))) %>%
  mutate(P_0.975 = rowQuantiles(as.matrix(.[,3:nsim]), probs=c(0.975))) %>% 
  mutate(Scenario = scenario.name)  %>% 
  dplyr::select(Mod_Year, Stage, Median_Pop, P_0.025, P_0.975, Scenario)

NTZ.F.Ages <- list()
NTZ.F.Ages[[1]] <- NTZ_Ages_S00
NTZ.F.Ages[[2]] <- F_Ages_S00

return(NTZ.F.Ages)

}


  
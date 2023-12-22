#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec movementfunc_cpp(const int AGE, const int MONTH, const int MaxCell, arma::mat AdultMove, arma::cube YearlyTotal) {
  
  // Check regarding indexing in rccp vs R, e.g. should Month=1 in R by Month=0 in rccp?
  int Cell_cl;
  int Cell_rw;
  double Pop;
    
  arma::vec All_Movers(MaxCell);
  arma::mat Pop2(MaxCell,MaxCell);
  
  for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) { 
    for (Cell_cl=0; Cell_cl<MaxCell; Cell_cl++) { 
      
      Pop = YearlyTotal(Cell_rw,MONTH,AGE);
      
      Pop2(Cell_rw,Cell_cl) = Pop2(Cell_rw,Cell_cl) + (AdultMove(Cell_rw,Cell_cl) * Pop); 
    } 
  } 
  
  for (Cell_cl=0; Cell_cl<MaxCell; Cell_cl++) { 
    for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) { 
      All_Movers(Cell_cl)=All_Movers(Cell_cl) + Pop2(Cell_rw,Cell_cl);
    }
  }
    
  return All_Movers;

}

// [[Rcpp::export]]
Rcpp::List mortalityfunc_cpp(const int AGE, const int MaxCell, const int MONTH, const int YEAR, const double NatMort,
                             arma::mat Weight, arma::cube Selectivity, arma::cube YearlyTotal, Rcpp::NumericVector TotalEffort) {

  int Cell_rw;
  //double NatMort;
  arma::vec tot_survived(MaxCell);
  arma::vec tot_died(MaxCell);
  arma::vec caught(MaxCell);
  arma::vec finite_f(MaxCell);
  arma::vec caught_weight(MaxCell);
  arma::vec Z(MaxCell);
  arma::vec YPR(MaxCell);
  
  // arma::vec tot_survived(MaxCell);
  // arma::vec sel_effort(MaxCell);
  // arma::vec temp_catch(MaxCell);
  // arma::vec Fish_catch(MaxCell);
  
  //Nat_Mort = 1-exp(-NatMort/12);
  for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
    tot_survived(Cell_rw) = YearlyTotal(Cell_rw,MONTH,AGE) * ((exp(-TotalEffort(Cell_rw) * Selectivity(AGE,MONTH,YEAR))) * exp(-NatMort/12));
    caught(Cell_rw) = YearlyTotal(Cell_rw,MONTH,AGE) * (1-(exp(-TotalEffort(Cell_rw) * Selectivity(AGE,MONTH,YEAR))));
    // tot_survived(Cell_rw) = YearlyTotal(Cell_rw,MONTH,AGE) * (1-((Effort(Cell_rw,MONTH,YEAR) * Selectivity(AGE,MONTH,YEAR)) + Nat_Mort));
    // caught(Cell_rw) = YearlyTotal(Cell_rw,MONTH,AGE) * (Effort(Cell_rw,MONTH,YEAR) * Selectivity(AGE,MONTH,YEAR));
    tot_died(Cell_rw) = YearlyTotal(Cell_rw,MONTH,AGE) - tot_survived(Cell_rw);
    finite_f(Cell_rw) = (exp(-TotalEffort(Cell_rw) * Selectivity(AGE,MONTH,YEAR)));
    caught_weight(Cell_rw) = caught(Cell_rw) * Weight(AGE, MONTH);
    
    Z(Cell_rw) = (TotalEffort(Cell_rw) * Selectivity(AGE,MONTH,YEAR)) + NatMort/12;
    YPR(Cell_rw) = tot_survived(Cell_rw)*((TotalEffort(Cell_rw) * Selectivity(AGE,MONTH,YEAR))/Z(Cell_rw)) * (1-(exp(-Z(Cell_rw)))) * Weight(AGE, MONTH);
  }


  //} //cell row loop
  
  return Rcpp::List::create(Rcpp::Named("tot_survived") = tot_survived,
                            Rcpp::Named("tot_died") = tot_died,
                            Rcpp::Named("Catch") = caught,
                            Rcpp::Named("finite_f") = finite_f,
                            Rcpp::Named("Catch_weight") = caught_weight,
                            Rcpp::Named("YPR") = YPR);
}


// [[Rcpp::export]]
Rcpp::List recruitmentfunc_cpp(const int MaxCell, const int MaxAge, const double BHa, const double BHb, const double PF, 
                              arma::mat Mature, arma::mat Weight, arma::vec Settlement, arma::cube YearlyTotal) {
  
  
  int AGE;
  int Cell_rw;
  double Fem_adults; 
  double tot_recs;
  arma::vec SB(MaxCell);
  arma::vec recs(MaxAge);
  arma::vec recs_variable(MaxAge);
  arma::vec settle_recs(MaxCell);
  arma::vec TotFemSB(MaxAge);
  
  //TotFemSB = 0;
  for (AGE=0; AGE<MaxAge; AGE++) { 
    for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) { 
      Fem_adults = YearlyTotal(Cell_rw,9,AGE) * PF;  // month=October Female adults of that age class
      SB(Cell_rw) = Fem_adults * (Mature(AGE,9)) * (Weight(AGE,9)); // Gives us biomass
    }
    TotFemSB(AGE) = sum(SB); // adding up across the cells - all biomass for one age group
    
    recs(AGE) = (TotFemSB(AGE) / (BHa + BHb*TotFemSB(AGE))); // calculate the number of recs from that age group
    recs_variable(AGE) = recs(AGE); 
  }
  //std::cout << "vec recruitmentfunc_cpp: TotFemSB " << TotFemSB << std::endl;
  
  tot_recs = (sum(recs_variable)) * (R::rlnorm(1, 0.6));
  //std::cout << "vec recruitmentfunc_cpp: recs " << recs << std::endl;
    
    for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) { 
      settle_recs(Cell_rw) = Settlement(Cell_rw) * tot_recs; 
    }
    
    return Rcpp::List::create(Rcpp::Named("settle_recs") = settle_recs,
                              Rcpp::Named("BH_recs") = recs,
                              Rcpp::Named("Fem_adults") = Fem_adults, // All these go once fixed
                              Rcpp::Named("SB") = SB,
                              Rcpp::Named("TotFemSB") = TotFemSB,
                              Rcpp::Named("tot_recs") = tot_recs);
    
    //return settle_recs;
}


// [[Rcpp::export]]
Rcpp::List effortfunc_cpp(const int MaxAge, const int NBR, const int MaxCell, const int MONTH, const int YEAR,        // MONTH and YEAR will become part of larger loop
                          arma::vec FuelPrice, arma::vec CellArea_log, arma::mat Fished,
                          arma::mat DistBR, arma::mat Spatial_q,
                          arma::cube YearlyTotal, arma::cube Select, arma::cube BoatDays){
  int Cell_rw;
  int Cell_col;
  int AGE;
  int ramp;
  double U;
  
  Rcpp::List Res;
  
  Rcpp::NumericVector RowU(NBR);
  //Rcpp::NumericVector Check(NBR);
  Rcpp::NumericVector Catchable(MaxCell);
  Rcpp::NumericVector Total_Effort(MaxCell);
  
  arma::mat CellU(MaxCell,NBR);
  arma::mat CellCoef(MaxCell,NBR);
  arma::mat BR_U(MaxCell,NBR);
  arma::mat Catchable_Age(MaxCell,MaxAge);
  arma::mat Effort(MaxCell,NBR);
  
  //Calculates the number of fish that can be caught in each cell in each month
  for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
    for(AGE=0; AGE<MaxAge; AGE++){
      
      Catchable_Age(Cell_rw,AGE) = YearlyTotal(Cell_rw,MONTH,AGE) * Select(AGE,MONTH,YEAR);
    }
    
    for(Cell_col=0; Cell_col<MaxAge; Cell_col++){
      Catchable(Cell_rw) = Catchable(Cell_rw)+Catchable_Age(Cell_rw, Cell_col);
    }
    
  } // End catchable loop
  
  //std::cout << "CharlottesModelfunc_cpp: " << "Catchable Done " << std::endl;
  
  // Calculates the coefficients for each of the cells
  for(Cell_col=0; Cell_col<NBR; Cell_col++){
    for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
      CellCoef(Cell_rw, Cell_col) = (FuelPrice(MONTH)*DistBR(Cell_rw, Cell_col)) + Catchable(Cell_rw);
    } // CHECK THIS WITH MATT
  } // End cell coefficients
  
  //std::cout << "CharlottesModelfunc_cpp: " << "Coefficients Done " << std::endl;
  
  // Calculate the utility of each cell
  for(ramp=0; ramp<NBR; ramp++){
    for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
      U = exp(CellCoef(Cell_rw,ramp)+CellArea_log(Cell_rw));
      CellU(Cell_rw, ramp) = U;
    }
  } // End cell utility
  
  //std::cout << "CharlottesModelfunc_cpp: " << "Cell utility Done " << std::endl;
  
  // Change values of no-take cells to 0 so they don't impact the sum utility
  for(ramp=0; ramp<NBR; ramp++){
    for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
      if(YEAR>24){
        CellU(Cell_rw, ramp) = Fished(Cell_rw,0) * CellU(Cell_rw, ramp);
      }
      if(YEAR>44 & YEAR<56){
        CellU(Cell_rw, ramp) = Fished(Cell_rw,1) * CellU(Cell_rw, ramp);
      }
      if(YEAR>56){
        CellU(Cell_rw, ramp) = Fished(Cell_rw,2) * CellU(Cell_rw, ramp);
      }
    }
  } // End setting NTZs
  
  //std::cout << "CharlottesModelfunc_cpp: " << "NTZs Done " << std::endl;
  
  // Sum utlity to get utility of going to all cells from each boat ramp
  for(ramp=0; ramp<NBR; ramp++){
    for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
      RowU(ramp) = RowU(ramp) + CellU(Cell_rw,ramp);
    }
  } // End summing utility
  
  //std::cout << "CharlottesModelfunc_cpp: " << "Summing Utility Done " << std::endl;
  
  // Work out utility of going to that cell from each boat ramp
  for (ramp=0; ramp<NBR; ramp++){
    for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
      BR_U(Cell_rw,ramp) = CellU(Cell_rw,ramp)/(RowU(ramp));
    }
  } // End working out utility from boat ramps
  
  //std::cout << "CharlottesModelfunc_cpp: " << "Calculating utility Done " << std::endl;
  
  // // Check everything sums to 0
  // for (ramp=0; ramp<NBR; ramp++){
  //   for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
  //     Check(ramp) = Check(ramp) + BR_U(Cell_rw,ramp);
  //   }
  // } // End check loop
  
  // Allocate trips to each of the cells
  for(ramp=0; ramp<NBR; ramp++){
    for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
      Effort(Cell_rw,ramp) = BR_U(Cell_rw,ramp)*BoatDays(MONTH,ramp,YEAR);
    }
  } // End allocating trips loop
  //std::cout << "CharlottesModelfunc_cpp: " << "Allocating Done " << std::endl;
  
  // Sum across ramps
  for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
    for(ramp=0; ramp<NBR; ramp++){
      Total_Effort(Cell_rw) = Total_Effort(Cell_rw)+Effort(Cell_rw, ramp);
    }
  } // End summing ramps loop
  
  //std::cout << "CharlottesModelfunc_cpp: " << "Summing ramps Done " << std::endl;
  
  // Add in catchability
  for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
    Total_Effort(Cell_rw) = Total_Effort(Cell_rw)*Spatial_q(Cell_rw, YEAR);
    
    // NaN (or Inf) won't compare to itself so we can use this to set it to 0
    if(Total_Effort[Cell_rw] != Total_Effort[Cell_rw] ){ 
      Total_Effort[Cell_rw] = 0; 
    } 
    
  } // End catchability loop
  //std::cout << "CharlottesModelfunc_cpp: " << "Catchability Done " << std::endl;
  
  return Rcpp::List::create(Rcpp::Named("Catchable_Fish") = Catchable,
                            Rcpp::Named("Catchable_Age") = Catchable_Age,
                            Rcpp::Named("CellCoef") = CellCoef,
                            Rcpp::Named("CellU") = CellU,
                            Rcpp::Named("RowU") = RowU,
                            //Rcpp::Named("Check") = Check,
                            Rcpp::Named("Effort") = Effort,
                            Rcpp::Named("Total_Effort") = Total_Effort
  );
}


// [[Rcpp::export]]
Rcpp::List RunModelfunc_cpp(const int YEAR, const int MaxAge, const int MaxYear, const int MaxCell, const int NBR, const double NatMort, const double BHa, const double BHb, const double PF, 
                            arma::vec Settlement, arma::vec FuelPrice, arma::vec CellArea_log, arma::mat Fished,
                                   arma::mat AdultMove, arma::mat Mature, arma::mat Weight, arma::mat DistBR, arma::mat Spatial_q,
                                   arma::cube YearlyTotal, arma::cube Select, arma::cube BoatDays) {
  
  // int Year;
  int MONTH;
  int AGE;
  int Cell_rw;
  // int Cell_col;
  // int ramp;
  // double U;
  
  Rcpp::List Res;
  Rcpp::List recruits;
  
  Rcpp::NumericVector Fish_catch(MaxCell);
  Rcpp::NumericVector tot_mort(MaxCell);
  Rcpp::NumericVector caught(MaxCell);
  Rcpp::NumericVector caught_weight(MaxCell);
  Rcpp::NumericVector tot_survived(MaxCell);
  Rcpp::NumericVector finite_f(MaxCell);
  Rcpp::NumericVector tot_died(MaxCell);
  Rcpp::NumericVector Weight_catch(MaxCell);
  Rcpp::NumericVector All_Movers(MaxCell);
  Rcpp::NumericVector settle_recs(MaxCell);
  Rcpp::NumericVector BH_recs(MaxCell);
  Rcpp::NumericVector Fem_SB(MaxAge);
  Rcpp::NumericVector Z(MaxCell);
  Rcpp::NumericVector YPR(MaxCell);
  Rcpp::NumericVector YPR_month(MaxCell);
  Rcpp::NumericVector TotalEffort(MaxCell);
  
  arma::cube month_catch(MaxCell, 12, MaxAge);
  arma::cube month_catch_weight(MaxCell, 12, MaxAge);
  arma::cube age_catch(MaxCell, 12, MaxAge);
  arma::cube tot_mort_all(MaxCell, 12, MaxAge);
  arma::cube age_dead(MaxCell, 12, MaxAge);
  arma::cube age_survived(MaxCell, 12, MaxAge);
  arma::cube month_YPR(MaxCell, 12, MaxAge);
  
  // for (Year=0; Year<MaxYear; Year++) { 
    for (MONTH=0; MONTH<12; MONTH++) { 
      // std::cout << " Month " << MONTH << std::endl; 
      
      // move fish
      for (AGE=0; AGE<MaxAge; AGE++) { 
        
        All_Movers = movementfunc_cpp(AGE, MONTH, MaxCell, AdultMove, YearlyTotal);

        for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
          YearlyTotal(Cell_rw,MONTH,AGE) = All_Movers(Cell_rw);
        } // Cell_rw
        
      } // AGE
      
      //std::cout << "CharlottesModelfunc_cpp: " << "Movement Done " << std::endl;
      
      // Calcualte the effort for the next month here 
      Res = effortfunc_cpp(MaxAge, NBR, MaxCell, MONTH, YEAR, FuelPrice, CellArea_log, Fished,
                           DistBR, Spatial_q, YearlyTotal, Select, BoatDays);
      
      TotalEffort = Res["Total_Effort"];
      
      //std::cout << "CharlottesModelfunc_cpp: " << "Effort Done " << TotalEffort <<std::endl;
      
      // kill fish
      for (AGE=0; AGE<MaxAge; AGE++) {
        
        if (MONTH==11) {
          // if (AGE<=1) {
            if (AGE<(MaxAge-1)) {
              
              Res = mortalityfunc_cpp(AGE, MaxCell, MONTH, YEAR, NatMort, Weight, Select, YearlyTotal, TotalEffort);
              
              tot_survived = Res["tot_survived"];
              tot_died = Res["tot_died"];
              caught = Res["Catch"];
              caught_weight = Res["Catch_weight"];
              YPR_month = Res["YPR"];
              
              
              for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
                YearlyTotal(Cell_rw,0,AGE+1) = tot_survived(Cell_rw);
              }
            } // AGE<MaxAge
          // } // AGE>1
        } else if (MONTH<11){ // Month==11
          
          Res = mortalityfunc_cpp(AGE, MaxCell, MONTH, YEAR, NatMort, Weight, Select, YearlyTotal, TotalEffort);
          
          tot_survived = Res["tot_survived"];
          tot_died = Res["tot_died"];
          caught = Res["Catch"];
          caught_weight = Res["Catch_weight"];
          YPR_month = Res["YPR"];
          
          for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
            YearlyTotal(Cell_rw,MONTH+1,AGE) = tot_survived(Cell_rw);
          }
        } // if
        // setting other outputs 
        for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
          month_catch(Cell_rw,MONTH,AGE) = caught(Cell_rw);
          month_catch_weight(Cell_rw,MONTH,AGE) = caught_weight(Cell_rw);
          age_dead(Cell_rw,MONTH,AGE) = tot_died(Cell_rw);
          age_survived(Cell_rw,MONTH,AGE) = tot_survived(Cell_rw);
          month_YPR(Cell_rw, MONTH, AGE) = YPR_month(Cell_rw);
        }
      } // AGE
      
      //std::cout << "CharlottesModelfunc_cpp: " << "Mortality Done " << std::endl;
      
      // recruit fish
      if (MONTH==9) { // October
        recruits = recruitmentfunc_cpp(MaxCell, MaxAge, BHa, BHb, PF, Mature, Weight, Settlement, YearlyTotal);
        
        settle_recs = recruits["settle_recs"];
        BH_recs = recruits["BH_recs"];
        Fem_SB = recruits["TotFemSB"];
        // std::cout << "CharlottesModelfunc_cpp: settle_recs " << settle_recs << std::endl;
        
        for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
          YearlyTotal(Cell_rw,0,0) = settle_recs(Cell_rw);
        }
      }
     // std::cout << "CharlottesModelfunc_cpp: " << "Recruitment Done " << std::endl;
    } // Month
 // } // Year
 
  return Rcpp::List::create(Rcpp::Named("YearlyTotal") = YearlyTotal,
                            Rcpp::Named("settle_recs") = settle_recs,
                            Rcpp::Named("BH_recs") = BH_recs,
                            Rcpp::Named("month_catch") = month_catch,
                            Rcpp::Named("month_catch_weight") = month_catch_weight,
                            Rcpp::Named("age_died") = age_dead,
                            Rcpp::Named("age_survived") = age_survived,
                            Rcpp::Named("Monthly_YPR_age") = month_YPR,
                            Rcpp::Named("Fem_SB") = Fem_SB);
  
} // End function






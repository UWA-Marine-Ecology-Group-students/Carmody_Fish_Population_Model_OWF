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
                             arma::mat Weight, arma::cube Selectivity, arma::cube YearlyTotal, arma::cube Effort) {

  int Cell_rw;
  arma::vec tot_survived(MaxCell);
  arma::vec tot_died(MaxCell);
  arma::vec caught(MaxCell);
  arma::vec finite_f(MaxCell);
  
  // arma::vec tot_survived(MaxCell);
  // arma::vec sel_effort(MaxCell);
  // arma::vec temp_catch(MaxCell);
  // arma::vec Fish_catch(MaxCell);
  
 
  for(Cell_rw=0; Cell_rw<MaxCell; Cell_rw++){
    tot_survived(Cell_rw) = YearlyTotal(Cell_rw,MONTH,AGE) * ((exp(-Effort(Cell_rw,MONTH,YEAR) * Selectivity(AGE,MONTH,YEAR))) * exp(-NatMort/12));
    caught(Cell_rw) = YearlyTotal(Cell_rw,MONTH,AGE) * (1-(exp(-Effort(Cell_rw,MONTH,YEAR) * Selectivity(AGE,MONTH,YEAR))));
    tot_died(Cell_rw) = YearlyTotal(Cell_rw,MONTH,AGE) - tot_survived(Cell_rw);
    finite_f(Cell_rw) = (exp(-Effort(Cell_rw,MONTH,YEAR) * Selectivity(AGE,MONTH,YEAR)));
  }
    //tot_survived(Cell_rw) = YearlyTotal(Cell_rw,MONTH,AGE) - tot_died(Cell_rw);

  //} //cell row loop
  
  return Rcpp::List::create(Rcpp::Named("tot_survived") = tot_survived,
                            Rcpp::Named("tot_died") = tot_died,
                            Rcpp::Named("Catch") = caught,
                            Rcpp::Named("finite_f") = finite_f);
}


// [[Rcpp::export]]
Rcpp::List recruitmentfunc_cpp(const int MaxCell, const int MaxAge, const double BHa, const double BHb, const double PF, 
                              arma::mat Mature, arma::mat Weight, arma::vec Settlement, arma::cube YearlyTotal) {
  
  
  int AGE;
  int Cell_rw;
  double Fem_adults; 
  double TotFemSB;
  double tot_recs;
  arma::vec SB(MaxCell);
  arma::vec recs(MaxAge);
  arma::vec recs_variable(MaxAge);
  arma::vec settle_recs(MaxCell);
  
  TotFemSB = 0;
  for (AGE=0; AGE<MaxAge; AGE++) { 
    for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) { 
      Fem_adults = YearlyTotal(Cell_rw,9,AGE) * PF;  // month=October Female adults of that age class
      SB(Cell_rw) = Fem_adults * (Mature(AGE,9)) * (Weight(AGE,9)); // Gives us biomass
    }
    TotFemSB = sum(SB); // adding up across the rows - all biomass for one age group
    
    recs(AGE) = (TotFemSB / (BHa + BHb*TotFemSB)); // calculate the number of recs from that age group
    recs_variable(AGE) = recs(AGE); //*R::rlnorm(0.0, 0.6); // same thing but with variability
  }
  //std::cout << "vec recruitmentfunc_cpp: TotFemSB " << TotFemSB << std::endl;
  
  tot_recs = (sum(recs_variable) * R::rlnorm(0.0, 0.6));
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
Rcpp::List RunModelfunc_cpp(const int YEAR, const int MaxAge, const int MaxYear, const int MaxCell, const double NatMort, const double BHa, const double BHb, const double PF, 
                                   arma::mat AdultMove, arma::mat Mature, arma::mat Weight, arma::vec Settlement,
                                   arma::cube YearlyTotal, arma::cube Select, arma::cube Effort) {
  
  // int Year;
  int MONTH;
  int AGE;
  int Cell_rw;
  Rcpp::List Res;
  Rcpp::List recruits;
  Rcpp::NumericVector Fish_catch(MaxCell);
  Rcpp::NumericVector tot_mort(MaxCell);
  Rcpp::NumericVector caught(MaxCell);
  Rcpp::NumericVector tot_survived(MaxCell);
  Rcpp::NumericVector finite_f(MaxCell);
  Rcpp::NumericVector tot_died(MaxCell);
  Rcpp::NumericVector Weight_catch(MaxCell);
  Rcpp::NumericVector All_Movers(MaxCell);
  Rcpp::NumericVector settle_recs(MaxCell);
  Rcpp::NumericVector BH_recs(MaxCell);
  
  arma::cube month_catch(MaxCell, 12, MaxAge);
  arma::cube age_catch(MaxCell, 12, MaxAge);
  arma::cube tot_mort_all(MaxCell, 12, MaxAge);
  arma::cube age_dead(MaxCell, 12, MaxAge);
  arma::cube age_survived(MaxCell, 12, MaxAge);
  
  // for (Year=0; Year<MaxYear; Year++) { 
    for (MONTH=0; MONTH<12; MONTH++) { 
      // std::cout << " Month " << MONTH << std::endl; 
      // move fish
      for (AGE=0; AGE<MaxAge; AGE++) { 
        
        All_Movers = movementfunc_cpp(AGE, MONTH, MaxCell, AdultMove, YearlyTotal);
        // std::cout << "CharlottesModelfunc_cpp: AGE " << AGE << "All_Movers " << All_Movers << std::endl;

        for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
          YearlyTotal(Cell_rw,MONTH,AGE) = All_Movers(Cell_rw);
        } // Cell_rw
        
      } // AGE
      // kill fish
      for (AGE=0; AGE<MaxAge; AGE++) {
        
        if (MONTH==11) {
          // if (AGE<=1) {
            if (AGE<(MaxAge-1)) {
              
              Res = mortalityfunc_cpp(AGE, MaxCell, MONTH, YEAR, NatMort, Weight, Select, YearlyTotal, Effort);
              
              tot_survived = Res["tot_survived"];
              tot_died = Res["tot_died"];
              caught = Res["Catch"];
              
              for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
                YearlyTotal(Cell_rw,0,AGE+1) = tot_survived(Cell_rw);
              }
            } // AGE<MaxAge
          // } // AGE>1
        } else if (MONTH<11){ // Month==11
          
          Res = mortalityfunc_cpp(AGE, MaxCell, MONTH, YEAR, NatMort, Weight, Select, YearlyTotal, Effort);
          
          tot_survived = Res["tot_survived"];
          tot_died = Res["tot_died"];
          caught = Res["Catch"];
          
          for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
            YearlyTotal(Cell_rw,MONTH+1,AGE) = tot_survived(Cell_rw);
          }
        } // if
        // setting other outputs 
        for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
          month_catch(Cell_rw,MONTH,AGE) = caught(Cell_rw);
          age_dead(Cell_rw,MONTH,AGE) = tot_died(Cell_rw);
          age_survived(Cell_rw,MONTH,AGE) = tot_survived(Cell_rw);
        }
      } // AGE
      // recruit fish
      if (MONTH==9) { // October
        recruits = recruitmentfunc_cpp(MaxCell, MaxAge, BHa, BHb, PF, Mature, Weight, Settlement, YearlyTotal);
        
        settle_recs = recruits["settle_recs"];
        BH_recs = recruits["BH_recs"];
        // std::cout << "CharlottesModelfunc_cpp: settle_recs " << settle_recs << std::endl;
        
        for (Cell_rw=0; Cell_rw<MaxCell; Cell_rw++) {
          YearlyTotal(Cell_rw,0,0) = settle_recs(Cell_rw);
        }
      }
    } // Month
 // } // Year
 
  return Rcpp::List::create(Rcpp::Named("YearlyTotal") = YearlyTotal,
                            Rcpp::Named("settle_recs") = settle_recs,
                            Rcpp::Named("BH_recs") = BH_recs,
                            Rcpp::Named("month_catch") = month_catch,
                            Rcpp::Named("age_died") = age_dead,
                            Rcpp::Named("age_survived") = age_survived);
  
} // End function


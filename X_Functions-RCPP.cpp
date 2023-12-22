#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

// // [[Rcpp::export]]
// Rcpp::List test1(arma::vec avec, arma::vec avec2) {
// // simple example, return named list object
//   
//   arma::vec bvec(2); 
//   arma::vec bvec2(2); 
//   
//   bvec = avec;
//   bvec2 = avec2;
//   return Rcpp::List::create(Rcpp::Named("bvec") = bvec,
//                             Rcpp::Named("bvec2") = bvec2);
// }
// 
// // [[Rcpp::export]]
// Rcpp::List test2(arma::vec avec, arma::vec avec2) {
// // simple example, access named list objects from another function, and return named list objects
//   
//   Rcpp::NumericVector bvec(2); 
//   Rcpp::NumericVector bvec2(2); 
//   Rcpp::List Res;
// 
//   Res = test1(avec,avec2);
//   bvec = Res["bvec"];  
//   bvec2 = Res["bvec2"];
//   
//   return Rcpp::List::create(Rcpp::Named("bvec") = bvec,
//                             Rcpp::Named("bvec2") = bvec2);
// }

// [[Rcpp::export]]
arma::vec movementfunc_cpp(const int Age, const int Month, const int Max_Cell, arma::mat Adult_Move, arma::cube Population) {
  
  // Check regarding indexing in rccp vs R, e.g. should Month=1 in R by Month=0 in rccp?
  int Cell_cl;
  int Cell_rw;
  double Pop;
    
  arma::vec All_Movers(Max_Cell);
  arma::mat Pop2(Max_Cell,Max_Cell);
  
  for (Cell_rw=0; Cell_rw<Max_Cell; Cell_rw++) { 
    for (Cell_cl=0; Cell_cl<Max_Cell; Cell_cl++) { 
      
      Pop = Population(Cell_rw,Month,Age);
      
      Pop2(Cell_rw,Cell_cl) = Pop2(Cell_rw,Cell_cl) + (Adult_Move(Cell_rw,Cell_cl) * Pop); 
    } 
  } 
  
  for (Cell_cl=0; Cell_cl<Max_Cell; Cell_cl++) { 
    for (Cell_rw=0; Cell_rw<Max_Cell; Cell_rw++) { 
      All_Movers(Cell_cl)=All_Movers(Cell_cl) + Pop2(Cell_rw,Cell_cl);
    }
  }
    
  return All_Movers;

}

// [[Rcpp::export]]
Rcpp::List mortalityfunc_cpp(const int YEAR, const int Age, const int Max_Cell, const int Month, const double Nat_Mort,
                 arma::cube Select, arma::cube Population, arma::cube Effort) {

  int Cell_rw;
  Rcpp::NumericVector tot_survived(Max_Cell);
  Rcpp::NumericVector sel_effort(Max_Cell);
  Rcpp::NumericVector temp_catch(Max_Cell);
  Rcpp::NumericVector Fish_catch(Max_Cell);
  
  // arma::vec tot_survived(Max_Cell);
  // arma::vec sel_effort(Max_Cell);
  // arma::vec temp_catch(Max_Cell);
  // arma::vec Fish_catch(Max_Cell);
  
  for (Cell_rw=0; Cell_rw<Max_Cell; Cell_rw++) { 
    tot_survived(Cell_rw) = Population(Cell_rw,Month,Age) * exp(-Select(Age,Month,Year) * Effort(Cell_rw,Month,Year))*
                              exp(-(Nat_Mort/12.0)); 
    
    // Calculate catch
    sel_effort(Cell_rw) = Select(Age,Month,Year) * Effort(Cell_rw,Month,Year);
    temp_catch(Cell_rw) = sel_effort(Cell_rw) / ((Nat_Mort/12.0) + sel_effort(Cell_rw));
    Fish_catch(Cell_rw) = tot_survived(Cell_rw) * temp_catch(Cell_rw);
  }


  return Rcpp::List::create(Rcpp::Named("Fish_catch") = Fish_catch,
                            Rcpp::Named("tot_survived") = tot_survived);
}


// [[Rcpp::export]]
arma::vec recruitmentfunc_cpp(const int Max_Cell, const int MaxAge, const double BHa, const double BHb, const double PF, 
                              arma::vec Mature, arma::mat Weight, arma::vec settlement, arma::cube Population) {
  
  
  int Age;
  int Cell_rw;
  double Fem_adults; 
  double TotFemSB;
  double recs;
  arma::vec tot_recs(MaxAge);
  arma::vec SB(Max_Cell);
  arma::vec settle_recs(Max_Cell);
  
  TotFemSB = 0;
  for (Age=0; Age<MaxAge; Age++) { 
    for (Cell_rw=0; Cell_rw<Max_Cell; Cell_rw++) { 
      Fem_adults = Population(Cell_rw,9,Age) * PF;  // month=October
      SB(Cell_rw) = Fem_adults * Mature(Age) * (Weight(Age,9) / 1000); // convert weight in g to kg
      TotFemSB = TotFemSB + SB(Cell_rw);
      
      // if (Age < 3) {
      //   if (Cell_rw < 10) {
      //     //std::cout << "vec recruitmentfunc_cpp: Age " << Age << " Cell_rw " << Cell_rw << " Fem_adults " << Fem_adults << std::endl;
      //     std::cout << "vec recruitmentfunc_cpp: Age " << Age << " Cell_rw " << Cell_rw << " Fem_adults " << Fem_adults << std::endl;
      //   }
      // }
    }
  }
  //std::cout << "vec recruitmentfunc_cpp: TotFemSB " << TotFemSB << std::endl;
  
  recs = TotFemSB / (BHa + BHb * TotFemSB); // This gives us males and females to go into the next generation
  //std::cout << "vec recruitmentfunc_cpp: recs " << recs << std::endl;
    
    for (Cell_rw=0; Cell_rw<Max_Cell; Cell_rw++) { 
      settle_recs(Cell_rw) = settlement(Cell_rw) * recs; // Gives us spawning biomass in each age group at the end of the year, hence the x 12+1 as it starts at 1 not zero
    }
    
    return settle_recs;
}


// [[Rcpp::export]]
Rcpp::List RunModelfunc_cpp(const int YEAR, const int MaxAge, const int MaxYear, const int Max_Cell, const double Nat_Mort, const double BHa, const double BHb, const double PF, 
                                   arma::mat Adult_Move, arma::vec Mature, arma::mat Weight, arma::vec settlement,
                                   arma::cube Population, arma::cube Select, arma::cube Effort) {
  
  int Year;
  int Month;
  int Age;
  int Cell_rw;
  double total_population;
  
  Rcpp::List Res;
  
  Rcpp::NumericVector Fish_catch(Max_Cell);
  Rcpp::NumericVector tot_survived(Max_Cell);
  Rcpp::NumericVector All_Movers(Max_Cell);
  Rcpp::NumericVector settle_recs(Max_Cell);
  
  for (Year=0; Year<MaxYear; Year++) { 
    
    for (Month=0; Month<12; Month++) { 
      // move fish
      for (Age=0; Age<MaxAge; Age++) { 
        
        All_Movers = movementfunc_cpp(Age, Month, Max_Cell, Adult_Move, Population);
        //std::cout << "CharlottesModelfunc_cpp: Age " << Age << "All_Movers " << All_Movers << std::endl;            
        
        for (Cell_rw=0; Cell_rw<Max_Cell; Cell_rw++) {
          Population(Cell_rw,Month,Age) = All_Movers(Cell_rw);
        } // Cell_rw
        
      } // Age
      
      // kill fish
      for (Age=0; Age<MaxAge; Age++) {
        if (Month==11) {
          if (Age<=1) {
            if (Age<MaxAge) {
              Res = mortalityfunc_cpp(Age, Max_Cell, Month, Year, Nat_Mort, Select, Population, Effort);
              Fish_catch = Res["Fish_catch"];
              tot_survived = Res["tot_survived"];
              
              for (Cell_rw=0; Cell_rw<Max_Cell; Cell_rw++) {
                Population(Cell_rw,0,Age+1) = tot_survived(Cell_rw);
              }
            } // Age<MaxAge
          } // Age>1
        } else { // Month=11
          //tot_survived = mortalityfunc_cpp(Age, Max_Cell, Month, Year, Nat_Mort, Select, Population, Effort);
          Res = mortalityfunc_cpp(Age, Max_Cell, Month, Year, Nat_Mort, Select, Population, Effort);
          Fish_catch = Res["Fish_catch"];
          tot_survived = Res["tot_survived"];
          for (Cell_rw=0; Cell_rw<Max_Cell; Cell_rw++) {
            Population(Cell_rw,Month+1,Age) = tot_survived(Cell_rw);
          }
        } // if
      } // Age
      
      // recruit fish
      if (Month==9) { // October
        settle_recs = recruitmentfunc_cpp(Max_Cell, MaxAge, BHa, BHb, PF, Mature, Weight, settlement, Population);
        // std::cout << "CharlottesModelfunc_cpp: settle_recs " << settle_recs << std::endl;
        
        for (Cell_rw=0; Cell_rw<Max_Cell; Cell_rw++) {
          Population(Cell_rw,0,0) = settle_recs(Cell_rw);
        }
      }
      
    } // Month
  } // Year
  return Rcpp::List::create(Rcpp::Named("Fish_catch") = Fish_catch,
                            Rcpp::Named("tot_survived") = tot_survived,
                            Rcpp::Named("settle_recs") = settle_recs,
                            Rcpp::Named("Population") = Population);
  
} // End function









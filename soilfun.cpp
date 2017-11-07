#include <Rcpp.h>
using namespace Rcpp;

/// COPIED FROM EcoHydro2D Github

// changes by Anneka: hb in mm! not cm anymore

//
// Soil data as a function
// Ecohydrology WIMEK project
// Willem Vervoort September 2007
// To simplify soil input
// ########################################
// [[Rcpp::export]]
List  Soil_cpp(std::string stype) {
  // double psi_sh = -10.0;
  
  // default = Medium Light Clay
  double n = 0.418; // porosity
  // more soil variables for evaporation & losses
  double K_s = 3.51; // cm/day
  double b = 13.48; // neurotheta LMC
  //  double nvg = 1.089;
  // double  avg  = 0.0591;
  double s_fc = 0.364/n; // Field capacity
  double psi_s_bar = -1.5E-3; // This is the bubbling pressure
  double hb = psi_s_bar*(-1e4);
  // double spec_y = 0.054; //Johnson 1967 says 0.05, specific yield. 
  double h1bar  = -psi_s_bar;
  
  
  if (stype == "L Med Clay Stony") {
    // Medium Light Clay
    n = 0.318; // porosity
    // more soil variables for evaporation & losses
    K_s = 3.51; // cm/day
    b = 13.48; // neurotheta LMC
    // nvg = 1.089;
    // avg = 0.0591;
    s_fc = 0.264/n; // Field capacity
    psi_s_bar = -1.5E-3; // This is the bubbling pressure
    hb = psi_s_bar*(-1e4);
    // spec_y = 0.054; //Johnson 1967 says 0.05, specific yield. 
    h1bar  = -psi_s_bar;
    
  }
  
  if (stype == "S Clay Loam") {
    // Sandy Clay Loam
    n = 0.367; // porosity
    // more soil variables for evaporation & losses
    K_s = 52.08; // cm/day
    b = 6.4069; // neurotheta sandy clay loam
    // avg = 0.0521;
    // nvg = 1.237;
    s_fc = 0.2677/n; // Field capacity
    psi_s_bar = -1.2E-3;
    hb = psi_s_bar*(-1e4);
    // spec_y = 0.07;  //difference Fc and por Johnson 1967 says 0.07 
    h1bar  = -psi_s_bar;
  }
  
  if (stype == "Loamy Sand") {
    // Loamy sand
    n = 0.37; // porosity
    // more soil variables for evaporation & losses
    K_s = 175.3; // cm/day
    b = 4.5206;
    // avg = 0.0641;
    // nvg = 1.344;
    s_fc = 0.2098/n; // Field capacity
    
    psi_s_bar = -0.66E-3; // This is the bubbling pressure
    // spec_y = 0.17;  // changed to 0.1 to increase rise in gw, not difference por and fc
     h1bar  = -psi_s_bar;
  }
  
  if (stype == "H Clay") {
    // Medium Heavy Clay
    n = 0.4473; // porosity
    // more soil variables for evaporation & losses
    K_s = 2.82; // cm/day
    b = 16.1501; // neurotheta Medium heavy clay
    // avg = 0.0613;
    // nvg = 1.086;
    s_fc = 0.3936/n; // Field capacity
    
    psi_s_bar = -1.4e-3;
    
    hb = psi_s_bar*-1.0e5;
    // spec_y = 0.05;  // difference por and fc
     h1bar  = -psi_s_bar;
  }
  
  if (stype == "M Clay") {
    // Medium Clay
    n =0.4391; // porosity
    // more soil variables for evaporation & losses
    K_s = 6.04; // cm/day
    b =  13.5127;
    // avg = 0.0507;
    //doouble nvg = 1.088;
    s_fc = 0.3818/n; // Field capacity
    
    psi_s_bar = -1.75E-3; // This is the bubbling pressure
    hb = psi_s_bar*-1.0e5;
    // spec_y = 0.05; // difference por and fc
    h1bar  = -psi_s_bar;
    
    
  }
  
  if (stype == "C Sand") {
    // Coarse Sand
    n = 0.368; // porosity
    // more soil variables for evaporation & losses
    K_s = 182.68; // cm/day
    b = 4.1152;
    // avg = 0.0712;
    // nvg = 1.392;
    s_fc = 0.1895/n; // Field capacity
    
    psi_s_bar = -0.61E-3; // This is the bubbling pressure
    hb = psi_s_bar*-1.0e5; //mm
    // spec_y = 0.27;  // difference por and fc
    h1bar  = -psi_s_bar;
  }
  
  
  // Other derived parameters
  // double s_h = pow(psi_sh/psi_s_bar,-1/b); // soil moisture at hygroscopic point
  // double beta = 2*b+4; //page 714 Laio et al., 2001a
  // 
  // // Define parameters for Eagleson function
  // // Beta parameters
  // double  beta1 = 2+3/b;
  // // alpha parameter
  // double a1 = 1+(3/2)/(beta1-1);
  // 
  
  // double h1bar  = -psi_s_bar;
  
  return(Rcpp::List::create(Rcpp::Named("n") = n,
                            Rcpp::Named("K_s") = K_s,
                            Rcpp::Named("b") = b,
                            Rcpp::Named("hb") = hb,
                            Rcpp::Named("psi_s_bar") = psi_s_bar,
                            Rcpp::Named("s_fc") = s_fc,
                            Rcpp::Named("h1bar") = h1bar));
}


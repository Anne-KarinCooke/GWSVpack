 #include <cstdlib>
#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;


// just a simple function to make sure the TauDEM output from the GeoTiff.R file has the right format to be used in Rcpp
// [[Rcpp::export]]
mat sub1( mat x, uword e) {
  x.shed_col(e-1);
  x.shed_row(e-1);
  return x;
}

// Infiltration
// [[Rcpp::export]]
double Infil(double h, double P, double alpha_i, double k, double W0){
  
  double I=alpha_i*h*((P+k*W0)/(P+k));
  return I;
}
// Overland flow/Runoff
// [[Rcpp::export]]
double OF(double h, double cn, double Mn, double slope){
  
  double qq = (cn/Mn) * (pow(h,1.666667)) * sqrt(slope);
  return qq;
}
// plant water uptake
// [[Rcpp::export]]
double WU(double M, double P, double gmax, double k1 ) {
  
  double Wu = gmax*(M/(M+k1))*P;
  return Wu;
}
// Plant growth
// [[Rcpp::export]]
double Gr(double M, double P, double c, double gmax, double k1){
  
  double Gro = c*WU(M,P,gmax,k1);
  return Gro;
}

//plant mortality function
// [[Rcpp::export]]
double Mo(double P, double d, double sigmaP, double CM){
  
  double Mort = P + (((d*P)-P)/(1+pow(CM,sigmaP)));
  
  return Mort;
  
}
// Capillary rise and drainage
// [[Rcpp::export]]
double L_n(double M, double Z, double n, double Zr, double b, double hb, double K_s, double psi_s_bar){
  
  double s=M/(n*Zr); 
  double psi = hb * pow(s,-b);
  // double s_fc = pow((Z/hb),(-1/b));
  double m = 2 + 3/b;
  double qf = (pow((Z/hb),(-m))-(pow((psi/hb),-m)))/((1+pow((psi/hb),-m))+(m-1)*pow((Z/hb),-m));
  double flux = K_s * qf;
  
  return flux;
  
}
// Runon
// [[Rcpp::export]]
mat Surface(int ro, int co, int border, mat flowdir, mat flowdirTable, mat qq, mat filler){
  
  
  const double pi = 3.141593; 
 
  mat destination(ro, co,fill::zeros); 
  
  
  
  int a;
  int x;
  int y;
  int ii;
  int jj;
  
  for (ii=1; ii< (ro-1); ii++) {
    // 
    for (jj=1; jj< (co-1); jj++ ){    
      
      for (a=0; a < 8; a++) {
        
        if (flowdir(ii,jj) == flowdirTable(0,a))
          
        {  
          
          x = flowdirTable(1,a);
          y = flowdirTable(2,a);
          
          
          destination(ii+x,jj+y) = qq(ii,jj);
        }
        
        if ((flowdir(ii,jj) >= flowdirTable(0,a)) && (flowdir(ii,jj) <= flowdirTable(0,a+1))) 
        { 
          
          x = flowdirTable(1,a);
          y = flowdirTable(2,a);
          
          
          destination(ii+x,jj+y) = qq(ii,jj) * filler(ii,jj)* (1.0 - ((flowdir(ii,jj) - flowdirTable(0,a+1))/(pi/4)));
          
        }
        
      }
      
    }
  }
  for (int kk=0; kk < co; kk++){
    destination((ro-border),kk) = destination((ro-border),kk) + destination(border,kk);
    destination(border,kk) = destination((border),kk) + destination((ro-border),kk);
  }
  
  for (int ll=0; ll < ro; ll++){
    destination(ll,(border)) = destination(ll,(border)) + destination(ll,(co-border));
    destination(ll,(co-border)) = destination(ll,(co-border)) + destination(ll,(border));
  }
  return destination;
  
  
  
}

// simply the distances from all other cell to one cell
// [[Rcpp::export]]
mat Distances(int ro, int co, double kk, double ll, double dx){  ///kk and ll are the i and j of the certain cell around with the distance matrix is calculated
  
  mat dist(ro, co,fill::zeros); 
  
  int ii;
  int jj;
  
  
  for (ii=0; ii< (ro); ii++) {
    for (jj=0; jj< (co); jj++ ){  
      
      dist(ii,jj) = sqrt(pow((ii-kk),2) + pow((jj-ll),2))*dx; //*dx for output in m
      
    }
  }
  return dist;
} 
// plant Interference function (competition/facilitation)
// [[Rcpp::export]]
double interference(int border, int ro, int co, double kk, double ll, mat Psub, double dx, double b1, double b2, double q1, double q2, double pi = 3.141593){  // doubel L
  
  mat interf(ro, co,fill::zeros);
  mat w(ro, co,fill::zeros);
  
  int ii;
  int jj;
  
  
  for (ii=border; ii< (ro-border); ii++) {
    for (jj=border; jj< (co-border); jj++){
      
      
      
      mat dist = Distances(ro, co, kk, ll, dx);
      
      
      w(ii,jj) = b1 *exp(-pow((dist(ii,jj)/q1),2)) - b2 *exp(-pow((dist(ii,jj)/q2),2));
      
    }
  }
  double sum =0.0;
  for (ii=border; ii< (ro-border); ii++) {
    for (jj=border; jj< (co-border); jj++){
      
      interf(ii,jj) = w(ii,jj) * Psub(ii,jj);
      
      sum = accu(interf);
    }
  }
  return sum;
}

// function for seed diffusion directions (all directions)
// [[Rcpp::export]] 
mat write_DiffdirectionTable() {  // this works
  
  arma::mat DiffdirTable(2,8, fill::zeros);
  
  DiffdirTable(0,0) = 0.0;
  DiffdirTable(0,1) = -1.0;
  DiffdirTable(0,2) = -1.0;
  DiffdirTable(0,3) = -1.0;
  DiffdirTable(0,4) = 0.0;
  DiffdirTable(0,5) = 1.0;
  DiffdirTable(0,6) = 1.0;
  DiffdirTable(0,7) = 1.0;
  
  
  DiffdirTable(1,0) = 1.0;
  DiffdirTable(1,1) = 1.0;
  DiffdirTable(1,2) = 0.0;
  DiffdirTable(1,3) = -1.0;
  DiffdirTable(1,4) = -1.0;
  DiffdirTable(1,5) = -1.0;
  DiffdirTable(1,6) = 0.0;
  DiffdirTable(1,7) = 1.0;
  
  return DiffdirTable;
  
}

// // seed Diffusion
// [[Rcpp::export]] 
mat seedDiffusionLoss(int ro, int co, mat Medium , double Dp, double timeincr){
  
  mat diffloss(ro, co,fill::zeros);
  
  
  int ii;
  int jj;
  
  for (ii=1; ii< (ro-1); ii++) {
    
    for (jj=1; jj< (co-1); jj++ ){
      
      
      diffloss(ii,jj) = Dp * Medium(ii,jj) * timeincr; 
      
      
    }
  }
  
  
  
  return diffloss;
  
}
// [[Rcpp::export]] 
mat seedDiffusionGain(int ro, int co, int border, mat DiffdirTable, mat Medium , double Dp, double timeincr){
  
  mat diffgain(ro, co,fill::zeros);
  float divide = 1.0/8.0;
  // border = border + 1; 
  
  int a;
  int x;
  int y;
  int ii;
  int jj;
  
  for (ii=1; ii< (ro-1); ii++) {
    for (jj=1; jj< (co-1); jj++ ){
      for (a=0; a < 7; a++) {
        
        x = DiffdirTable(0,a);
        y = DiffdirTable(1,a);
        
        diffgain(ii+x,jj+y) = Dp  * Medium(ii,jj) *timeincr * divide;
      }
    }
  }
  for (jj=0; jj < co; jj++){
    diffgain((ro-border),jj) = diffgain((ro-border),jj) + diffgain(0,jj);
    diffgain(0,jj) = diffgain(0,jj) + diffgain((ro-border),jj);
  }
  
  for (ii=0; ii < ro; ii++){
    diffgain(ii,0) = diffgain(ii,0) + diffgain(ii,(co-border));
    diffgain(ii,(co-border)) = diffgain(ii,(co-border)) + diffgain(ii,0);
  }
  
  return diffgain;
}


// [[Rcpp::export]]
List  Veg_cpp() {
  
  double k = 12.0;//Saco et al, 2013
  double W0 = 0.2;//Saco et al, 2013
  double gmax = 0.05;//Saco et al, 2013
  double c = 10.0;//Saco et al, 2013
  double k1 = 5.0;//Saco et al, 2013
  double d = 0.24;//Saco et al, 2013 //fraction of plant mortality
  
  
  
  return(Rcpp::List::create(Rcpp::Named("k") = k,
                            Rcpp::Named("W0") = W0,
                            Rcpp::Named("gmax") = gmax,
                            Rcpp::Named("c") = c,
                            Rcpp::Named("k1") = k1,
                            Rcpp::Named("d") = d));
  
}  
// salinity parameter list function 
// [[Rcpp::export]]
List Salt_cpp(std::string stype) {
  
  // default, NO SALT from nowhere
  double ConcConst = 0.0; //ConcConst is the concentration of the salt in the infiltrating water in g/l
  double CMgw = 0.0; //CMgw is the goundwater salt concentration  in g/l
  double f = 1; //f is the soil salt leaching efficiency (whether some salt is retained)
  
  
  if (stype == "Groundwater") {
    
    ConcConst = 0.0;
    CMgw = 0.01;
    f = 1;
    
  }
  
  if (stype == "Rain") {
    
    ConcConst = 0.01;
    CMgw = 0.0;
    f = 1;
    
  }
  
  if (stype == "Both") {
    
    ConcConst = 0.01;
    CMgw = 0.01;
    f = 1;
    
  }
  
  if (stype == "None") {
    
    ConcConst = 0.0;
    CMgw = 0.0;
    f = 1;
    
  }
  
  
  return(Rcpp::List::create(Rcpp::Named("ConcConst") = ConcConst,
                            Rcpp::Named("CMgw") = CMgw,
                            Rcpp::Named("f") = f));
  
}

// takes elevation, calls and applies TauDEM, outputs flowdir and slope for every cell
// [[Rcpp::export]]
int call_Taudem(arma::mat B){
  
  B.save("B.txt",arma::raw_ascii); 
  system("R CMD BATCH GeoTiff.R");
  return 0;
}
// [[Rcpp::export]]
mat flowdir_load(){
  arma::mat flowdir_new;
  flowdir_new.load("flowdir_new.txt");
  mat flo = sub1(flowdir_new,1);
  
  return flo;
}
// [[Rcpp::export]]
mat slope_load(){
  
  arma::mat slp_matrix_new;
  slp_matrix_new.load("slp_matrix.txt");
  mat slp = sub1(slp_matrix_new,1);
  
  return slp;
}

// translates flowdir from TauDEM (0 to 2pi) into the indeces of the corresponding cells
// [[Rcpp::export]] 
mat write_flowdirTable() {
  
  arma::mat flowdirTable(3,9, fill::zeros);
  
  flowdirTable(0,0) = 0.0;
  flowdirTable(0,1) = 0.7853982;
  flowdirTable(0,2) = 1.570796;
  flowdirTable(0,3) = 2.356195;
  flowdirTable(0,4) = 3.141593; 
  flowdirTable(0,5) = 3.926991;
  flowdirTable(0,6) = 4.712389;
  flowdirTable(0,7) = 5.497788;
  flowdirTable(0,8) = 6.283186;
  
  flowdirTable(1,0) = 0.0;
  flowdirTable(1,1) = -1.0;
  flowdirTable(1,2) = -1.0;
  flowdirTable(1,3) = -1.0;
  flowdirTable(1,4) = 0.0;
  flowdirTable(1,5) = 1.0;
  flowdirTable(1,6) = 1.0;
  flowdirTable(1,7) = 1.0;
  flowdirTable(1,8) = 0.0;
  
  
  flowdirTable(2,0) = 1.0;
  flowdirTable(2,1) = 1.0;
  flowdirTable(2,2) = 0.0;
  flowdirTable(2,3) = -1.0;
  flowdirTable(2,4) = -1.0;
  flowdirTable(2,5) = -1.0;
  flowdirTable(2,6) = 0.0;
  flowdirTable(2,7) = 1.0;
  flowdirTable(2,8) = 1.0;
  
  return flowdirTable;
  
}

// Just a function to find a condition and produce a submatrix (used in following function)
// [[Rcpp::export]]
arma::vec matrix_sub(arma::mat y, double B) {
                         
  arma::vec v = y.elem( find( y > B )); 
  return v;
}
// calculates difference in height from one cell to the next highest cell of its surrounding cells
// [[Rcpp::export]]
mat diff_next_highest(int rows, int cols, mat elev){
  
  arma::mat diff_next_highest_cell(rows, cols, fill::zeros);
  
  int i;
  int j;
  for (i=1; i< (rows-1); i++) {
    
    for (j=1; j< (cols-1); j++ ){
      
      
      mat y = elev.submat(i-1, j-1, i+1, j+1); //D8
      
      
      vec rest = matrix_sub(y,elev(i,j));
      
      NumericVector rest1 = as<NumericVector>(wrap(rest));
      
      double minimum = min(rest1);
      
      diff_next_highest_cell(i,j) =  minimum- elev(i,j);
      
      
    }
  }
  
  return diff_next_highest_cell;
}
// adding slope to elevation grid
// [[Rcpp::export]]
mat add_Slope_to_Elev(mat elev, int rows, double gslp, double ext){  
  
  int iii;
  elev.row(0) = elev.row(0) + (gslp * ext);
  
  for (iii=1; iii< rows; iii++) {
    
    elev.row(iii) = elev.row(iii) + (gslp * (ext-((ext/rows) * (iii-1))));
  }
  return elev;
}
// [[Rcpp::export]]
mat Z_matrix(mat elev, int ro, int co, double Z){
  
  mat m(ro, co, fill::zeros);
  int ii;
  int jj;
  
  for (ii = 0; ii< ro; ii++){
    
    for (jj = 0; jj< co; jj++){
      
      m(ii,jj) = (elev(ii,jj)*1000.0) + Z;
      
    }
  }
  return m;
  
}

// the core of the model
// [[Rcpp::export]]
List TheFunctionMulti(int border, Rcpp::List soilpar,
                  Rcpp::List dims, NumericVector Rain,
                  Rcpp::List fixedInput, Rcpp::List simInput
) {
  
  //general
  int i = 0;
  int j = 0;
  int t = 1;
  int tt = 0;
  
  int deltat = fixedInput["deltat"];
  // temporal discretization
  float timeincr = 1.0/deltat;
  
  // spatial dimensions of the domain
  int rows = dims["rows"];
  int cols = dims["cols"];
  // duration of simulation
  int time = dims["time"];
  // hillslope
  double gslp = simInput["gslp"]; 
  // extension of the domain
  double ext = dims["ext"];
  double dx = ext/rows; 
  
  
  //salt
  double ConcConst_in = simInput["ConcConst"];
  double f_in = fixedInput["f"];
  double CMgw_in = simInput["CMgw"];
  
  //soil
  double n_in = soilpar["n"];
  double b_in = soilpar["b"];
  double K_s_in = soilpar["K_s"];
  K_s_in = K_s_in*10.0;  // corrected to be in mm/d
  double hb_in = soilpar["hb"];
  double psi_s_bar_in = soilpar["psi_s_bar"];
  double h1bar_in = soilpar["h1bar"];
  double s_fc_in = soilpar["s_fc"];
  // infiltration rate
  double alpha_i = fixedInput["alpha_i"]; 
  
  //vegetation
  double k_in = fixedInput["k"];
  double gmax_in = fixedInput["gmax"];
  double c_in = fixedInput["c"];
  double d_in = fixedInput["d"];
  double W0_in = fixedInput["W0"];
//  double k1_in = simInput["k1"];
  double Zr_in = fixedInput["Zr"];
  
  // half saturation constant
  double k1_inA = simInput["k1A"];
  double k1_inB = simInput["k1B"];
  double k1_inC = simInput["k1C"];
  
  // mortality factor, half life?
  // double d_inA = simInput["dA"];
  // double d_inB = simInput["dB"];
  // double d_inC = simInput["dC"];
  
  
  
  
  // plant sensitivity to salinity
 
  // plant sensitivity to salinity species A
  double sigmaPA = simInput["sigmaPA"];
  // plant sensitivity to salinity species B
  double sigmaPB = simInput["sigmaPB"];
  // plant sensitivity to salinity  species C
  double sigmaPC = simInput["sigmaPC"];
  
  
  //groundwater
  // groundwater depth
  double Z = simInput["Z"];
  
  // further parameters
  // Kinematic wave (runoff) paramters
  double cn = fixedInput["cn"];  // conversion factor
  double Mn = fixedInput["Mn"]; // Mannings n
  
  // ecosystem carrying capacity of plant biomass density
  double P0 =  fixedInput["P0"];

  // seed dispersal and diffusion parameters
  double c1 = fixedInput["c1"]; // [1/mm] 
  double c02 = fixedInput["c02"];  //[m/d] 
  //seed diffusivity
  double Dp = fixedInput["Dp"]; 
  
  // interference parameters, competition and facilitation
  
  double b1 = fixedInput["b1"];
  double b2 = fixedInput["b2"];

  // SPECIES A interference parameters, competition and facilitation
  // double b1A = simInput["b1A"];
  // double b2A = simInput["b2A"];
  double q1A = simInput["q1A"];
  double q2A = q1A/0.6;
  
  // SPECIES B interference parameters, competition and facilitation
  // double b1B = simInput["b1B"];
  // double b2B = simInput["b2B"];
  double q1B = simInput["q1B"];
  double q2B = q1B/0.6;
  
  // SPECIES C interference parameters, competition and facilitation
  // double b1C = simInput["b1C"];
  // double b2C = simInput["b2C"];
  double q1C = simInput["q1C"];
  double q2C = q1C/0.6;
  double zeta = fixedInput["zeta"]; // relative importance of non-local effects vs local effects (impact of interference)
  
  
  // //// Creating the DEM
  arma::mat elev1;
  elev1.load("B1.txt"); // DEM input
  mat elev2 = sub1(elev1,1);
  mat elev = add_Slope_to_Elev(elev2, rows, gslp, ext); 
  
  // matrix with elevation differences
  mat diff = diff_next_highest(rows, cols, elev);
  // flow directions
  int bla = call_Taudem(elev);
  mat flowdir_new = flowdir_load();
  //slope of each cell
  mat slope = slope_load();
  
  // Creating the groundwater depth matrix
  mat Zras = Z_matrix(elev, rows, cols,  Z);
  
  
  
  /// storage cubes (arrays for the daily timesteps t:
  
  arma::cube h = arma::zeros(rows, cols, time);
  arma::cube q = arma::zeros(rows, cols, time);
  arma::cube runon = arma::zeros(rows, cols, time);
  arma::cube In = arma::zeros(rows, cols, time);
//  arma::cube P = arma::zeros(rows, cols, time);
//  arma::cube Wu = arma::zeros(rows, cols, time);
  arma::cube M = arma::zeros(rows, cols, time);
  arma::cube flux = arma::zeros(rows, cols, time);
  arma::cube CM = arma::zeros(rows, cols, time);
  arma::cube SmI = arma::zeros(rows, cols, time);
  arma::cube SmM = arma::zeros(rows, cols, time);
  arma::cube Svir = arma::zeros(rows, cols, time);
  arma::cube mb = arma::zeros(rows, cols, time);
  arma::cube mb_salt = arma::zeros(rows, cols, time);
  arma::cube qsd = arma::zeros(rows, cols, time);
  arma::cube runonsd = arma::zeros(rows, cols, time);
  arma::cube seep = arma::zeros(rows, cols, time);
  arma::cube Smh = arma::zeros(rows, cols, time);
  arma::cube Ch = arma::zeros(rows, cols, time);
  arma::cube runonSubsSalt = arma::zeros(rows, cols, time);
  arma::cube salt_runon = arma::zeros(rows, cols, time);
  arma::cube Subsrunon = arma::zeros(rows, cols, time);
  
  /// storage cubes (arrays) for the subdaily time steps tt:
  
  /// HYDROLOGY
  // overland flow depth
  arma::cube h_sub = arma::zeros(rows, cols, deltat);
  // runoff
  arma::cube q_sub = arma::zeros(rows, cols, deltat);
  //runon
  arma::cube runon_sub = arma::zeros(rows, cols, deltat);
  // infiltration
  arma::cube I_sub = arma::zeros(rows, cols, deltat);
  //water mass balance
  arma::cube mb_sub = arma::zeros(rows, cols, deltat);
  //salt mass balance
  arma::cube mb_salt_sub = arma::zeros(rows, cols, deltat);
  //soil moisture
  arma::cube M_sub = arma::zeros(rows, cols, deltat);
  // vertical water flux, capillary rise and drainage
  arma::cube flux_sub = arma::zeros(rows, cols, deltat);
  // Seepage
  arma::cube seep_sub = arma::zeros(rows, cols, deltat);
  
  
  /// VEGETATION
  // /// plant biomass density
  // arma::cube P_sub = arma::zeros(rows, cols, deltat);
  // // plant water uptake
  // arma::cube WU_sub = arma::zeros(rows, cols, deltat);
  // // plant growth
  // arma::cube Gr_sub = arma::zeros(rows, cols, deltat);
  // // plant mortality
  // arma::cube Mo_sub = arma::zeros(rows, cols, deltat);
  // // advective seed transport with runoff 
  // arma::cube qsd_sub = arma::zeros(rows, cols, deltat);
  // // advective seed transport with runon
  // arma::cube runonsd_sub = arma::zeros(rows, cols, deltat);
  
  // of all species together
  
  arma::cube P_sub = arma::zeros(rows, cols, deltat);
  arma::cube P = arma::zeros(rows, cols, time);
  
  arma::cube WU_sub = arma::zeros(rows, cols, deltat);
  arma::cube Wu = arma::zeros(rows, cols, time);
  
  
  // different species
  arma::cube P_subA = arma::zeros(rows, cols, deltat);
  arma::cube P_subB = arma::zeros(rows, cols, deltat);
  arma::cube P_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube Gr_subA = arma::zeros(rows, cols, deltat);
  arma::cube Gr_subB = arma::zeros(rows, cols, deltat);
  arma::cube Gr_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube Mo_subA = arma::zeros(rows, cols, deltat);
  arma::cube Mo_subB = arma::zeros(rows, cols, deltat);
  arma::cube Mo_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube WU_subA = arma::zeros(rows, cols, deltat);
  arma::cube WU_subB = arma::zeros(rows, cols, deltat);
  arma::cube WU_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube qsd_subA = arma::zeros(rows, cols, deltat);
  arma::cube qsd_subB = arma::zeros(rows, cols, deltat);
  arma::cube qsd_subC = arma::zeros(rows, cols, deltat);
  

  // advective seed transport with runon
  arma::cube runonsd_subA = arma::zeros(rows, cols, deltat);
  arma::cube runonsd_subB = arma::zeros(rows, cols, deltat);
  arma::cube runonsd_subC = arma::zeros(rows, cols, deltat);
  
  arma::cube qsdA = arma::zeros(rows, cols, time);
  arma::cube runonsdA = arma::zeros(rows, cols, time);
  arma::cube qsdB = arma::zeros(rows, cols, time);
  arma::cube runonsdB = arma::zeros(rows, cols, time);
  arma::cube qsdC = arma::zeros(rows, cols, time);
  arma::cube runonsdC = arma::zeros(rows, cols, time);
  
  arma::cube P_A = arma::zeros(rows, cols, time);
  arma::cube P_B = arma::zeros(rows, cols, time);
  arma::cube P_C = arma::zeros(rows, cols, time);
  
  arma::cube WU_A = arma::zeros(rows, cols, time);
  arma::cube WU_B = arma::zeros(rows, cols, time);
  arma::cube WU_C = arma::zeros(rows, cols, time);
  
  
  /// SALT
  // soil solution salt concentration
  arma::cube CM_sub = arma::zeros(rows, cols, deltat);
  // salt mass in infiltrating water
  arma::cube SmI_sub = arma::zeros(rows, cols, deltat);
  //soil solution salt mass
  arma::cube SmM_sub = arma::zeros(rows, cols, deltat);
  // virtual saturation
  arma::cube Svir_sub = arma::zeros(rows, cols, deltat);
  // salt transport with surface runon
  arma::cube salt_runon_sub = arma::zeros(rows, cols, deltat);
  // salt rising up with capillary rise
  arma::cube U_salt = arma::zeros(rows, cols, deltat);
  // salt leached due to drainage
  arma::cube L_salt = arma::zeros(rows, cols, deltat);
  // salt mass in h
  arma::cube Smh_sub = arma::zeros(rows, cols, deltat);
  // salt concentration in h
  arma::cube Ch_sub = arma::zeros(rows, cols, deltat);
  
  // actual model loop starts here
  
  for (t = 1; t< time; t++){
    
    
    for (i= border; i< (rows-border); i++) {
      
      for (j= border; j< (cols-border); j++ ){
        
        //initialise cubes at t= 0
        h(i,j,0) = 0.1;

        P(i,j,0) = 150.0;
        P_A(i,j,0) = 50.0;
        P_B(i,j,0) = 50.0;
        P_C(i,j,0) = 50.0;
        M(i,j,0) = 50.0;
        Svir(i,j,0) = 50.0;
        CM(i,j,0) = 0.0;
        SmI(i,j,0) = 0.0;
        SmM(i,j,0) = 0.0;
        Smh(i,j,0) = 0.0;
        Ch(i,j,0) = 0.0;
        seep(i,j,0) = 0.0;
        
        for (tt = 0; tt< (deltat-1); tt++){
          
          
          if(tt == 0) {
            h_sub(i,j,tt) = h(i,j,t-1);
            P_sub(i,j,tt) = P(i,j,t-1);
            P_subA(i,j,tt) = P_A(i,j,t-1);
            P_subB(i,j,tt) = P_B(i,j,t-1);
            P_subC(i,j,tt) = P_C(i,j,t-1);
            M_sub(i,j,tt) = M(i,j,t-1);
            CM_sub(i,j,tt) = CM(i,j,t-1);
            SmI_sub(i,j,tt) = SmI(i,j,t-1);
            SmM_sub(i,j,tt) = SmM(i,j,t-1);
            Svir_sub(i,j,tt) = Svir(i,j,t-1);
            Smh_sub(i,j,tt) = Smh(i,j,t-1);
            Ch_sub(i,j,tt) = Ch(i,j,t-1);
            seep_sub(i,j,tt) = seep(i,j,t-1);
          }
          if((tt==0) & (P(i,j,t-1) < 0.0)){
            P(i,j,t-1) = 0.0;
          }
          
          double Rain_in;
          
          // if ((Rain(t) > 0.0) & (tt == 0)){
          //   Rain_in = 10.0 * Rain(t);
          // } else {
          //   Rain_in = 0.0;
          // }
          Rain_in = timeincr * Rain(t);
          
          // adjust infiltration rate
          if(h_sub(i,j,tt) < (timeincr * K_s_in)) {
            alpha_i = 1.0;
          } else {
            alpha_i = 1.0-((h_sub(i,j,tt) - (timeincr * K_s_in))/h_sub(i,j,tt));
          }
          
          
          // Runoff
          q_sub(i,j,tt+1) =  timeincr * OF(h_sub(i,j,tt), cn, Mn, slope(i,j));
          
          
          // Runon
          mat runon_store(rows, cols, fill::ones);
          mat ones = mat(rows, cols, fill::ones);
          runon_store = Surface(rows, cols, border, flowdir_new, write_flowdirTable(), q_sub.slice(tt), ones);
          runon_sub(i,j,tt+1) = runon_store(i,j);
          
          
          // calculate water depth on soil
          h_sub(i,j,tt+1) =  h_sub(i,j,tt) + Rain_in
            - (timeincr * Infil(h_sub(i,j,tt),P_sub(i,j,tt), alpha_i, k_in, W0_in)) - q_sub(i,j,tt+1) + seep_sub(i,j,tt) + runon_sub(i,j,tt+1);
          
          //    Ponding check and rewrite of elevation mat
          if (((h_sub(i,j,tt+1) + (elev(i,j)*1000.0)) > (diff(i,j)*1000.0)) & (diff(i,j) > 0.0)){ // times 1000 for [mm]

            mat elev_substitute(rows,cols);
            // elev_substitute = elev;
            elev_substitute(i,j) =  elev(i,j) + diff(i,j);

            // calling TauDEM again
            int bla = call_Taudem(elev_substitute);
            // flow directions
            mat flowdir_adj = flowdir_load();
            //slope of each cell
            mat slope_adj = slope_load();

            q_sub(i,j,tt+1) = timeincr * OF(h_sub(i,j,tt+1), cn, Mn, slope_adj(i,j));

            runon_store = Surface(rows, cols, border, flowdir_adj, write_flowdirTable(), q_sub.slice(tt), ones);
            runon_sub(i,j,tt+1) = runon_store(i,j);

            h_sub(i,j,tt+1) =  h_sub(i,j,tt+1) + Rain_in
              - (timeincr * Infil(h_sub(i,j,tt+1),P_sub(i,j,tt), alpha_i, k_in, W0_in)) - q_sub(i,j,tt+1) + seep_sub(i,j,tt) + runon_sub(i,j,tt+1);

          }
          
          //Infiltration
          I_sub(i,j,tt+1) = timeincr * Infil(h_sub(i,j,tt+1), P_sub(i,j,tt) , alpha_i, k_in, W0_in); 
          
          
          // plant water uptake
          WU_subA(i,j,tt+1) = timeincr * WU(Svir_sub(i,j,tt), P_subA(i,j,tt) , gmax_in, k1_inA); 
          WU_subB(i,j,tt+1) = timeincr * WU(Svir_sub(i,j,tt), P_subB(i,j,tt) , gmax_in, k1_inB); 
          WU_subC(i,j,tt+1) = timeincr * WU(Svir_sub(i,j,tt), P_subC(i,j,tt) , gmax_in, k1_inC); 
          
          WU_sub(i,j,tt+1) = WU_subA(i,j,tt+1) + WU_subB(i,j,tt+1) + WU_subC(i,j,tt+1);
          
          // adjustment soil moisture
          M_sub(i,j,tt+1) = M_sub(i,j,tt) + I_sub(i,j,tt+1) - WU_sub(i,j,tt+1);
          
          //           
          // plant growth
          Gr_subA(i,j,tt) = timeincr * Gr(Svir_sub(i,j,tt), P_subA(i,j,tt), c_in, gmax_in, k1_inA); //, P0, sigmaP, CM_sub(i,j,tt));
          Gr_subB(i,j,tt) = timeincr * Gr(Svir_sub(i,j,tt), P_subB(i,j,tt), c_in, gmax_in, k1_inB);
          Gr_subC(i,j,tt) = timeincr * Gr(Svir_sub(i,j,tt), P_subC(i,j,tt), c_in, gmax_in, k1_inC);
          //  plant Mortality
          Mo_subA(i,j,tt) =timeincr * Mo(P_subA(i,j,tt), d_in, sigmaPA, CM_sub(i,j,tt));
          Mo_subB(i,j,tt) =timeincr * Mo(P_subB(i,j,tt), d_in, sigmaPB, CM_sub(i,j,tt));
          Mo_subC(i,j,tt) =timeincr * Mo(P_subC(i,j,tt), d_in, sigmaPC, CM_sub(i,j,tt));
          
          //Species A
          if((q_sub(i,j,tt+1) > c1) & (q_sub(i,j,tt+1)<c02))
          {
            qsd_subA(i,j,tt+1) = (1.0/c1) * q_sub(i,j,tt+1)*P_subA(i,j,tt); //(Saco, 2007)

          }
          
          if((q_sub(i,j,tt+1)*c1) > c02){ //(Saco, 2007)
            qsd_subA(i,j,tt) = c02 * P_subA(i,j,tt);
          }
          // Species B
          
          if((q_sub(i,j,tt+1) > c1) & (q_sub(i,j,tt+1)<c02))
          {
            qsd_subB(i,j,tt+1) = (1.0/c1) * q_sub(i,j,tt+1)*P_subB(i,j,tt); //(Saco, 2007)
     
          }
          
          if((q_sub(i,j,tt+1)*c1) > c02){ //(Saco, 2007)
            qsd_subB(i,j,tt) = c02 * P_subB(i,j,tt);
          }
          //Species C
          if((q_sub(i,j,tt+1) > c1) & (q_sub(i,j,tt+1)<c02))
          {
            qsd_subC(i,j,tt+1) = (1.0/c1) * q_sub(i,j,tt+1)*P_subC(i,j,tt); //(Saco, 2007)

          }
          
          if((q_sub(i,j,tt+1)*c1) > c02){ //(Saco, 2007)
            qsd_subA(i,j,tt) = c02 * P_subA(i,j,tt);
          }
          
 
          
          //seed runon Species A
          mat runonsd_storeA(rows, cols, fill::ones);
          runonsd_storeA = Surface(rows, cols, border,flowdir_new, write_flowdirTable(),qsd_subA.slice(tt),ones);
          runonsd_subA(i,j,tt) = runonsd_storeA(i,j);
          
          //seed runon Species B
          mat runonsd_storeB(rows, cols, fill::ones);
          runonsd_storeB = Surface(rows, cols, border,flowdir_new, write_flowdirTable(),qsd_subB.slice(tt),ones);
          runonsd_subB(i,j,tt) = runonsd_storeB(i,j);
          
          //seed runon Species C
          mat runonsd_storeC(rows, cols, fill::ones);
          runonsd_storeC = Surface(rows, cols, border, flowdir_new, write_flowdirTable(),qsd_subC.slice(tt),ones);
          runonsd_subC(i,j,tt) = runonsd_storeC(i,j);
          
          
          //Germination reduction due to salt
          double germ = (1.0-(M_sub(i,j,tt+1)/Svir_sub(i,j,tt)));
          // 
          // seed diffusion (into all directions due to wind/animals)
          mat seed_diff_lossA = seedDiffusionLoss(rows, cols,  P_subA.slice(tt),  Dp,  timeincr);
          mat seed_diff_gainA = seedDiffusionGain(rows, cols, border, write_DiffdirectionTable(), P_subA.slice(tt),  Dp,  timeincr);
          
          mat seed_diff_lossB = seedDiffusionLoss(rows, cols,  P_subB.slice(tt),  Dp,  timeincr);
          mat seed_diff_gainB = seedDiffusionGain(rows, cols, border, write_DiffdirectionTable(), P_subB.slice(tt),  Dp,  timeincr);
          
          mat seed_diff_lossC = seedDiffusionLoss(rows, cols,  P_subC.slice(tt),  Dp,  timeincr);
          mat seed_diff_gainC = seedDiffusionGain(rows, cols, border, write_DiffdirectionTable(), P_subC.slice(tt),  Dp,  timeincr);
          
          //  Plant biomass balance
          
          double P01A = 0.01 + P0*exp(-sigmaPA*CM_sub(i,j,tt));
          double P01B = 0.01 + P0*exp(-sigmaPB*CM_sub(i,j,tt));
          double P01C = 0.01 + P0*exp(-sigmaPC*CM_sub(i,j,tt));
          
        
          
          //A
          if (interference(border, rows, cols, i,j,P_sub.slice(tt), dx, b1, b2, q1A, q2A) > 0.0) {
            P_subA(i,j,tt+1) = P_subA(i,j,tt) + ((Gr_subA(i,j,tt)  +  runonsd_subA(i,j,tt)  
           + seed_diff_gainA(i,j) + zeta*interference(border, rows, cols, i,j,P_sub.slice(tt), dx, b1, b2, q1A, q2A))*(1.0-(P_sub(i,j,tt)/P01A)))- Mo_subA(i,j,tt) 
            - ((qsd_subA(i,j,tt) + seed_diff_lossA(i,j))*(P01A/P0));
          }
          else {
            P_subA(i,j,tt+1) = P_subA(i,j,tt) + (Gr_subA(i,j,tt)  +  runonsd_subA(i,j,tt)  
                                                 + seed_diff_gainA(i,j))*(1.0 -(P_sub(i,j,tt)/P01A))- Mo_subA(i,j,tt) 
            - ((qsd_subA(i,j,tt) + seed_diff_lossA(i,j) - zeta*interference(border, rows, cols, i,j,P_sub.slice(tt), dx, b1, b2, q1A, q2A))*(P01A/P0));
          }
          // B
          if (interference(border, rows, cols, i,j,P_sub.slice(tt), dx, b1, b2, q1B, q2B) > 0.0) {
            P_subB(i,j,tt+1) = P_subB(i,j,tt) + ((Gr_subB(i,j,tt)  +  runonsd_subB(i,j,tt)  
                                                    + seed_diff_gainB(i,j) + zeta*interference(border, rows, cols, i,j,P_sub.slice(tt), dx, b1, b2, q1B, q2B))*(1.0-(P_sub(i,j,tt)/P01B)))- Mo_subB(i,j,tt) 
            - ((qsd_subB(i,j,tt) + seed_diff_lossB(i,j))*(P01B/P0));
          }
          else {
            P_subB(i,j,tt+1) = P_subB(i,j,tt) + (Gr_subB(i,j,tt)  +  runonsd_subB(i,j,tt)  
                                                   + seed_diff_gainB(i,j))*(1.0 -(P_sub(i,j,tt)/P01B))- Mo_subB(i,j,tt) 
            - ((qsd_subB(i,j,tt) + seed_diff_lossB(i,j) - zeta*interference(border, rows, cols, i,j,P_sub.slice(tt), dx, b1, b2, q1B, q2B))*(P01B/P0));
          }
          
          // C
          if (interference(border, rows, cols, i,j,P_sub.slice(tt), dx, b1, b2, q1B, q2C) > 0.0) {
            P_subC(i,j,tt+1) = P_subC(i,j,tt) + ((Gr_subC(i,j,tt)  +  runonsd_subC(i,j,tt)  
                                                    + seed_diff_gainC(i,j) + zeta*interference(border, rows, cols, i,j,P_sub.slice(tt), dx, b1, b2, q1C, q2C))*(1.0-(P_sub(i,j,tt)/P01C)))- Mo_subC(i,j,tt) 
            - ((qsd_subC(i,j,tt) + seed_diff_lossC(i,j))*(P01C/P0));
          }
          else {
            P_subC(i,j,tt+1) = P_subC(i,j,tt) + (Gr_subC(i,j,tt)  +  runonsd_subC(i,j,tt)  
                                                   + seed_diff_gainC(i,j))*(1.0 -(P_sub(i,j,tt)/P01C))- Mo_subC(i,j,tt) 
            - ((qsd_subC(i,j,tt) + seed_diff_lossC(i,j) - zeta*interference(border, rows, cols, i,j,P_sub.slice(tt), dx, b1, b2, q1C, q2C))*(P01C/P0));
          }
          
          if (P_subA(i,j,tt+1) < 0.0){
            P_subA(i,j,tt+1) = 0.0;
          }
          if (P_subB(i,j,tt+1) < 0.0){
            P_subB(i,j,tt+1) = 0.0;
          }
          if (P_subC(i,j,tt+1) < 0.0){
            P_subC(i,j,tt+1) = 0.0;
          }
          
          // overall biomass of all species combined
          P_sub(i,j,tt+1) = P_subA(i,j,tt+1) + P_subB(i,j,tt+1) + P_subC(i,j,tt+1);
          
          
          //vertical water flux (capillary rise/drainage)
          flux_sub(i,j,tt) = L_n(M_sub(i,j,tt+1),Zras(i,j),n_in,Zr_in,b_in,hb_in,K_s_in,psi_s_bar_in);
          
          //           
          //adjustment for soil mositure balance
          M_sub(i,j,tt+1) = M_sub(i,j,tt+1) +  (timeincr * flux_sub(i,j,tt));
          
          
          // Seepage
          if(M_sub(i,j,tt+1) > (s_fc_in * n_in * Zr_in)){
            
            seep_sub(i,j,tt) = M_sub(i,j,tt+1)-(s_fc_in * n_in * Zr_in);
            // soil moisture adjustment due to seepage (if condition is met)
            M_sub(i,j,tt+1) = M_sub(i,j,tt+1) - seep_sub(i,j,tt);
          }
          
          // salt leaching
          if(flux_sub(i,j,tt) < 0.0 ){
            L_salt(i,j,tt) = f_in * CM_sub(i,j,tt) * (timeincr * flux_sub(i,j,tt));
          } else {
            L_salt(i,j,tt) = 0.0;
          }
          
          // salt upflow
          if(flux_sub(i,j,tt) > 0.0 ) {
            U_salt(i,j,tt) =  CMgw_in * (timeincr * flux_sub(i,j,tt));
          } else {
            U_salt(i,j,tt) = 0.0;
          }
          
          // # salt mass coming in with infiltration
          SmI_sub(i,j,tt+1) =  (I_sub(i,j,tt+1) * Ch_sub(i,j,tt));
          
          // #salt mass balance in soil
          SmM_sub(i,j,tt+1) = SmM_sub(i,j,tt) + SmI_sub(i,j,tt+1) + U_salt(i,j,tt) - L_salt(i,j,tt) - (CM_sub(i,j,tt) * seep_sub(i,j,tt));
          
          //  salt concentration in soil
          
          CM_sub(i,j,tt+1) = SmM_sub(i,j,tt+1)/M_sub(i,j,tt+1);
          
          
          // # Virtual saturation (Shah et al., 2012), here in [mm] to be in the same unit as M
          Svir_sub(i,j,tt+1) = n_in * Zr_in * (pow((h1bar_in * 0.1),(1.0/b_in))) *
            pow((h1bar_in * 0.1 * pow((M_sub(i,j,tt+1)/(n_in * Zr_in)),-b_in))+(3.6 * (CM_sub(i,j,tt+1)*(1.0/58.44))),(-1.0/b_in));
          
          //Stress-gradient-hypothesis
          double beta = Svir_sub(i,j,tt+1)/M_sub(i,j,tt+1);
          b2 = b2 * beta;
          b1 = 1.0 - b2;
          //           
          // surface salt transport with runon
          mat salt_runon_store(rows, cols, fill::ones);
          salt_runon_store = Surface(rows, cols, border, flowdir_new, write_flowdirTable(),q_sub.slice(tt), CM_sub.slice(tt+1));
          salt_runon_sub(i,j,tt) = salt_runon_store(i,j);
          
          //     Balance of salt mass in overland flow depth h
          Smh_sub(i,j,tt+1) = Smh_sub(i,j,tt) + (Rain_in * ConcConst_in) + (CM_sub(i,j,tt+1) * seep_sub(i,j,tt)) - (q_sub(i,j,tt+1)*Ch_sub(i,j,tt))
            +  salt_runon_sub(i,j,tt);
          
          // salt concentration in overland flow depth h
          Ch_sub(i,j,tt+1) = (Smh_sub(i,j,tt+1)/h_sub(i,j,tt+1));
          
                  
          // // salt balance
                   
          mb_salt_sub(i,j,tt) = SmI_sub(i,j,tt) + U_salt(i,j,tt)  - (CM_sub(i,j,tt)*seep_sub(i,j,tt)) - L_salt(i,j,tt);
          // # checking the water mass balance
          mb_sub(i,j,tt) = I_sub(i,j,tt) - WU_sub(i,j,tt) + (flux_sub(i,j,tt) * timeincr);
          
        }
        
        // Aggregating the subdaily to daily values
        
        h(i,j,t) = h_sub(i,j,(deltat-2));
        Wu(i,j,t) = WU_sub(i,j,(deltat-2));
        
        WU_A(i,j,t) = WU_subA(i,j,(deltat-2));
        WU_B(i,j,t) = WU_subB(i,j,(deltat-2));
        WU_C(i,j,t) = WU_subC(i,j,(deltat-2));
        
        P(i,j,t) = P_sub(i,j,(deltat-2));
        
        P_A(i,j,t) = P_subA(i,j,(deltat-2));
        P_B(i,j,t) = P_subB(i,j,(deltat-2));
        P_C(i,j,t) = P_subC(i,j,(deltat-2));
        
        
        M(i,j,t) = M_sub(i,j,(deltat-2));
        CM(i,j,t) = CM_sub(i,j,(deltat-2));
        SmI(i,j,t) = SmI_sub(i,j,(deltat-2));
        SmM(i,j,t) = SmM_sub(i,j,(deltat-2));
        Svir(i,j,t) = Svir_sub(i,j,(deltat-2));
        Ch(i,j,t) = Ch_sub(i,j,(deltat-2));
        Smh(i,j,t) = Smh_sub(i,j,(deltat-2));
        seep(i,j,t) = seep_sub(i,j,(deltat-2));
        
        
        double sumI = 0.0;
        double sumq = 0.0;
        double sumrunon = 0.0;
        double sumflux = 0.0;
        double summb = 0.0;
        double summb_salt = 0.0;
        
        
        double sumqsdA = 0.0;
        double sumqsdB = 0.0;
        double sumqsdC = 0.0;
        
        double sumrunonsdA = 0.0;
        double sumrunonsdB = 0.0;
        double sumrunonsdC = 0.0;
        
        
        double sumsalt_runon = 0.0;
        
        
        for(int tt = 0; tt < deltat; tt++)
        {
          sumI += I_sub(i,j,tt);
          sumrunon += runon_sub(i,j,tt);
          sumq += q_sub(i,j,tt);
          sumflux += flux_sub(i,j,tt) * timeincr;
          summb += mb_sub(i,j,tt);
          summb_salt += mb_salt_sub(i,j,tt);

          sumsalt_runon += salt_runon_sub(i,j,tt);
          
          sumqsdA += qsd_subA(i,j,tt);
          sumrunonsdA += runonsd_subA(i,j,tt);
          
          sumqsdB += qsd_subB(i,j,tt);
          sumrunonsdB += runonsd_subB(i,j,tt);
          
          sumqsdC += qsd_subC(i,j,tt);
          sumrunonsdC += runonsd_subC(i,j,tt);
        }
        
        q(i,j,t) = sumq;
        runon(i,j,t) = sumrunon;
        In(i,j,t) = sumI;
        flux(i,j,t) = sumflux;
        mb(i,j,t) = summb;
        mb_salt(i,j,t) = summb_salt;
        qsdA(i,j,t) = sumqsdA;
        qsdB(i,j,t) = sumqsdB;
        qsdC(i,j,t) = sumqsdC;
        runonsdA(i,j,t) = sumrunonsdA;
        runonsdB(i,j,t) = sumrunonsdB;
        runonsdC(i,j,t) = sumrunonsdC;
        salt_runon(i,j,t) = sumsalt_runon;
        
        
      }
    }
  }
  
  // Creating output list
  
  
  List fields;
  arma::field<arma::cube> f1(30);
  f1( 0 ) = h;
  f1( 1 ) = q;
  f1( 2 ) = In;
  f1( 3 ) = runon;
  f1( 4 ) = Wu;
  f1( 5 ) = WU_A;
  f1( 6 ) = WU_B;
  f1( 7 ) = WU_C;
  f1( 8 ) = P;
  f1( 9 ) = P_A;
  f1( 10 ) = P_B;
  f1( 11 ) = P_C;
  f1( 12 ) = flux;
  f1( 13 ) = M;
  f1( 14 ) = SmM;
  f1( 15 ) = CM;
  f1( 16 ) = mb;
  f1( 17 ) = Svir;
  f1( 18 ) = SmI;
  f1( 19 ) = Smh;
  f1( 20 ) = Ch;
  f1( 21 ) = qsdA;
  f1( 22 ) = qsdB;
  f1( 23 ) = qsdC;
  f1( 24 ) = runonsdA;
  f1( 25 ) = runonsdB;
  f1( 26 ) = runonsdC;
  f1( 27 ) = seep;
  f1( 28 ) = salt_runon;
  f1( 29 ) = mb_salt;
  fields["field<cube>"] = f1;
  
  
  List output = List::create(_["fields : field<cube>"] = f1);
  
  return output;
  
}



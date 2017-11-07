### Script for simulation


# necessities:

# packages for Rtools, Rcpp, RcppArma...
# install.packages("devtools")
library("devtools")
find_rtools() # check wheter it finds Rtools, otherwise Rcpp won't run
install.packages("Rcpp")
library("Rcpp")
install.packages("RcppArmadillo")
library("RcppArmadillo")
# # packages for the patch and diversity analysis 
library("SDMTools")
library("vegan") #



# ***************************************************************

# Rainfall function
Precip <- function(time,alpha,lambda,delta) {
  # generate a vector of times between rainfall events > delta
  f_P<-floor(rexp(time,lambda*exp(-delta/alpha))) # vector of times between rainfall occurrences (equation 4 & 8)
  # generate a binary vector from this (impulse function)
  binary.vec <- unlist(lapply(1:time,function(i) c(rep(0,f_P[i]),1)))
  R <- rexp(length(binary.vec),1/alpha)*binary.vec 
  return(R[1:time])
}

# duration of simulation
time <- 365 # 14600(40 years) # days, 

## Generate rainfall data
delta <- 0
alpha <- 0.91 #cm/event
lambda <- 0.16 #d/event
# Rain <- Precip(time,alpha,lambda,delta)
# Rain <- c(Rain*10) # to be in mm

# Read in actual rainfall data
TennantAirport <- read.csv("M:/Master thesis/tennant creek airport/IDCJAC0009_015135_1800_Data.csv", header= T, sep = ",")
Rain2 <-TennantAirport$Rainfall.amount..millimetres.
Rain2 <- na.omit(Rain2)
Rain2 <- c(Rain2)

# source the parameter list functions
sourceCpp("soilfun.cpp", rebuild = T)

# sOURCE THE MODEL
sourceCpp("Multi_Species_Model.cpp", rebuild=T)

# PARAMETERS THAT STAY THE SAME
# size of the grid
rows <- 20## rows of cells
cols <- 20 ## columns of cells
# EXTEND of PLOT in [m]
ext <- 20.0 
# "diverseInput"
deltat <- 12 # temporal discretization, subdaily timesteps
# infiltration
alpha_i <- 1.0
#Kinematic wave (runoff) paramters
cn <-0.001
Mn <- 0.04 # Manning's n
# soil moisture diffusivity
Dm <- 0.27 #Saco and Moreno-de las Heras,2013
# coefficient impact of plant interference
zeta <- 0.1
# ecosystem carrying capacity [g/m2]
P0 <-  500.0
# seed dispersal and diffusion
c1 <- 2.25 # [1/mm] Saco and Moreno-de las Heras 2013
c02 <- 0.0002   #[m/d] tranformed to [mm/deltat] Saco and Moreno-de las Heras 2013
#seed diffusivity
Dp <- 0.3 #Saco and Moreno-de las Heras,2013
b1 <- 0.5 # relative importance of facilitation [0 -1]
b2 <- 0.5 # relative importance of competition [0 - 1]
Zr <- 300 #mm root depth
f <-1 # leaching factor [0 - 1]

k = 12.0 #Saco et al, 2013
W0 = 0.2 #Saco et al, 2013
gmax = 0.05 #Saco et al, 2013
c = 10.0 #Saco et al, 2013 , plant growth factor
d = 0.02 #Saco et al, 2013 //fraction of plant mortality


# ***************************************************************
## SIMULATION

# declaring the input lists
dims = list(rows=rows,cols=cols,time=time, ext = ext)

fixedInput = list(deltat = deltat, alpha_i = alpha_i, cn = cn, Mn = Mn, P0 = P0, d=d,
                  c1 = c1, c02 = c02, Dp = Dp, zeta = zeta, f=f, Zr = Zr, b1 = b1, b2 = b2,
                  k=k, W0=W0, gmax=gmax, c=c, d=d)

# These simulation input parameters have been varied in the simulations

# groundwater depth
Z <- 1000 #mm
# hillslope of model domain
gslp <- 0.01
# groundwater salt concentration
CMgw <- 3.0 #g/L
# salt concentration in rainfall
ConcConst <-  0.01 #g/L

# water uptake half saturation constant for species A, B, C [mm]
k1A <- 2.5 #
k1B <- 2.5 #
k1C <- 2.5 #
# facilitation range for species A, B, C
q1A <- 0.1 #m
q1B <- 0.5 #m
q1C <- 1.0 #m
#salinity sensitivity for species A, B, C
sigmaPA <- 0.1 # L/g
sigmaPB <- 0.1 # L/g
sigmaPC <- 0.1 # L/g

# loading soil properties (Ksat, porosity etc), 
# other possiblities : L Med Clay Stony, S Clay Loam, Loamy Sand, H Clay, M Clay, C Sand
soilparS <- Soil_cpp("C Sand") 


simInput = list(Z =Z,k1A = k1A,k1B = k1B,k1C = k1C,
                q1A = q1A, q1B = q1B,q1C = q1C, sigmaPA= sigmaPA,sigmaPB= sigmaPB,sigmaPC= sigmaPC,
                gslp = gslp,CMgw = CMgw,
                ConcConst = ConcConst)


# calling the actual model
# adjust choice of "border" (outer cells that are not considered), input parameter lists for soil paramters, dimensions
#specific simulation parameters, rainfall input

resultsMulti <- TheFunctionMulti(border = 1, 
                                soilpar=soilparS, 
                                dims = dims,Rain = Rain2, 
                                fixedInput = fixedInput, 
                                simInput=simInput)


# loading soil properties (Ksat, porosity etc), 
# other possiblities : L Med Clay Stony, S Clay Loam, Loamy Sand, H Clay, M Clay, C Sand

soilparS <- Soil_cpp("C Sand") 

simInput = list(Z =Z,k1 = k1,
                q1 = q1,sigmaP= sigmaP,
                gslp = gslp,CMgw = CMgw,
                ConcConst = ConcConst)


# visualization
library(rasterVis)
coul = brewer.pal(8, "YlGn")
coul = colorRampPalette(coul)(100)

b <- resultsMulti
# creating raster images
qr<-brick(b$fields[[6]][2:19,2:19,364:365])
levelplot(qr,main="Scenario Nr ",sub="day X and Y", col.regions = coul, ylab="Plant biomass density [g m-2]")

# saving data
save(b, file= "resultsMulti.RData")

#************************************************************
# ANALYSIS
#************************************************************
# analysing the data (just some examples)

# vertical water flux
flux <- b$fields[[13]][2:19,2:19,1:time]
mean(flux)
sd(flux)
sum(flux)

# Plant bimass density (overall)
mean(b$fields[[9]][2:19,2:19,1:time])
sd(b$fields[[9]][2:19,2:19,1:time])

# Plant bimass density Species A
mean(b$fields[[10]][2:19,2:19,1:time])
sd(b$fields[[10]][2:19,2:19,1:time])

# Plant bimass density Species B
mean(b$fields[[11]][2:19,2:19,1:time])
sd(b$fields[[11]][2:19,2:19,1:time])

# Plant bimass density Species C
mean(b$fields[[12]][2:19,2:19,1:time])
sd(b$fields[[12]][2:19,2:19,1:time])

# soil mositure
min(b$fields[[14]][2:19,2:19,1:time])
max((b$fields[[14]][2:19,2:19,1:time]))
mean(b$fields[[14]][2:19,2:19,1:time])

# attention. Indices and numbering in C++ and R differ. What is #0 in C++ is #1 in R.
#according to this:

# # in R it is:
# 1 = h;
# 2 = q;
# 3 = In;
# 4 = runon;
# 5 = Wu;
# 6 = WU_A;
# 7 = WU_B;
# 8 = WU_C;
# 9 = P;
# 10 = P_A;
# 11 = P_B;
# 12 = P_C;
# 13 = flux;
# 14 = M;
# 15 = SmM;
# 16 = CM;
# 17 = mb;
# 18 = Svir;
# 19 = SmI;
# 20 = Smh;
# 21 = Ch;
# 22 = qsdA;
# 23 = qsdB;
# 24 = qsdC;
# 25 = runonsdA;
# 26 = runonsdB;
# 27 = runonsdC;
# 28 = seep;
# 29 = salt_runon;
# 30 = mb_salt;

#**************************************************************************
#Patch Analysis

# prepare the matrix of biomass of the day that shall be considered for pattern analysis
storage <- a$fields[[6]][2:19,2:19,365]

storage <- as.matrix(storage)
# 
for (i in 1:nrow(storage)){
  
  for (j in 1:ncol(storage)){
    
    if(storage[i,j] < 1){
      
      storage[i,j] <- 0
    }
    else{
      storage[i,j] <- 1
    }
    
  }
}

#Run patch statistics

PatchAnalysis<- ClassStat(storage, cellsize = 1, bkgd = 0, latlon = FALSE)

#read out patch statistics
PatchAnalysis$n.patches
PatchAnalysis$total.area
PatchAnalysis$mean.patch.area
PatchAnalysis$sd.patch.area
PatchAnalysis$patch.density
PatchAnalysis$mean.shape.index
PatchAnalysis$sd.shape.index


# Do Diversity analysis (Shannon and Simpson Indices)
# install.packages("vegan")
library("vegan")
div_mat <- matrix(0, nrow= 1, ncol=3)

div_mat[1,1]<-sum(storage)
div_mat[1,2]<-sum(storage)
div_mat[1,3]<-sum(storage)


shannon <- diversity(div_mat, index = "shannon", MARGIN = 1, base = exp(1))
simpson <- diversity(div_mat, index = "simpson", MARGIN = 1)
shannon
simpson

# End ******************************************************************************************
# **********************************************************************************************

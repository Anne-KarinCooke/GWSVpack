### Script for simulation


# necessities:

# packages for Rtools, Rcpp, RcppArma...
# install.packages("devtools")
library("devtools")
find_rtools() # check wheter it finds Rtools, otherwise Rcpp won't run
# install.packages("Rcpp")
library("Rcpp")
# install.packages("RcppArmadillo")
library("RcppArmadillo")
# # packages for the patch and diversity analysis 
# library("SDMTools")
# library("vegan") #

# ****************************************************************************************************************************************
# Rainfall function to generate Poisson-distributed random rainfall series
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
sourceCpp("Mono_Species_Model.cpp", rebuild=T)

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
# ecosystem carrying capacity, [g/m2]
P0 <-  500.0
# seed dispersal and diffusion
c1 <- 2.25 # [1/mm] Saco and Moreno-de las Heras 2013
c02 <- 0.0002   #[m/d] tranformed to [mm/deltat] Saco and Moreno-de las Heras 2013
#seed diffusivity
Dp <- 0.3 #Saco and Moreno-de las Heras,2013
b1 <- 0.5 # relative importance of facilitation [0 -1]
b2 <- 0.5 # relative importance of competition [0 - 1]
Zr <- 300 #mm root depth
f <-1 # leaching factor [ 0 -1]


k = 12.0 #Saco et al, 2013
W0 = 0.2 #Saco et al, 2013
gmax = 0.05 #Saco et al, 2013
c = 10.0 #Saco et al, 2013
d = 0.02 #Saco et al, 2013 //fraction of plant mortality


# ****************************************************************************************************************************************
## SIMULATION

# declaring the input lists
dims = list(rows=rows,cols=cols,time=time, ext = ext)

fixedInput = list(deltat = deltat, alpha_i = alpha_i, cn = cn, Mn = Mn, P0 = P0, d=d,
                  c1 = c1, c02 = c02, Dp = Dp, zeta = zeta, f=f, Zr = Zr, b1 = b1, b2 = b2,
                  k=k, W0=W0, gmax=gmax, c=c, d=d)

# These simulation input parameters have been varied in the simulations

# groundwater depth
Z <- 3000 #mm
# water uptake half saturation constant
k1 <- 2.5 #mm
# facilitation range
q1 <- 0.1 #m
#salinity sensitivity
sigmaP <- 0.1 # L/g
# hillslope of model domain
gslp <- 0.01
# groundwater salt concentration
CMgw <- 3.0 #g/L
# salt concentration in rainfall
ConcConst <-  0.01 #g/L

# loading soil properties (Ksat, porosity etc), 
# other possiblities : L Med Clay Stony, S Clay Loam, Loamy Sand, H Clay, M Clay, C Sand

soilparS <- Soil_cpp("C Sand") 

simInput = list(Z =Z,k1 = k1,
                q1 = q1,sigmaP= sigmaP,
                gslp = gslp,CMgw = CMgw,
                ConcConst = ConcConst)

# calling the actual model
# adjust choice of "border" (outer cells that are not considered), input parameter lists for soil paramters, dimensions

#specific simulation parameters, rainfall input
results <- TheFunction1(border = 1, 
                         soilpar=soilparS, 
                         dims = dims,Rain = Rain2, 
                         fixedInput = fixedInput, 
                         simInput=simInput)


# visualization
library(rasterVis)
coul = brewer.pal(8, "YlGn")
coul = colorRampPalette(coul)(100)

a <- results
# creating raster images
qr<-brick(a$fields[[6]][2:19,2:19,364:365])
levelplot(qr,main="Scenario Nr ",sub="day X and Y", col.regions = coul, ylab="Plant biomass density [g m-2]")

# saving data
save(a, file= "results.RData")

#************************************************************************************************************************
# ANALYSIS
#************************************************************************************************************************
# analysing the data (just some examples)

# vertical water flux
flux <- a$fields[[7]][2:19,2:19,1:time]
mean(flux)
sd(flux)
sum(flux)

# Plant bimass density
mean(a$fields[[6]][2:19,2:19,1:time])
sd(a$fields[[6]][2:19,2:19,1:time])
max(a$fields[[6]][2:19,2:19,1:time])

# soil mositure
min(a$fields[[8]][2:19,2:19,1:time])
max((a$fields[[8]][2:19,2:19,1:time]))
mean(a$fields[[8]][2:19,2:19,1:time])

# attention. Indices and numbering in C++ and R differ. What is #0 in C++ is #1 in R.
#according to this:

# in C++ the number of the output is:

# f1( 0 ) = h;
# f1( 1 ) = q;
# f1( 2 ) = In;
# f1( 3 ) = runon;
# f1( 4 ) = Wu;
# f1( 5 ) = P;
# f1( 6 ) = flux;
# f1( 7 ) = M;
# f1( 8 ) = SmM;
# f1( 9 ) = CM;
# f1( 10 ) = mb;
# f1( 11 ) = Svir;
# f1( 12 ) = SmI;
# f1( 13 ) = Smh;
# f1( 14 ) = Ch;
# f1( 15 ) = qsd;
# f1( 16 ) = runonsd;
# f1( 17 ) = seep;
# f1( 18 ) = salt_runon;
# f1( 19 ) = mb_salt;


# in R it is:

# 1 = h;
# 2 = q;
# 3 = In;
# 4 = runon;
# 5 = Wu;
# 6 = P;
# 7 = flux;
# 8 = M;
# 9 = SmM;
# 10 = CM;
# 11 = mb;
# 12 = Svir;
# 13 = SmI;
# 14 = Smh;
# 15 = Ch;
# 16 = qsd;
# 17 = runonsd;
# 18 = seep;
# 19 = salt_runon;
# 20 = mb_salt;

#*****************************************************************************************************************
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


# End *************************************************************************************************************************************
# ****************************************************************************************************************************************

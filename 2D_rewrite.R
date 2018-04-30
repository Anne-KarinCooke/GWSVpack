## =======================================================================
## Ecohydroligical model
## dispersing in 2 dimensions
## =======================================================================

require("deSolve")
#Prerequisits
#Plant functions
WU <- function(M,P,gmax,k1){ 
  
  WU=gmax*(M/((M+k1)))*P 
  
  return(WU)
}
Gr <- function(M,P,cGr,gmax,k1) { 
  
  Gr = cGr*WU(M,P,gmax,k1)
  
  return(Gr)
}
# Mo <- function(P,M,Svir,pars) {
#   
#   Mo = P*(d*(M/Svir))
#   
#   return(Mo)
#   
# }
Mo <- function(P,d) {

  Mo = P*d

  return(Mo)

}

# water functions
Infil <- function(h,P,alpha_i,k, W0){
  I=alpha_i*h*(P+k*W0)/(P+k)
  return (I)
}

L_n <- function(M,Z,Zr,psi_s_bar,b,K_s) {
  hb <- - psi_s_bar*10^5 # (mm?)
  #soilpar$s_fc <- (Z/hb)^(-1/soilpar$b)
  s=M/(n*Zr)
  psi = hb*s^-b
  m=2 + 3/b#5.64  # in Salvucci's paper it is called n, but I called it m here to not confuse it with porosity
  q <-((Z/hb)^(-m)-(psi/hb)^(-m))/(1+(psi/hb)^(-m)+(m-1)*(Z/hb)^(-m))#/(soilpar$n*Zr) # # Salvucci
  # q<- (s^(2*soilpar$b+3))-(1+((3/2)/(m)))*(hb/(Z))^m  # Eagelson
  flux <- -K_s*q
  return(flux)
}


## ==================
## Model definitions
## ==================

lvmod2D <- function (time, state, pars, N, Da, dx) { # Rain, add that later
  NN <- N*N
  # Plants, P biomass density
  P <- matrix(nrow = N, ncol = N,state[1:NN])
  #Water , M for moisture in the soil
  # M <- matrix(nrow = N, ncol = N,state[(NN+1):(2*NN)]) # I did not get why there should be a second layer?
  M <- matrix(nrow = N, ncol = N,state[1:NN])
  
  h <- matrix(nrow = N, ncol = N,state[1:NN])

  
  with (as.list(pars), {
    
    # change of biomass density P
    dP <- Gr(M,P,cGr,gmax,k1) - Mo(P,d)
    #some simple uniform rain for now
    Prec <- 10 #  Rain[pmax(1,ceiling(t))]
    # change oh h (h = overland flow depth)
    dh <- Prec - Infil(h,P,alpha_i,k, W0)
    # infiltration coefficient adjustment
    alpha_i <- ifelse(h < K_s, 1, (1-(h-K_s)/h))
    # soil moisture
    dM <- Infil(h,P,alpha_i,k, W0) - WU(M,P,gmax,k1) + L_n(M,Z,Zr,psi_s_bar,b,K_s)
    
    zero <- rep(0, N)
    
    ## 1. Fluxes in x-direction; zero fluxes near boundaries
    FluxP<- -Da * rbind(zero,(P[2:N,] - P[1:(N-1),]), zero)/dx
    # FluxM <- -Da * rbind(zero,(M[2:N,] - M[1:(N-1),]), zero)/dx
    Fluxh <- -Da * rbind(zero,(h[2:N,] - h[1:(N-1),]), zero)/dx
    
    ## Add flux gradient to rate of change
    dP    <- dP - (FluxP[2:(N+1),] - FluxP[1:N,])/dx
    # dM    <- dM - (FluxM[2:(N+1),] - FluxM[1:N,])/dx
    dh    <- dh - (Fluxh[2:(N+1),] - Fluxh[1:N,])/dx
    
    ## 2. Fluxes in y-direction; zero fluxes near boundaries
    FluxP <- -Da * cbind(zero,(P[,2:N] - P[,1:(N-1)]), zero)/dx
    # FluxM <- -Da * cbind(zero,(M[,2:N] - M[,1:(N-1)]), zero)/dx # Da should be different for plants and water!!!
    Fluxh <- -Da * cbind(zero,(h[,2:N] - h[,1:(N-1)]), zero)/dx
    
    ## Add flux gradient to rate of change
    dP    <- dP - (FluxP[,2:(N+1)] - FluxP[,1:N])/dx
    # dM    <- dM - (FluxM[,2:(N+1)] - FluxM[,1:N])/dx
    dh    <- dh - (Fluxh[,2:(N+1)] - Fluxh[,1:N])/dx
    
    return(list(c(as.vector(dP), as.vector(dM)))) #, as.vector(dh)
  })
}


## ===================
## Model applications
## ===================
# pars    <- c(rIng   = 0.2,    # /day, rate of ingestion
#              rGrow  = 1.0,    # /day, growth rate of prey
#              rMort  = 0.2 ,   # /day, mortality rate of predator
#              assEff = 0.5,    # -, assimilation efficiency
#              K      = 5 ,      # mmol/m3, carrying capacity

pars    <- c(gmax = 0.05,
             k1 = 5.0,
             k =12,
             cGr = 10.0,
             d = 0.24,
             W0 = 0.2,
             alpha_i = 1.0,
             Zr = 400,
             n<-0.367, # porosity
             K_s<-52.08*10, # mm/day
             # campbell's b
             b<-6.4069, # neurotheta sandy clay loam;
             s_fc<-0.2677/n, # Field capacity
             # This is the bubbling pressure
             psi_s_bar<--1.2E-3, #
             h1bar =  -psi_s_bar,
             hb = psi_s_bar*-10^5, # mm
             Z = 1000 # [mm] actual groundwater depth 
             ) 

  Rain <- abs(rnorm(200))*2
  Rain <- ifelse(Rain<2,0,Rain)

R  <- 100                      # total length of surface, m 20
N  <- 100                      # number of boxes in one direction 50
dx <- R/N                     # thickness of each layer
Da <- 0.05                    # m2/d, dispersion coefficient

NN <- N*N                     # total number of boxes

## initial conditions
yini    <- rep(0, 2*N*N)
cc      <- c((NN/2):(NN/2+1)+N/2, (NN/2):(NN/2+1)-N/2)
yini[cc] <- yini[NN+cc] <- 50

# # Initial conditions: # ind/m2
# Water <- rep(1, times = numboxes)
# Water[25:40] <- 50
# state <- c(Water = Water) # initialise state variables

## solve model (5000 state variables...  use Cash-Karp Runge-Kutta method
times   <- seq(0, 100, by = 1)
out <- ode.2D(y = yini, times = times, func = lvmod2D, parms = pars,
              dimens = c(N, N), names = c("P", "M"),
              N = N, dx = dx, Da = Da, method = rkMethod("rk45ck")) # Rain=Rain, add that later

diagnostics(out)
summary(out)


# testing odesolve for water balance

# using example from Soetaert et al. 
# https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf
require(deSolve)

# example on page 3 APHIDS
Aphid <- function(t, APHIDS, parameters) {
  deltax <- c (0.5, rep(1, numboxes - 1), 0.5)
  Flux <- -D * diff(c(0, APHIDS, 0)) / deltax
  dAPHIDS <- -diff(Flux) / delx + APHIDS * r
  
    # the return value
    list(dAPHIDS )
 } # end

D <- 0.3 # m2/day diffusion rate
r <- 0.01 # /day net growth rate
delx <- 1 # m thickness of boxes
numboxes <- 60
# distance of boxes on plant, m, 1 m interval
Distance <- seq(from = 0.5, by = delx, length.out = numboxes)

# Initial conditions: # ind/m2
APHIDS <- rep(0, times = numboxes)
APHIDS[30:31] <- 1
state <- c(APHIDS = APHIDS) # initialise state variables

# Run for 200 timesteps
times <-seq(0, 200, by = 1)
 print(system.time(
  out <- ode.1D(state, times, Aphid, parms = 0, nspec = 1, names = "Aphid")
 ))
 
 head(out[,1:5])
 
 image(out, method = "filled.contour", grid = Distance,
       xlab = "time, days", ylab = "Distance on plant, m",
       main = "Aphid density on a row of plants")
 
 # Now rewrite this for the water balance
 # Lorenz function example
 WB <- function(t, Rain, state, parameters) {
   with(as.list(c(state, parameters)),{
     P <- Rain[pmax(1,ceiling(t))]
     #https://stackoverflow.com/questions/21557634/using-a-time-series-of-parameters-to-solve-ode-in-r
     # rate of change
       dWater <- -a*Water + P*b
         # return the rate of change
         list(c(dWater))
       }) 
   }
 
 times <- seq(0, 100, by = 0.01)
 a <- 0.2
 b <- 0.4
 parameters <- list(a = a, b = b)
 Rain <- abs(rnorm(100))*2
 Rain <- ifelse(Rain<2,0,Rain)
 
 state <- c(Water = 1)
 
 out <- ode(y = state, times = times, func = WB, parms = parameters, Rain = Rain)
 head(out)
 
 plot(out[,"time"],out[,"Water"])
 
 
 # Include Rain and ET
 
 WB <- function(t, Water, parameters, Rain) {
   P <- Rain[pmax(1,ceiling(t))]
   
   deltax <- c (0.5, rep(1, numboxes - 1), 0.5)
   Flux <- -D * diff(c(0, Water, 0)) / deltax
   dWater <- -diff(Flux) / delx + r*P - b*Water
   
   # the return value
   list(dWater )
 } # end
 
 # some random rainfall
 Rain <- abs(rnorm(200))*2
 Rain <- ifelse(Rain<2,0,Rain)
 
 D <- 0.3 # m2/day diffusion rate
 r <- 0.5 # infiltration rate
 b <- 0.15 # actualET rate
 delx <- 1 # m thickness of boxes
 numboxes <- 60
 Distance <- seq(from = 0.5, by = delx, length.out = numboxes)
 
 # Initial conditions: # ind/m2
 Water <- rep(1, times = numboxes)
 Water[25:40] <- 5
 state <- c(Water = Water) # initialise state variables
 
 times <-seq(0, 200, by = 1)
 print(system.time(
   out <- ode.1D(state, times, WB, parms = 0, Rain = Rain, nspec = 1, names = "Water")
 ))
 
 image(out, method = "filled.contour", grid = Distance,
       xlab = "time, days", ylab = "Distance, m",
       main = "Water")
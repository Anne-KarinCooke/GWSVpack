# testing odesolve for water balance

# using example from Soetaert et al.
# https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf
require(deSolve)
ode1    <- function (y, times, func, parms,
                     method = c("lsoda","lsode","lsodes","lsodar","vode","daspk",
                                "euler", "rk4", "ode23", "ode45", "radau",
                                "bdf", "bdf_d", "adams", "impAdams", "impAdams_d",
                                "iteration"),
                     ...)  {
  if (is.null(method)) method <- "lsoda"
  if (is.list(method)) {
    #  is() should work from R 2.7 on ...
    #   if (!is(method, "rkMethod"))
    if (!"rkMethod" %in% class(method))
      stop("'method' should be given as string or as a list of class 'rkMethod'")
    out <- rk(y, times, func, parms, method = method, ...)
  } else if (is.function(method))
    out <- method(y, times, func, parms,...)
  else if (is.complex(y))
    out <- switch(match.arg(method),
                  vode  = zvode(y, times, func, parms, ...),
                  bdf  = zvode(y, times, func, parms, mf = 22, ...),
                  bdf_d = zvode(y, times, func, parms, mf = 23, ...),
                  adams = zvode(y, times, func, parms, mf = 10, ...),
                  impAdams = zvode(y, times, func, parms, mf = 12, ...),
                  impAdams_d = zvode(y, times, func, parms, mf = 13, ...)
    )
  else
    out <- switch(match.arg(method),
                  lsoda = lsoda(y, times, func, parms, ...),
                  vode  = vode(y, times, func, parms, ...),
                  lsode = lsode(y, times, func, parms, ...),
                  lsodes= lsodes(y, times, func, parms, ...),
                  lsodar= lsodar(y, times, func, parms, ...),
                  daspk = daspk(y, times, func, parms, ...),
                  euler = rk(y, times, func, parms, method = "euler", ...),
                  rk4   = rk(y, times, func, parms, method = "rk4", ...),
                  ode23 = rk(y, times, func, parms, method = "ode23", ...),
                  ode45 = rk(y, times, func, parms, method = "ode45", ...),
                  radau = radau(y, times, func, parms, ...),
                  bdf  = lsode(y, times, func, parms, mf = 22, ...),
                  bdf_d = lsode(y, times, func, parms, mf = 23, ...),
                  adams = lsode(y, times, func, parms, mf = 10, ...),
                  impAdams = lsode(y, times, func, parms, mf = 12, ...),
                  impAdams_d = lsode(y, times, func, parms, mf = 13, ...),
                  iteration = iteration(y, times, func, parms, ...)
    )
  
  return(out)
}

ode.1DD    <- function (y, times, func, parms, nspec = NULL,
                       dimens = NULL, method = c("lsoda","lsode",
                                                 "lsodes","lsodar","vode","daspk",
                                                 "euler", "rk4", "ode23", "ode45","radau",
                                                 "bdf", "adams", "impAdams", "iteration"),
                       names = NULL, bandwidth = 1,
                       restructure = FALSE, ...){
  # check input
  if (is.character(method)) method <- match.arg(method)
  islsodes <- FALSE
  if (is.character(method))
    if (method=="lsodes") islsodes <- TRUE
  
  if (is.null(method)) method <- "lsoda"
  
  if (any(!is.na(pmatch(names(list(...)), "jacfunc"))))
    stop ("cannot run ode.1D with jacfunc specified - remove jacfunc from call list")
  
  if (is.null(nspec) && is.null(dimens))
    stop ("cannot run ode.1D: nspec OR dimens should be specified")
  
  #  if (islsodes && bandwidth != 1)
  #    stop ("cannot combine 'method = lsodes' with 'bandwidth' not = 1")
  
  iscomplex <- is.complex(y)
  
  N     <- length(y)
  if (is.null(nspec)  )
    nspec <- N/dimens
  if (N %% nspec != 0    )
    stop ("cannot run ode.1D: nspec is not an integer fraction of number of state variables")
  
  if (! is.null(names) && length(names) != nspec)
    stop("length of 'names' should equal 'nspec'")
  
  # Use ode.band if implicit method with nspec=1
  if (is.character(method))
    if( nspec == 1 & method %in% c("lsoda","lsode","lsodar","vode","daspk","radau")) {
      out <- ode.band(y, times, func, parms, nspec = nspec,
                      method = method, bandup = nspec * bandwidth,
                      banddown = nspec * bandwidth, ...)
      attr(out,"ynames") <- names
      if (is.null(dimens)) dimens <- N/nspec
      attr (out, "dimens") <- dimens
      attr (out, "nspec") <- nspec
      
      return(out)
    }
}

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
delx <- 10 # m thickness of boxes
numboxes <- 60
# distance of boxes on plant, m, 1 m interval
Distance <- seq(from = 0.5, by = delx, length.out = numboxes)

# Initial conditions: # ind/m2
APHIDS <- rep(0, times = numboxes)
APHIDS[20:31] <- 10 #  AC
state <- c(APHIDS = APHIDS) # initialise state variables

# Run for 200 timesteps
times <-seq(0, 100, by = 1) #  AC
 print(system.time(
  out <- ode.1DD(state, times, Aphid, parms = 0, nspec = 1, names = "Aphid")
 ))

 head(out[,1:5])

 image(out, method = "filled.contour", grid = Distance,
       xlab = "time, days", ylab = "Distance on plant, m",
       main = "Aphid density on a row of plants")
 
 # Now rewrite this for the water balance
 # Lorenz function example
 # WB <- function(t, Rain, state, parameters) {
 #   with(as.list(c(state, parameters)),{
 #     P <- Rain[pmax(1,ceiling(t))]
 #     #https://stackoverflow.com/questions/21557634/using-a-time-series-of-parameters-to-solve-ode-in-r
 #     # rate of change
 #       dWater <- -a*Water + P*b
 #         # return the rate of change
 #         list(c(dWater))
 #       })
 #   }
 # 
 # times <- seq(0, 100, by = 0.01)
 # a <- 0.2
 # b <- 0.4
 # parameters <- list(a = a, b = b)
 # Rain <- abs(rnorm(100))*2
 # Rain <- ifelse(Rain<2,0,Rain)
 # 
 # state <- c(Water = 1)
 # 
 # 
 # 
 # out <- ode(y = state, times = times, func = WB, parms = parameters, Rain = Rain)
 # head(out)
 # 
 # plot(out[,"time"],out[,"Water"])
 # 
 # 
 # # Include Rain and ET
 # 
 # WB <- function(t, Water, parameters, Rain) {
 #   P <- Rain[pmax(1,ceiling(t))]
 # 
 #   deltax <- c (0.5, rep(1, numboxes - 1), 0.5)
 #   Flux <- -D * diff(c(0, Water, 0)) / deltax
 #   dWater <- -diff(Flux) / delx + r*P - b*Water
 # 
 #   # the return value
 #   list(dWater )
 # } # end
 # 
 # # some random rainfall
 # Rain <- abs(rnorm(200))*2
 # Rain <- ifelse(Rain<2,0,Rain)
 # 
 # D <- 0.3 # m2/day diffusion rate
 # r <- 0.5 # infiltration rate
 # b <- 0.15 # actualET rate
 # delx <- 1 # m thickness of boxes
 # numboxes <- 60
 # Distance <- seq(from = 0.5, by = delx, length.out = numboxes)
 # 
 # # Initial conditions: # ind/m2
 # Water <- rep(1, times = numboxes)
 # Water[25:40] <- 50
 # state <- c(Water = Water) # initialise state variables
 # 
 # times <-seq(0, 200, by = 1)
 # print(system.time(
 #   out <- ode.1D(state, times, WB, parms = 0, Rain = Rain, nspec = 1, names = "Water")
 # ))
 # 
 # image(out, method = "filled.contour", grid = Distance,
 #       xlab = "time, days", ylab = "Distance, m",
 #       main = "Water")
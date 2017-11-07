#*******************************************************************************************************************************
# DEM generation

#for random field DEM generation (surface elevation variance and corr.length)
# install.packages('geoR')
library("geoR")

# Elev1
A <- c(0.0001, 5) # [m], standard deviation: 1 cm, correlation length: 5m
# generating a Gaussian random field with geoR-package function grf()
# grid 5 times finer than needed
sim1 <-grf((((rows*5)*(cols*5))), grid = "reg", cov.pars = A)
elev_data1 <- matrix(sim1$data, nrow = (rows*5),ncol=(cols*5))

# Elev2
B <- c(0.01, 20) # m, standard deviation: 10 cm, correlation length:  20 m
# generating a Gaussian random field with geoR-package function grf()
# grid 5 times finer than needed
sim2 <-grf((((rows*5)*(cols*5))), grid = "reg", cov.pars = B)
# Taking every 5th element
elev_data2 <- matrix(sim2$data, nrow = (rows*5),ncol=(cols*5))

image(sim2)
# Taking every 5th element
elev_data1 <- elev_data1[seq(1,nrow(elev_data1),5),seq(1,ncol(elev_data1),5)]
elev_data2<- elev_data2[seq(1,nrow(elev_data2),5),seq(1,ncol(elev_data2),5)]

# writing the tables that will be used as input in model
write.table(elev_data1, 'B1.txt', sep=',')
write.table(elev_data2, 'B2.txt', sep=',')
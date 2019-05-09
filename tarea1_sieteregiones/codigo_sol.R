library(sp)
library(spData)

library(spgwr)
library(adespatial)
library(raster)
library(tmap)
library(lmtest)
library(RColorBrewer)
library(classInt)
library(spdep)
library(sphet)
library(pgirmess)
library(olsrr)
library(rgdal)
library(spatialreg)

library(Matrix)
library(maptools)
library(lattice)
library(rgdal)
library(ade4)
library(ggplot2)



#Ingresar matriz de pesos


R1 <- c(0	,1,	0,	0,	0,	0,	0)
R2 <- c(1	,0,	1,	0	,0	,0	,0)
R3 <- c(0,	1,	0,	1,	0,	0,	0)
R4 <- c(0	,0,	1,	0,	1,	0,	0)
R5 <- c(0,	0,	0,	1,	0,	1,	0)
R6 <- c(0,	0,	0,	0,	1,	0,	1)
R7  <- c(0,	0,	0	,0	,0	,1,	0)
matrizpesos <- matrix(scan(), ncol=7)







### Cargar Librerias 

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


## 


col.poly <- readOGR("/home/viviana/Escritorio/UD/Geoestadistica/EJERCICIOS y R/localidades_bogota/localidades.shp") 
coords <- coordinates(col.poly)
summary(col.poly)

#Reina1
nb.FOQ = poly2nb(col.poly, queen=TRUE, row.names=col.poly$DESEMPLEO)
plot(col.poly, border="black", main="Reina orden 1")
plot(nb.FOQ, coords, add=TRUE, col="green")

#Torre1
nb.RK = poly2nb(col.poly, queen=FALSE, row.names=col.poly$DESEMPLEO)
plot(col.poly, border="black", main="Torre orden 1")
plot(nb.RK, coords, add=TRUE, col="red")

#Reina2
o.nbQ2 <- read.gal("/home/viviana/Escritorio/UD/Geoestadistica/EJERCICIOS y R/Q2localidades.gal") 
plot(col.poly, border="black", main="Reina orden 2")
plot(o.nbQ2, coords, add=TRUE, col="green")

#I de Moran Global -- FALTA

#I de Moran Local
wQ1 <-  nb2listw(o.nbQ2, style='B')
moran.test(col.poly$DESEMPLEO,wQ1, randomisation=FALSE)
lmoran <- localmoran(col.poly$DESEMPLEO,wQ1)
print(lmoran)


#I de Moran Bivariado --FALTA
library(spdep)
datosMoranB <- readOGR("/home/viviana/Escritorio/UD/Geoestadistica/EJERCICIOS y R/parcial1_localidades/localidades_bogota/localidades.shp") 
localidades <- poly2nb(datosMoranB)
a.lw <- nb2listw(localidades, style="W")
W <- as.matrix(as_dgRMatrix_listw(a.lw))

moran.bi(datosMoranB$DESEMPLEO,datosMoranB$ROBOS,a.lw,zero.policy =T)



source("/home/viviana/Escritorio/UD/Geoestadistica/EJERCICIOS y R/Funciones/moran.bi.R")
W <- as.matrix(as_dgRMatrix_listw(a.lw))
moran.bi(columbus$CRIME,columbus$INC,a.lw,zero.policy =T)


#Dispersograma 

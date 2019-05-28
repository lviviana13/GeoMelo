#######################################################################
############             ############
#######################################################################

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

#Cargar Shape File
ohlung <- readOGR("/home/viviana/Escritorio/UD/Geoestadistica/GeoMelo/PROYECTO/datos/ohiolung/ohlung.shp")
coords <- coordinates(ohlung)

theme_set(theme_update())
qplot(ohlung$POPMW88,geom="histogram", stat="count",col=I("cYan"),fill=I("cyan"), alpha=I(0.4), main="", xlab="LMW88C",ylab="Frecuencia")

###Esto Corre pero no sé que hace
rn <- sapply(slot(ohlung, "polygons"), function(x) slot(x, "ID"))
k1 <- knn2nb(knearneigh(coords))
all.linked <- max(unlist(nbdists(k1, coords)))
col.nb.0.all <- dnearneigh(coords, 0, all.linked, row.names=rn)
x11()
summary(col.nb.0.all, coords)
plot(ohlung, border="blue")
plot(col.nb.0.all, coords, add=TRUE)
title(main=paste("Distance based neighbours 0-", format(all.linked), " distance units", sep=""))


# Identificaci?n vecinos
par(mai=c(0,0,0,0))
set.seed(1)
plot(ohlung, col=sample(rainbow(500))[1:49])
xy <- coordinates(ohlung)
points(xy, cex=3, pch=20, col='white')
text(ohlung, "POLYID", cex=0.6)
detach(package:raster)


#######################################################
# 1. ACCESO A LOS DATOS Y CALCULO PESOS ESPACIALES
#######################################################
col_nbq1 <- poly2nb(ohlung) 
col.lags <- nblag(col_nbq1, 9)

col.poly <- ohlung

# col.poly <- readShapePoly("D:/GeoDA/columbus.shp")       # Seleccionar el directorio previamente de GeoDA
brks <- quantile(col.poly$CRIME, seq(0, 1, 1/max(col.poly$CRIME)))
cols <- grey((length(brks):2)/length(brks))
# my.colors <- rev(c("red", "green", "blue","limegreen","brown","darkgrey","yellow","purple","orange","cyan","pink", "white"))    
# asigna colores diferentes, pero sin quantil
l2 = list("SpatialPolygonsRescale", layout.north.arrow(), offset =  c(6.1,13.5), scale = 0.5)   # Layout
l3 = list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(6.6,13.5),  scale = 0.5, fill=c("transparent","black"))
l4 = list("sp.text", c(6.6,13.7), "0")
l5 = list("sp.text", c(7.1,13.7), "500 m")
spplot(col.poly["CRIME"], scales=list(draw=TRUE), col.regions=cm.colors(20),  sp.layout=list(l2,l3,l4,l5))
# rainbow(n, start=.7, end=.1), heat.colors(n), terrain.colors(n), topo.colors(n), cm.colors(n)





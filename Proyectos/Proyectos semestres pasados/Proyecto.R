library(zoo)
library(splines)
library(carData)
library(car)
library(sandwich)
library(RcmdrMisc)
library(sp)
library(rgdal)
library(geoR)
install.packages("stringi")
library(maptools)
library(Matrix)
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
library(rCarto)
library(ggplot2)

source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/northarrow.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/scalebar.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/moran.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/geary.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/moranbi.test.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/randomize_vector.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/moran.cluster.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/moran.bi.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/moran.cluster.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/getis.cluster.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/localmoran.bi.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/moranbi.plot.R") 
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/quantile.e.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/sp.na.omit.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/correlogram.d.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/sp.correlogram.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/spcorrelogram.bi.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/moran1.bi.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/moranbi1.test.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/geary.bi.R")
source("C:/Users/eljul/OneDrive/Documents/R/Funciones-20181009T202642Z-001/Funciones/test.W.R")

#BASE DE DATOS
setwd("C:/Users/eljul/OneDrive/Desktop/GEOESTADISTICA/PROYECTO/houston_hom")
houston <-readShapePoly("C:/Users/eljul/OneDrive/Desktop/GEOESTADISTICA/PROYECTO/houston_hom/hou_hom.shp")
x11()
BASE_DATOS <- data.frame(houston$HR8187,houston$HR8791,houston$HR9195,houston$HC8187,houston$HC8791,houston$HC9195,houston$PO8187,houston$PO8791,houston$PO9195,houston$PE77,houston$PE82,houston$PE87,houston$RDAC80,houston$RDAC85,houston$RDAC90)
names(BASE_DATOS)<-c("Homicidio 8187", "Homicidio 8791", "Homicidio 9195", "Cantidad de hom 8187", "Cantidad de hom 8791", "Cantidad de hom 9195", "Poblacion 8187", "Poblacion 8791","Poblacion 9195","Gastos en policia 77", "Gastos en policia 82","Gastos en policia 87","Privacion de rescursos 80", "Privacion de rescursos 85", "Privacion de rescursos 90")

#DESCRIPCION DE LAS VARIABLES
summary(houston$HR9195)
sd(houston$HR9195)
summary(houston$HR8187)
sd(houston$HR8187)
summary(houston$HR8791)
sd(houston$HR8791)
summary(houston$HC8187)
sd(houston$HC8187)
summary(houston$HC8791)
sd(houston$HC8791)
summary(houston$HC9195)
sd(houston$HC9195)
summary(houston$PO8791)
sd(houston$PO8791)
summary(houston$PO8187)
sd(houston$PO8187)
summary(houston$PO9195)
sd(houston$PO9195)
summary(houston$PE77)
sd(houston$PE77)
summary(houston$PE82)
sd(houston$PE82)
summary(houston$PE87)
sd(houston$PE87)
summary(houston$RDAC80)
sd(houston$RDAC80)
summary(houston$RDAC85)
sd(houston$RDAC85)
summary(houston$RDAC90)
sd(houston$RDAC90)

boxplot(houston$HR8187, col="blue")
title(main="Diagrama de caja (HR8187)")
shapiro.test(t(houston$HR8187))
#Qplot variable respuesta
qqnorm(houston@data$HR8187, col="blue")
qqline(houston@data$HR8187)

#diagrama de dispersion
ggplot(houston@data, aes(x=HR8187, y=HR9195)) +   geom_point()
ggplot(houston@data, aes(x=HR8791, y=HR9195)) +   geom_point()
ggplot(houston@data, aes(x=HC8187, y=HC9195)) +   geom_point()
ggplot(houston@data, aes(x=HC8791, y=HC9195)) +   geom_point()
ggplot(houston@data, aes(x=PO8187, y=PO9195)) +   geom_point()
ggplot(houston@data, aes(x=PO8791, y=PO9195)) +   geom_point()
ggplot(houston@data, aes(x=PE77, y=PE87)) +   geom_point()
ggplot(houston@data, aes(x=PE82, y=PE87)) +   geom_point()
ggplot(houston@data, aes(x=RDAC80, y=RDAC90)) +   geom_point()
ggplot(houston@data, aes(x=RDAC85, y=RDAC90)) +   geom_point()

#########PESOS##############
#CONTIGUIDAD EFECTO TORRE
coordenadas <- coordinates(houston)
houston_nbr <- poly2nb(houston, queen=F)
#EFECTO REINA
houston_nbq1 <- poly2nb(houston)                    #ORDEN 1
plot(houston)
plot(houston_nbq1, coordenadas, add=T, col=4, lwd=2)
title(main="Efecto reina orden 1")
houston_nbq2 <- nblag(houston_nbq1, 2)                  # Orden 2
houston_nbq3 <- nblag(houston_nbq1, 3)                 # Orden 3
houston_nbq4 <- nblag(houston_nbq1, 4)                 # Orden 4
houston_nbq5 <- nblag(houston_nbq1, 5)                 # Orden 5
#EFECTO REINA ORDEN 2
plot(houston_nbq2[2], coordenadas, add=T, col=2, lwd=2)
plot(houston)
#CORRELOGRAMA DE MORAN 
sp.cr <- sp.correlogram(houston_nbq1, houston@data$HR8187, order=5, method="corr", style="W", zero.policy=T) 
cor.s <- sp.correlogram(houston_nbq1, houston@data$HR8187, order=5, method="I", style="W", zero.policy=T)    
plot(cor.s)
x11()
#K VECINOS MAS CERCANOS 
IDs <- row.names(as(houston, "data.frame"))
hous_kn1<-knn2nb(knearneigh(coordenadas, k=1), row.names=IDs)
hous_kn2<-knn2nb(knearneigh(coordenadas, k=2), row.names=IDs)
hous_kn3<-knn2nb(knearneigh(coordenadas, k=3), row.names=IDs)
hous_kn4<-knn2nb(knearneigh(coordenadas, k=4), row.names=IDs)
hous_kn5<-knn2nb(knearneigh(coordenadas, k=5), row.names=IDs)

par(mar = c(0, 0, 0, 0), pty = "s")
plot(houston)
plot(hous_kn1, coordenadas, add=T, col=2, lwd=2)
plot(houston)
plot(hous_kn5, coordenadas, add=T, col=3, lwd=2)
par(mar = c(0, 0, 0, 0), pty = "s")
plot(houston)
plot(hous_kn4, coordenadas, add=T, col=2, lwd=2)

# CORRELOGRAMA MORAN K VECINOS 
IDs <- row.names(houston@data)                                                                
hous.kn5<-knn2nb(knearneigh(coordenadas, k=5), row.names=IDs)                                                       
sp.cr <- sp.correlogram(hous.kn5, houston$HR8187, order=5, method="corr", style="W", zero.policy=T) 
cor <- sp.correlogram(hous.kn5, houston$HR8187, order=5, method="I", style="W", zero.policy=T)      
plot(cor)

#DISTANCIA DE GABRIEL
gabrielnb=graph2nb(gabrielneigh(coordenadas),sym=TRUE)
plot(houston,border="black")
plot(gabrielnb,coordenadas,add=T,col="purple")
title(main="Gráfica de Gabriel")

#DISTANCIA DE DELAUNAY
trinb=tri2nb(coordenadas)
plot(houston,border="gray")
plot(trinb,coordenadas,add=T,col="yellow")
title(main="Triangulación Delaunay")

#MATRICES DE VECINDAD
hous_nbq1 <- poly2nb(houston)               # Efecto Reina
colt.lags <- nblag(hous_nbq1, 10)
q.lw1  <- nb2listw(colt.lags[[1]], style="W",zero.policy =T)
q.lw2  <- nb2listw(colt.lags[[2]], style="W",zero.policy =T)
q.lw3  <- nb2listw(colt.lags[[3]], style="W",zero.policy =T)
q.lw4  <- nb2listw(colt.lags[[4]], style="W",zero.policy =T)
q.lw5  <- nb2listw(colt.lags[[5]], style="W",zero.policy =T)

hous_kn1<-knn2nb(knearneigh(coordenadas, k=1), row.names=IDs)
hous_kn2<-knn2nb(knearneigh(coordenadas, k=2), row.names=IDs)
hous_kn3<-knn2nb(knearneigh(coordenadas, k=3), row.names=IDs)
hous_kn4<-knn2nb(knearneigh(coordenadas, k=4), row.names=IDs)
hous_kn5<-knn2nb(knearneigh(coordenadas, k=5), row.names=IDs)

#I-MORAN  , Z.I  Y P_value  para cada peso espacial
HR8187 <- as.data.frame(houston)$HR8187
y <- c(HR8187)

a.lwR1 <- nb2listw(hous_nbq1, style="W",zero.policy =T)
a.lwbR1 <- nb2listw(hous_nbq1, style="B",zero.policy = T)
moran.test(y, a.lwR1, zero.policy=T)

a.lwR2 <- nb2listw(colt.lags[[2]], style="W",zero.policy =T)
a.lwbR2 <- nb2listw(colt.lags[[2]], style="B",zero.policy = T)
moran.test(y, a.lwR2, zero.policy=T)

a.lwR3 <- nb2listw(colt.lags[[3]], style="W",zero.policy =T)
a.lwbR3 <- nb2listw(colt.lags[[3]], style="B",zero.policy = T)
moran.test(y, a.lwR3, zero.policy=T)

a.lwR4 <- nb2listw(colt.lags[[4]], style="W",zero.policy =T)
a.lwbR4 <- nb2listw(colt.lags[[4]], style="B",zero.policy = T)
moran.test(y, a.lwR4, zero.policy=T)

a.lwR5 <- nb2listw(colt.lags[[5]], style="W",zero.policy =T)
a.lwbR5 <- nb2listw(colt.lags[[5]], style="B",zero.policy = T)
moran.test(y, a.lwR5, zero.policy=T)


a.lwK1 <- nb2listw(hous_kn1, style="W",zero.policy =T)
a.lwbK1 <- nb2listw(hous_kn1, style="B",zero.policy = T)
moran.test(y, a.lwK1, zero.policy=T)

a.lwK2 <- nb2listw(hous_kn2, style="W",zero.policy =T)
a.lwbK2 <- nb2listw(hous_kn2, style="B",zero.policy = T)
moran.test(y, a.lwK2, zero.policy=T)

a.lwK3 <- nb2listw(hous_kn3, style="W",zero.policy =T)
a.lwbK3 <- nb2listw(hous_kn3, style="B",zero.policy = T)
moran.test(y, a.lwK3, zero.policy=T)

a.lwK4 <- nb2listw(hous_kn4, style="W",zero.policy =T)
a.lwbK4 <- nb2listw(hous_kn4, style="B",zero.policy = T)
moran.test(y, a.lwK4, zero.policy=T)

a.lwK5 <- nb2listw(hous_kn5, style="W",zero.policy =T)
a.lwbK5 <- nb2listw(hous_kn5, style="B",zero.policy = T)
moran.test(y, a.lwK5, zero.policy=T)

a.lwDG <- nb2listw(gabrielnb, style="W",zero.policy =T)
a.lwbDG <- nb2listw(gabrielnb, style="B",zero.policy = T)
moran.test(y, a.lwDG, zero.policy=T)

a.lwTN <- nb2listw(trinb, style="W",zero.policy =T)
a.lwbTN <- nb2listw(trinb, style="B",zero.policy = T)
moran.test(y, a.lwTN, zero.policy=T)

#AIC MEJOR MATRIZ DE PESOS
col_nbr1 <- poly2nb(houston,queen=FALSE)
colr.lags <- nblag(col_nbr1, 10)
r.lw1  <- nb2listw(colr.lags[[1]], style="W",zero.policy =T)
r.lw2  <- nb2listw(colr.lags[[2]], style="W",zero.policy =T)

summary(test.W(houston@data$HR8187,colr.lags[[2]],coordenadas))
summary(test.W(houston@data$HR8187,colr.lags[[3]],coordenadas))
summary(test.W(houston@data$HR8187,colr.lags[[4]],coordenadas))
summary(test.W(houston@data$HR8187,colr.lags[[5]],coordenadas))
summary(test.W(houston@data$HR8187,hous_nbq1))
summary(test.W(houston@data$HR8187,trinb))
summary(test.W(houston@data$HR8187,gabrielnb))
summary(test.W(houston@data$HR8187,hous_kn1))
summary(test.W(houston@data$HR8187,hous_kn2))
summary(test.W(houston@data$HR8187,hous_kn3))
summary(test.W(houston@data$HR8187,hous_kn4))
summary(test.W(houston@data$HR8187,hous_kn5))
detach(package:adespatial)

#MAPAS DE LOS VECINOS 
op=par(mfrow=c(2,2))

plot(houston)
plot(houston_nbq5[[5]], coordenadas, add=T, col="cyan", lwd=2)
title(main="Reina Orden 3")

plot(houston)
plot(hous_kn5, coordenadas, add=T, col="orange", lwd=2)
title(main="K-vecinos K-5")

plot(houston,border="gray")
plot(trinb,coordenadas,add=T,col="green")
title(main="Triangulación Delaunay")

plot(houston,border="gray")
plot(gabrielnb,coordenadas,add=T,col="red")
title(main="Gráfica de Gabriel")
par(op)

# Correlograma Moran a partir de matriz contiguidad espacial                                                  
sp.cr <- sp.correlogram(a.lwTN, houston@data$HR8187, order=5, method="corr", style="W", zero.policy=T) 
cor.s <- sp.correlogram(a.lwTN, houston@data$HR8187, order=5, method="I", style="W", zero.policy=T)    
plot(cor.s)

######ANALISIS EXPLORATORIO##########
plotvar <- round(houston@data$HR8187,1)    
nclr <- 5
plotclr <- brewer.pal(nclr,"YlOrBr")
class <- classIntervals(plotvar, nclr, style="quantile")
houcode <- findColours(class, plotclr)
plot(houston)   
plot(houston, col=houcode, add=T)
title(main="Tasa de Homicidios", sub="Cuantiles (Igual-Frecuencia) Intervalos de Clase")
mtext("Cuantiles (Igual-Frecuencia) Intervalos de Clase", side=1)
legend(locator(1), legend=names(attr(houcode, "table")), fill=attr(houcode, "palette"), cex=0.8, bty="n")
northarrow(c(10.1,14.5),0.2,cex=0.6)
scalebar(c(6.1,10.7),length=1,unit="m",division.cex=0.6)

###I MORAN Y G GETIS GLOBALES VAR.DEPENDIENTE
moran.test(y, a.lwTN , zero.policy=T)
library(xtable)
library(stargazer)

HR8187<-houston@data$HR8187
W <- as.matrix(as_dgRMatrix_listw(a.lwTN))
Y <- HR8187
WY <- W%*%Y
WY1<-(Y-mean(Y))/sd(Y)
X11()
moranbi.plot(HR8187,WY1,quiet =F,zero.policy =F,listw=a.lwTN)
title(main="Grafico de dispersión de Moran")

# Getis-Ord global G statistic
globalG.test(y, listw=a.lwTN)

##I MORAN Y G GETIS GLOBALES VAR.INDEPENDIENTES ##
moran.test(houston$HC8187, a.lwTN , zero.policy=T)
moran.test(houston$PO8187, a.lwTN, zero.policy=T)
moran.test(houston$PE77, a.lwTN, zero.policy=T)
moran.test(houston$RDAC80, a.lwTN, zero.policy=T)

#I DE MORAN LOCALES 
moranbi1.test(x=houston@data$HR8187,y=houston@data$HC8187,a.lwTN,zero.policy =T,randomisation =T,
              alternative="two.sided",adjust.n=TRUE)
moranbi.plot(HR8187,houston@data$HC8187,quiet =F,zero.policy =F,listw=a.lwTN)
set.seed(123)
MBCH1 <- moranbi.test(HR8187,houston@data$HC8187,a.lwTN,N=999,zero.policy =T,graph=T)

moranbi1.test(x=houston@data$HR8187,y=houston@data$PO8187,a.lwTN,zero.policy =T,randomisation =T,
              alternative="two.sided",adjust.n=TRUE)
moranbi.plot(HR8187,houston@data$PO8187,quiet =F,zero.policy =F,listw=a.lwTN)

moranbi1.test(x=houston@data$HR8187,y=houston@data$PE77,a.lwTN,zero.policy =T,randomisation =T,
              alternative="two.sided",adjust.n=TRUE)
moranbi.plot(HR8187,houston@data$PE77,quiet =F,zero.policy =F,listw=a.lwTN)

moranbi1.test(x=houston@data$HR8187,y=houston@data$RDAC80,a.lwTN,zero.policy =T,randomisation =T,
              alternative="two.sided",adjust.n=TRUE)
moranbi.plot(HR8187,houston@data$RDAC80,quiet =F,zero.policy =F,listw=a.lwTN)

# LISA Cluster Map
moran.cluster(houston@data$HR8187, a.lwTN, zero.policy = T, houston, significant=T)
moran.cluster(houston@data$HC8187, a.lwTN, zero.policy = T, houston, significant=T)
moran.cluster(houston@data$PO8187, a.lwTN, zero.policy = T, houston, significant=T)
moran.cluster(houston@data$PE77, a.lwTN, zero.policy = T, houston, significant=T)
moran.cluster(houston@data$RDAC80, a.lwTN, zero.policy = T, houston, significant=T)

LMCI <- localmoran.bi(houston@data$HR8187, houston@data$RDAC80, a.lwTN, zero.policy =T)
LMCH <- localmoran.bi(houston@data$HR8187, houston@data$HC8187, a.lwTN, zero.policy =T)
LMCJ<- localmoran.bi(houston@data$HR8187, houston@data$PO8187, a.lwTN, zero.policy =T)
LMCK <- localmoran.bi(houston@data$HR8187, houston@data$PE77, a.lwTN, zero.policy =T)
xtable(LMCI)

# Getis Cluster Map
getis.cluster(houston@data$HR8187, a.lwTN, zero.policy = T, houston, significant=T)
X11()

#  Mapeando outliers e influyentes locales   ##
x1<-(y-mean(y))/sd(y)
mp<-moran.plot(x1, a.lwTN, main="Gráfico de Dispersión de Moran")  
# Cambiar por "x1", para la estandarización
infl1 <- apply(mp$is.inf, 1, any)
lhx1 <- cut(houston@data$HR8187, breaks=c(min(houston@data$HR8187), mean(houston@data$HR8187), 
                                       max(houston@data$HR8187)), labels=c("L", "H"), OWNHlude.lowest=TRUE)
wx1 <- lag(a.lwTN, houston$HR8187)
lhwx1 <- cut(wx1, breaks=c(min(wx1), mean(wx1), max(wx1)), labels=c("L", "H"), OWNHlude.lowest=TRUE)
lhlh1 <- interaction(lhx1, lhwx1, infl1, drop=TRUE)
cols1 <- rep(1, length(lhlh1))
cols1[lhlh1 == "L.L.TRUE"] <- 2
cols1[lhlh1 == "L.H.TRUE"] <- 3
cols1[lhlh1 == "H.L.TRUE"] <- 4
cols1[lhlh1 == "H.H.TRUE"] <- 5
plot(houston, col=houcode[cols1])          # gray.colors(4, 0.95, 0.55, 2.2)[cols])
legend(locator(1), legend=c("None", "LL", "LH", "HL", "HH"), fill=attr(houcode, "palette"), bty="n", cex=1, y.intersp=0.8)
title(main="OUTLIERS", sub="Cuantiles (Igual-Frecuencia) Intervalos de Clase")

#Ajuste Modelo Clásico

HR8187<-houston$HR8187
HC8187<-houston$HC8187
PO8187<-houston$PO8187
PE77<-houston$PE77
RDAC80<-houston$RDAC80
X<- c(-95.16776,-95.65264,-94.61727,-94.17066,-95.99597,-93.85106,-95.42231,-94.61541,-95.13589,-96.51482,-93.74313,-94.02381,-94.83122,-96.97963,-95.92894,-94.37574,-95.57348,-96.30147,-97.60054,-95.16754,-95.98265,-96.62163,-95.50280,-96.96913,-94.39044,-93.35198,-92.80772,-94.81729,-97.31930,-96.40587,-95.99345,-93.89652,-94.17504,-96.92798,-95.39997,-92.31682,-96.28334,-97.63011,-93.15898,-96.53472,-94.60937,-95.77508,-97.49746,-96.93384,-96.22918,-95.45939,-95.05608,-97.36174,-96.57651,-95.98546,-96.97217,-96.67527)
Y<- c(31.84179,31.81864,31.61987,31.39499,31.30229,31.34771,31.32299,31.25833,31.09167,31.3619,30.78774,30.74747,30.79483,30.79765,30.97228,30.77024,30.74642,30.66814,30.65359,30.58425,30.54904,30.49999,30.30396,30.31724,30.33164,30.23505,30.27607,30.15462,30.10444,30.21878,30.01187,30.12769,29.88312,29.88161,29.86115,29.85527,29.88682,29.84065,29.88179,29.62126,29.74417,29.53185,29.46171,29.38628,29.27753,29.19447,29.41036,29.08584,28.94879,28.86984,28.79591,28.50108)
Xo<-X
Yo<-Y
houston$WXX <- lag.listw(a.lwTN, Xo)
houston$WYY <- lag.listw(a.lwTN, Yo)
WX<-houston$WXX 
WY<-houston$WYY

model1 <- lm(HR8187 ~ HC8187 + PO8187+PE77+RDAC80 , data=houston)
summary(model1)
houston$lm_resi <- residuals(model1)
x11()
#MAPA RESIDUALES 
spplot(houston["lm_resi"], col.regions = rev(terrain.colors(20)))
lm.morantest(model1, a.lwTN, alternative = "two.sided")
Z.Ir <- (0.110059425-(-0.026377907))/sqrt(0.006038176) 
acf(houston$lm_resi)
pacf(houston$lm_resi)# Validación supuestos
bptest(model1)
resettest(model1)
raintest(model1)
shapiro.test(residuals(model1))
vif(model1)
AIC(model1)

#MULTIPLICADORES DE LAGRANGE
errorsarlm.morantest(col.error.sm, a.lwTN)
lm.LMtests(model1,a.lwTN,test = c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA"))

#3. AJUSTE DEL MODELO SAR 
datos.hou <- as.data.frame(houston)
d.hou <- as.data.frame(cbind(datos.hou$CNTY_FIPS, datos.hou$HC8187, datos.hou$PO8187, datos.hou$PE77, datos.hou$RDAC80, datos.hou$HR8187))
names(d.hou) <- c("CODIGO", "HC8187", "PO8187", "PE77", "RDAC80", "HR8187")

attach(d.hou)
Y <- HR8187
X <- cbind(HC8187, PO8187,PE77, RDAC80)
X <- cbind(1, X)
X[1:5,]
cardnb <- card(a.lwTN$neighbours)                   # cuenta los vecinos de cada poligono
str(sd1 <- 1/sqrt(cardnb))

glist <- vector(mode = "list", length = length(a.lwTN$neighbours))
for (i in 1:length(a.lwTN$neighbours)) glist[[i]] <- cardnb[i]*a.lwTN$weights[[i]]
for (i in 1:length(a.lwTN$neighbours)) {
  inb <- a.lwTN$neighbours[[i]]
  icd <- cardnb[i]
  if (icd > 0) {
    for (j in 1:icd) {
      glist[[i]][j] <- sd1[i] * glist[[i]][j] * sd1[inb[j]]
    }
  }
}

W <- as_dgRMatrix_listw(similar.listw(a.lwTN))
str(W)
W@X
W <- as.matrix(W)

I <- as_dsCMatrix_I(nrow(X))
# I <- as_dsCMatrix_I(dim(W)[1])
I@x                                              # Slot "x":
I@i                                              # Slot "i":
I@p                                              # Slot "p":
I@Dim                                            # Slot "Dim":
I <- as.matrix(I)

lm.null <- lm(Y ~ X)                         # Modelo sin intercepto
lm.w <- lm.fit(X, WY)
e.null <- lm.null$residuals
e.w <- lm.w$residuals
e.a <- t(e.null) %*% e.null
e.b <- t(e.w) %*% e.null
e.c <- t(e.w) %*% e.w

sar.lag.mix.f.sM <- function(rho, W, I, e.a, e.b, e.c, n, tmpmax){
  SSE <- e.a - 2*rho*e.b + rho*rho*e.c
  s2 <- SSE/n
  Jacobian <- log(det(chol((I - rho * W), LINPACK=TRUE))^2)
  # tmpmax=(3*((sum(card(a.lw$neighbours))) + n))))^2)
  ret <- (Jacobian - ((n/2)*log(2*pi)) - (n/2)*log(s2) - 
            (1/(2*s2))*SSE)
  ret
}

opt <- optimize(sar.lag.mix.f.sM, interval = c(-1,1), maximum =  TRUE, 
                tol = .Machine$double.eps^0.5, W = W, I = I, e.a = e.a, e.b =  e.b, e.c = e.c, n = nrow(X))

rho <- opt$maximum

lm.lag <- lm((Y - rho * WY) ~ X)
summary(lm.lag)
lm.morantest(lm.lag, a.lwTN)

#SPATIAL LAG
col.lag.sm <- lagsarlm(HR8187 ~ HC8187+PE77+RDAC80+Xo+Yo, data=as.data.frame(houston), listw=a.lwTN, method="eigen") #HR8187 ~1
summary(col.lag.sm, Nagelkerke=T, correlation=TRUE)
predict.sarlm(col.lag.sm)
AIC(col.lag.sm)
deviance.sarlm(col.lag.sm)
residuals.sarlm(col.lag.sm)
coef.sarlm(col.lag.sm)
fitted.sarlm(col.lag.sm)
hetero.plot <- function(model) {
  plot(residuals(model) ~ fitted(model))
  abline(h=0, lty="dotted")
  lines(lowess(fitted(model), residuals(model)), col="red")
}
hetero.plot(col.lag.sm)
hist(residuals(col.lag.sm), col="red")
# Validación supuestos
lm.morantest(col.lag.sm, a.lwTN, alternative = "two.sided")
bptest.sarlm(col.lag.sm)
resettest(col.lag.sm)
raintest(col.lag.sm)
shapiro.test(residuals(col.lag.sm))
vif(col.lag.sm)

# Modelo Spatial Error:
col.error.sm <- errorsarlm(HR8187 ~ HC8187+PE77+RDAC80, data=as.data.frame(houston), listw=a.lwTN)
summary(col.error.sm, correlation=TRUE, Nagelkerke=T)
NK(col.error.sm,HR8187)
predict.sarlm(col.error.sm)
AIC(col.error.sm)
deviance.sarlm(col.error.sm)
residuals.sarlm(col.error.sm)
coef.sarlm(col.error.sm)
fitted.sarlm(col.error.sm)
bptest.sarlm(col.error.sm)
hetero.plot(col.error.sm)
shapiro.test(residuals(col.error.sm))
residuals.sarlm(col.error.sm)

# SARAR
col.sarar <- sacsarlm(HR8187 ~ HC8187+PE77+RDAC80, data=as.data.frame(houston), listw=a.lwTN, method="eigen")    # HR8187 ~1
summary(col.sarar, correlation=TRUE, Nagelkerke=T)
predict.sarlm(col.lag.sm)
AIC(col.sarar)
deviance.sarlm(col.sarar)
residuals.sarlm(col.sarar)
coef.sarlm(col.sarar)
fitted.sarlm(col.sarar)
bptest.sarlm(col.sarar)
hetero.plot(col.sarar)
shapiro.test(col.sarar$residuals)

# Durbin spatial 
houslagsd<-lagsarlm(HR8187 ~ HC8187+PE77+RDAC80+Xo+Yo, data=as.data.frame(houston), listw=a.lwTN, type="mixed") 
summary(houslagsd, correlation=TRUE)
predict.sarlm(houslagsd)
AIC(houslagsd)
deviance.sarlm(houslagsd)
residuals.sarlm(houslagsd)
coef.sarlm(houslagsd)
fitted.sarlm(houslagsd)
bptest.sarlm(houslagsd)
shapiro.test(houslagsd$residuals)
moran.test(residuals(houslagsd), a.lwTN, style="W")
pseudoR2.glm <- cor(exp(predict(houslagsd)),houston$HR8187)^2

# Modelo SLX
W <- as_dgRMatrix_listw(similar.listw(a.lwTN))
WY <- W%*%Yo
WX <- W%*%Xo
WHC8187<- W%*%HC8187
WPE77<- W%*%PE77
WRDAC80<- W%*%RDAC80
WPO8187<- W%*%PO8187


############################
##  Análisis de Impactos  ##
############################
# Se rezagan las variables, para ello se copia la información en houston1
houston1 <- as.data.frame(houston)
houston1$WX <- lag.listw(a.lwTN, Xo)
houston1$WY <- lag.listw(a.lwTN, Yo)
houston1$WHC8187 <- lag.listw(a.lwTN, houston1$HC8187)
houston1$WPE77 <- lag.listw(a.lwTN, houston1$PE77)
houston1$WRDAC80 <- lag.listw(a.lwTN, houston1$RDAC80)
col.lags <- nblag(hous_nbq1, 9)
a.lw1 <- nb2listw(col.lags[[1]], style="W")
a.lw2 <- nb2listw(col.lags[[2]], style="W")
a.lw3 <- nb2listw(col.lags[[3]], style="W")
a.lw4 <- nb2listw(col.lags[[4]], style="W")
a.lw5 <- nb2listw(col.lags[[5]], style="W")
a.lw6 <- nb2listw(col.lags[[6]], style="W",zero.policy =T)
a.lw7 <- nb2listw(col.lags[[7]], style="W",zero.policy =T)
houston1$W1HR8187 <- lag.listw(a.lw1,houston1$HR8187)
houston1$W2HR8187 <- lag.listw(a.lw2,houston1$HR8187)
houston1$W3HR8187<- lag.listw(a.lw3,houston1$HR8187)
houston1$W4HR8187 <- lag.listw(a.lw4,houston1$HR8187)

# Modelo SLX
col.slx2 <- lmSLX(HR8187~ HC8187 + PE77 + RDAC80+ Xo + Yo, data=houston1, listw=a.lwTN)
summary(col.slx2, Nagelkerke=T)
col.slx2 <- lmSLX(HR8187~ HC8187, data=houston1, listw=a.lwTN)
summary(col.slx2, Nagelkerke=T)
NK(col.slx2,HR8187)
predict.sarlm(col.slx2)
AIC(col.slx2)
deviance.sarlm(col.slx2)
residuals.sarlm(col.slx2)
coef.sarlm(col.slx2)
fitted.sarlm(col.slx2)
bptest(col.slx2)
hetero.plot(col.slx2)
shapiro.test(residuals(col.slx2))
hist(residuals(col.slx2), col="red")
dwtest(col.slx2)
############################
##  Análisis de Impactos  ##
############################

col.new <- houston
col.new@data <- houston1

# Cambiando la tasa de criminalidad en u
col.new@data[col.new@data$NEIG == "38","hr8187"] <- 10.0

# Los valores de las predicciones originales
orig.pred <- as.data.frame(predict(col.slx2))

W <- as(as_dgRMatrix_listw(a.lwTN), "CsparseMatrix")
trMatc <- trW(W, type="mult")
trMC <- trW(W, type="MC")

impacts(col.slx2, listw=a.lwTN)

# Nuestro modelo que trata las tasas de criminalidad CRIME (Robo de viviendas y vehículos por cada 1000 hogares): en terminos de los impactos 
# nos dice que un aumento del 100% en HOVAL (valor de la vivienda en miles de dolares) conduce a una caida del  34.9% (en promedio) en la tasa 
# de CRIME (Robo de viviendas y vehículos por cada 1000 hogares).

########################################################################
##########       Regresión Geográficamente Ponderada       #############
########################################################################

adapt <- gwr.sel(HR8187~ HC8187 + PE77 + RDAC80+ Xo + Yo+ WX + WY, data=houston1, coords=cbind(Xo,Yo)) 
gwr_fit <- gwr(HR8187~ HC8187 + PE77 + RDAC80+ Xo + Yo+ WX + WY, data=houston1, coords=cbind(Xo,Yo), adapt = 1, hatmatrix = TRUE)
gwr_fit
gwr1 <- gwr.sel(HR8187~ HC8187 + PE77 + RDAC80+ Xo + Yo+ WX + WY, data=houston1, coords=cbind(Xo,Yo), gweight=gwr.bisquare)
# Creando un objecto con el valor de Quasi-global R2
globalR2 <- (1 - (gwr_fit$results$rss/gwr_fit$gTSS))

results<-as.data.frame(gwr_fit$SDF)
head(results)
x11()
head(gwr_fit$SDF)
# Create a colour palette 
lm.palette <- colorRampPalette(c("red","orange", "white"), space = "rgb")
colours = c("dark blue", "blue", "red", "dark red")
houston$Intercept_coef <- gwr_fit$SDF@data[,2]
houston$HC8187_coef <- gwr_fit$SDF@data$HC8187
houston$PE77_coef <- gwr_fit$SDF@data$PE77
houston$RDAC80_coef <- gwr_fit$SDF@data$RDAC80
houston$localR2 <- gwr_fit$SDF@data$localR2
houston$pred <- gwr_fit$SDF@data$pred
spplot(houston, "HC8187_coef", col.regions=lm.palette(20), cex=0.6, main = "HC8187 Coef") 
spplot(houston, "localR2", key.space = "right", col.regions = lm.palette(20), cuts = 7, 
       main = "HR8187_localR2")
spplot(houston, "PE77_coef", key.space = "right", col.regions = lm.palette(20), cuts = 7, 
       main = "PE77_coef")
spplot(houston, "RDAC80_coef", key.space = "right", col.regions = lm.palette(20), cuts = 7, 
       main = "RDAC80_coef")
spplot(houston, "Intercept_coef", key.space = "right", col.regions = lm.palette(20), cuts = 7, 
       main = "Intercept_coef")

names(gwr_fit$SDF)

#Not unsurprisingly, the relationship is significant at a very high confidence level.
# But, what about the residuals ? Is there a geographical patterning to where the model over- or under-predicts? Let's look! 

# Makes the data into a "mappable" format 
map = houston

# Lists the "mappable" data 
names(map) 

#calculate t-value
t = gwr_fit$SDF$HC8187 / gwr_fit$SDF$HC8187_se  
map@data$t = t 
colours2=c("green","red","green") 

#estimated GWR t- values, red indicates a relationship that is not significant
spplot(map, "t", col.regions=lm.palette(20), cex=c(0.6,1,0.6), main = "t - Value") 

## predicted and standard error
map@data$pred <- gwr_fit$SDF@data$pred
map@data$pred.se <- gwr_fit$SDF@data$pred.se
map@data$localR2 <- gwr_fit$SDF@data$localR2
spplot(map, "pred",col.regions=lm.palette(20), main = "Predicted Value")
spplot(map, "pred.se", col.regions=lm.palette(20), main = "Standard Error")
spplot(map, "localR2", col.regions=lm.palette(20), main = "Local R2")

#######Global tests of geographical weighted regressions ###########
##Four related test statistics for comparing OLS and GWR models based on papers by Brunsdon, 
#Fotheringham and Charlton (1999) and Leung et al (2000), and a development from the GWR book (2002).


##Brunsdon, Fotheringham & Charlton (1999) ANOVA
BFC99.gwr.test(gwr_fit)
#Brunsdon, Fotheringham & Charlton (2002, pp. 91-2)
BFC02.gwr.test(gwr_fit)

# anova
anova(gwr_fit)

col.poly2 <- houston1
coordinates(col.poly2) <- c("Xo", "Yo")
xx <- gwr(HR8187~ HC8187 + PE77 + RDAC80+ Xo + Yo+ WX + WY, col.poly2, bandwidth = 9.454624, hatmatrix=TRUE)
xx

# local R2
map@data$localR2b <- xx$SDF@data$localR2
spplot(map, "localR2b", key.space="right", col.regions = lm.palette(20), cuts=7, main = "Local R2")



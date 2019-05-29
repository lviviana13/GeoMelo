library(spdep) ## LIB DATOS AREA
library(sp)
library(maptools)
library(rgdal)
library(ggplot2)
library(tmap)
library(RColorBrewer)
library(spatial)
library(RANN)
library(adespatial)
library(ade4)
library(classInt)
library(pg)

set.seed(123)
set.ZeroPolicyOption(T)
scotlip<-readOGR("D:/Loca/finalcancer3.shp")
scotlip2<-readOGR("D:/Loca/finalcancer4.shp")

scotlip<-readOGR("D:/GEOESTADISTICA/Trabajo Analisis Datos Espaciales/scot_final.shp")
coords <- coordinates(scotlip)
summary(scotlip)
plot(scotlip)
head(scotlip@data)

##### CONTIGUIDAD TORRE - REINA ####

##  ORDEN 1
scot.nbt.1 <- poly2nb(scotlip, queen=F,row.names = scotlip$DISTRICT )      # Torre orden 1
scot.nbq1 <- poly2nb(scotlip,row.names = scotlip$DISTRICT)        ####### Reina orden 1

sp.correlogram(scot.nbq1,scotlip$TCLH1,order = 10,method = "I",style = "W")
x11()
plot(scotlip)
plot(scot.nbq1, coords, add=T, col=2, lwd=1)

###### GRAFICOS BASICOS

l2 = list("SpatialPolygonsRescale", layout.north.arrow(), offset =  c(150000,1100000), scale = 70000)   # Layout
l3 = list("SpatialPolygonsRescale", layout.scale.bar, offset = c(300000,400000), scale=0.001, fill=c("transparent","black"))
l4 = list("sp.text", c(990000,1000450), "0",cex=0.8)
l5 = list("sp.text", c(991000,1000450), "1000 m",cex=0.8)
l6 <- list("sp.text",coords,  cex=0.7)
myCols <- adjustcolor(colorRampPalette(brewer.pal(n=9, 'Reds'))(100), .95)

scotlip$
tt<-spplot(scotlip["TCLH1"], scales=list(draw=TRUE), col.regions=myCols,  
       sp.layout=list(l2,l4,l5),main="N° casos de cancer\n/100.000 hab")

tt + scale_x_continuous(labels = comma)

# bubble plot equal-frequency class intervals
plotvar <- round(scotlip$TC1,3)
nclr <- 5
plotclr <- brewer.pal(nclr,"YlOrBr")
#plotclr <- plotclr[nclr:1] # reorder colors if appropriate
max.symbol.size=6
min.symbol.size=1
  class <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
symbol.size <- ((plotvar-min(plotvar))/
                  (max(plotvar)-min(plotvar))*(max.symbol.size-min.symbol.size)
                +min.symbol.size)
Data <- scotlip@data
plot(scotlip)
col.poly.cntr <- coordinates(scotlip)
points(col.poly.cntr, pch=16, col=colcode, cex=symbol.size)
points(col.poly.cntr, cex=symbol.size)
text(8.1,15.2, "Area: Equal-Frequency Class Intervals")
legend(locator(1), legend=names(attr(colcode, "table")), fill=attr(colcode, "palette"), cex=0.8, bty="n")



quantile.e(x=col.poly@data$CRIME,ic=5,digits=1,style="quantile",border=col.poly,size=0.3, xn=10.1,yn=14.5,xsb=6.1,ysb=10.7,length=1,unit="m")

coordinates(scotlip) <- c("COORD-X","COORD-Y")
bubble(scotlip,"TC1",maxsize=2,xlim = col.poly@bbox[,1],	ylim = col.poly@bbox[,2], axes=F) 
#plot(columbus,add=T)


#### K vecinos mas cercanos; con distancias #####

# IDs <- row.names(as(ken, "data.frame"))
x11()
par(mfrow=c(1,3))
scot_kn1<-knn2nb(knearneigh(coords,k=1), row.names=scotlip$DISTRICT) ## k=1 vecino + cercano
plot(scotlip)
plot(scot_kn1, coords, add=T, col=2, lwd=1)

scot_kn2<-knn2nb(knearneigh(coords,k=2), row.names=scotlip$DISTRICT) ## k=2 vecinos + cercanos
plot(scotlip)
plot(scot_kn2, coords, add=T, col=2, lwd=1)

scot_kn3<-knn2nb(knearneigh(coords,k=3), row.names=scotlip$DISTRICT) ## k=3 vecinos + cercanos
plot(scotlip)
plot(scot_kn3, coords, add=T, col=2, lwd=1)

scot_kn4<-knn2nb(knearneigh(coords,k=4), row.names=scotlip$DISTRICT) ## k=4 vecinos + cercanos
plot(scotlip)
plot(scot_kn4, coords, add=T, col=2, lwd=1)

scot_kn5<-knn2nb(knearneigh(coords,k=5), row.names=scotlip$DISTRICT) ## k=5 vecinos + cercanos
plot(scotlip)
plot(scot_kn5, coords, add=T, col=2, lwd=1)
x11()
scot_kn6<-knn2nb(knearneigh(coords,k=6), row.names=scotlip$DISTRICT) ## k=6 vecinos + cercanos
plot(scotlip)
plot(scot_kn6, coords, add=T, col=2, lwd=1)

scot_kn7<-knn2nb(knearneigh(coords,k=7), row.names=scotlip$DISTRICT) ## k=7 vecinos + cercanos
plot(scotlip)

scot_kn8<-knn2nb(knearneigh(coords,k=8), row.names=scotlip$DISTRICT) ## k=8 vecinos + cercanos
plot(scotlip)

scot_kn9<-knn2nb(knearneigh(coords,k=9), row.names=scotlip$DISTRICT) ## k=9 vecinos + cercanos
plot(scotlip)
#============================================================
# Criterios basados en gráficas
#============================================================
x11()
op=par(mfrow=c(1,1))
trinb=tri2nb(coords)
plot(scotlip,border="dimgray")
plot(trinb,coords,add=T,col="blue")
title(main="Triangulación Delaunay")

soinb=graph2nb(soi.graph(trinb,coords))
plot(scotlip,border="dimgray")
plot(soinb,coords,add=T,col="firebrick4")
title(main="Esfera de influencia")

gabrielnb=graph2nb(gabrielneigh(coords),sym=TRUE)
plot(scotlip,border="dimgray")
plot(gabrielnb,coords,add=T,col="green4")
title(main="Gráfica de Gabriel")

relativenb=graph2nb(relativeneigh(coords),sym=TRUE)
plot(scotlip,border="dimgray")
plot(relativenb,coords,add=T,col="skyblue4")
title(main="Vecinos relativos")
par(op)

################################
##### Umbral (Distancias)  #####
################################

Dist <- unlist(nbdists(knn6, coords))
summary(Dist)
max_k1 <- max(Dist)1

scot_kd1<-dnearneigh(coords, d1=0, d2=20000, row.names=scotlip$DISTRICT)                   
scot_kd2<-dnearneigh(coords, d1=0, d2=50000, row.names=scotlip$DISTRICT)                     
scot_kd3<-dnearneigh(coords, d1=0, d2=100000, row.names=scotlip$DISTRICT)                   
scot_kd4<-dnearneigh(coords, d1=0, d2=220000, row.names=scotlip$DISTRICT)

x11()
par(mar = c(0, 0, 0, 0), pty = "s")
par(mfrow=c(1,3))
plot(scotlip,main="Distancia entre\n 0 y 20.000 m")
plot(scot_kd1, coords, add=T, col=2, lwd=1)

plot(scotlip,main="Distancia entre\n 0 y 50.000 m")
plot(scot_kd2, coords, add=T, col=6, lwd=1)

plot(scotlip,main="Distancia entre\n 0 y 200.000 m")
plot(scot_kd4, coords, add=T, col=4, lwd=1)

source("D:/GEOESTADISTICA/Datos de area/Datos Área/Funciones/correlogram.d.R")
# Correlograma Moran a partir de matriz distancias                                         
corD <- correlogram.d(coords,scotlip$TCLH1,method="Moran",nbclass=3)           
corD$res                                                                                   
plot(corD$res,main="") 



par(mfrow=c(1,1))
x1<-((scotlip$TCLH1)-mean(scotlip$TCLH1))/sd(scotlip$TCLH1) 
moran.plot(X1,kn6, main="Gráfico de Dispersión de Moran", ylim=c(-2,2))
x11()
moran.plot(x1,kn6, main="Dispersograma de Moran", ylim=c(-2,2),xlab="Z (TCLH1)",
           ylab ="W Z (TCLH1)")
set.seed(127)
mc<-moran.mc(scotlip$TCLH1, kn6, nsim=1000, zero.policy=F)     # Se rechaza Ho:rho=0, es decir hay autocorrelación espacial.
mc
x11()
plot(mc, main="Simulación de Monte Carlo",xlab="TCLH1")

moran.bi(scotlip3$TCLH1,scotlip$AFF,kn6,zero.policy =  )
source("D:/GEOESTADISTICA/Datos de area/Datos Área/Funciones/moranbi1.test.R")
moranbi1.test(scotlip3$TCLH1,scotlip$AFF,kn6,randomisation = T,zero.policy =T ,alternative ="two.sided" ,
              rank = F)


# LISA Cluster Map
moran.cluster(scotlip$TCLH1, knn6, zero.policy = T, scotlip, significant=T)

## PRUEBAS I- DE MORAN
source("D:/GEOESTADISTICA/Datos de area/Datos Área/Funciones/test.W.R")
reina1<-nb2listw(scot.nbq1,style="W")
moran.test(scotlip$TCLH1,reina1,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,scot.nbq1))
set.seed(123)
kn1<-nb2listw(scot_kn1,style="W")
moran.test(scotlip$TCLH1,kn1,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,scot_kn1))

kn2<-nb2listw(scot_kn2,style="W")
moran.test(scotlip$TCLH1,kn2,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,scot_kn2))

kn3<-nb2listw(scot_kn3,style="W")
moran.test(scotlip$TCLH1,kn3,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,scot_kn3))

kn4<-nb2listw(scot_kn4,style="W")
moran.test(scotlip$TCLH1,kn4,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,scot_kn4))

kn5<-nb2listw(scot_kn5,style="W")
moran.test(scotlip$TCLH1,kn5,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,scot_kn5))

kn6<-nb2listw(scot_kn6,style="W")
moran.test(scotlip$TCLH1,kn6,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,scot_kn6))

kn7<-nb2listw(scot_kn7,style="W")
moran.test(scotlip$TCLH1,kn7,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,scot_kn7))

delaunay<-nb2listw(trinb,style="W")
moran.test(scotlip$TCLH1,delaunay,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,trinb))

esfera<-nb2listw(soinb,style="W")
moran.test(scotlip$TCLH1,esfera,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,soinb))
W.esf<-nb2listw(soinb,style="W")

gabriel<-nb2listw(gabrielnb,style="W")
moran.test(scotlip$TCLH1,gabriel,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,gabrielnb))

vecinosrel<-nb2listw(relativenb,style="W")
moran.test(scotlip$TCLH1,vecinosrel,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,relativenb))

scot_rango1<-nb2listw(scot_kd4,style="W")
moran.test(scotlip$TCLH1,scot_rango1,alternative = "two.sided")
summary(test.W(scotlip$TCLH1,relativenb))

scores.listw(kn6)
cand<-listw.candidates(coords,style ="W")
listw.explore()
listw.select(scotlip$TCLH1,cand)


moran.plot(as.vector(scale(scotlip$TCLH1)),kn6,xlab="z",ylab="Wz",main="Torre")

# Correlograma Moran k vecinos  
x11()
cor.k1 <- sp.correlogram(scot_kn1, scotlip$TCLH1 , order=2, method="I", style="W", zero.policy=T)      
 plot(cor.k1, main="Correlograma de Moran para TC1")
 cor.k1
 X11()

 cor.k1 <- sp.correlogram(scot_kn6, scotlip$LnTCLH1  , order=5, method="I", style="W", zero.policy=T)      
 plot(cor.k1, main="Correlograma de Moran para lNCOORD_Y") 
 
 cor.k6 <- sp.correlogram(scot_kn6, scotlip$TCLH1 , order=6, method="I", style="W", zero.policy=T)      
 plot(cor.k6, main="Correlograma de Moran para TCLH1\n usando criterio kn6")
 cor.k6

 
 ## MATRICES ORDENES MAYORES
x11() 
 par(mfrow=c(1,2))
lags<-nblag(scot_kn6,10) 

kn6_orden2 <-nb2listw(lags[[2]], style="W", zero.policy=T)
kn6_orden3 <-nb2listw(lags[[3]], style="W", zero.policy=T)
scot_kn6_orden4 <-nb2listw(lags[[4]], style="W", zero.policy=T)
scot_kn6_orden7 <-nb2listw(lags[[7]], style="W", zero.policy=T)
scot_kn6_orden8 <-nb2listw(lags[[8]], style="W", zero.policy=T)
scot_kn6_orden7 <-nb2listw(lags[[7]], style="W", zero.policy=T)

plot(scotlip)
plot(scot_kn6_orden2, coords, add=T, col=4, lwd=1)

moran.test(scotlip$TCLH1,scot_kn6_orden2,alternative = "two.sided")
moran.test(scotlip$TCLH1,scot_kn6_orden3,alternative = "two.sided")

plot(scotlip)
plot(scot_kn6_orden4, coords, add=T, col=4, lwd=1)

sp.correlogram(scot_kn6, scotlip$COORD_Y , order=7, method="I", style="W", zero.policy=T)       


 sp.correlogram(scot_kn7, scotlip$TCLH1 , order=7, method="I", style="W", zero.policy=T)       
  
 esf_inf <- sp.correlogram(soinb, scotlip$TCLH1 , order=7, method="I", style="W", zero.policy=T)      
 plot(esf_inf, main="Correlograma de Moran para TCLH1\n usando criterio Esfera de Influencia")
 esf_inf
    sp.correlogram
 
source("D:/GEOESTADISTICA/Datos de area/Datos Área/Funciones/correlogram.d.R") 
correlogram.d(coords,cor.k1,method = "Moran") 
# Correlograma Moran a partir de matriz distancias                                         
corD <- correlogram.d(coords,scotlip$CANCER ,method="Moran",nbclass=5)           
corD$res                                                                                   
plot(corD$res,main="")  

k1 <- knn2nb(knearneigh(coords))
all.linked <- max(unlist(nbdists(k1, coords)))
col.nb.0.all <- dnearneigh(coords, 0, all.linked, row.names=scotlip$DISTRICT)
summary(col.nb.0.all, coords)
plot(scotlip, border="blue")
plot(col.nb.0.all, coords, add=TRUE)
#==========================================================================
############### SELECCIÓN DE MATRICES DE VECINDAD POR (PCNM) ##############
#==========================================================================
#principal coordinates of neighbour matrices (PCNM, Borcard and Legendre (2002))
col_nbq1 <- poly2nb(col.poly)               # Efecto Reina

scotlq2 <- nb2listw(scot_kn2, style="W",zero.policy=T)
lag.listw(scotlq2,scotlip$Tasacan1)

library(adespatial)
w_cand<-listw.candidates(coords,style = "W")
listw.select(scotlip$TCLH1,w_cand)

list.w
listw.explore() ## IMPORTANTE

summary(test.W(scotlip$TC1,scot_kn1))
W1<-test.W(scotlip$TC1,scot_kn1)
names(W1$best)
summary(test.W(scotlip$TC1,scot_kn2))
summary(test.W(scotlip$TC1,scot_kn3))
summary(test.W(scotlip$TC1,scot_kn4))
summary(test.W(scotlip$TC1,scot_kn5))
summary(test.W(scotlip$TC1,scot_kn6))
summary(test.W(scotlip$TC1,scot_kn7))
summary(test.W(scotlip$TC1,scot_kn8))
summary(test.W(scotlip$TC1,scot_kn9))
summary(test.W(scotlip$TC1,trinb))
summary(test.W(scotlip$TC1,soinb))
summary(test.W(scotlip$TCLH1,gabrielnb))
summary(test.W(scotlip$TCLH1,relativenb))
summary(test.W(scotlip$TCLH1,scot_kd4))

summary(test.W(scotlip$TCLH1,scot_kn7))
ortho.AIC(scotlip$TCLH1,scot_kn7)
?ortho.AIC
moran.plot(x, a.lw, main="Gráfico de Dispersión de Moran", ylim=c(-2,2))      # Cambiar por "x1", para la estandarización
mp$infmat
set.seed(127)
mc<-moran.mc(x, a.lw, nsim=1000, zero.policy=F)     # Se rechaza Ho:rho=0, es decir hay autocorrelación espacial.
mc

## SILUMACION MONTE-CARLO
set.seed(1234)
scot.lknn6 <- nb2listw(scot_kn6, style="W",zero.policy=T)
bperm<-moran.mc(scotlip$TC1,scot.lknn6,nsim=1000,zero.policy = F)
bperm
var(bperm$res[1:999])
summary(bperm$res[1:999])
hist(bperm$res,freq=T,breaks = 20,xlab = "I de Moran Simulado")
abline(v=0,col="red")

X11()
p<-boxplot(TCLH1 ~ REGION, data= scotlip[-c(6,8,11), ],col=c("orange","darkolivegreen3","indianred"),las=2,
        ylab="TC1",outcol="midnightblue",par(cex.axis=0.7),par(cex.lab=1))
points(p$group, p$out, type = "p", pch=3)

p2<-boxplot(TC1 ~ AFF, data= scotlip[-c(6,8,11), ],col=c("orange","darkolivegreen3","indianred"),las=2,
           ylab="TC1",xlab="AFF",outcol="midnightblue",par(cex.axis=0.7),par(cex.lab=1))
points(p2$group, p2$out, type = "p", pch=3)

par(mfrow=c(1,1))
source("D:/GEOESTADISTICA/Datos de area/Datos Área/Funciones/correlogram.d.R")
cord <- correlogram.d(coords,scotlip$TC1,method="Moran",nbclass=20)



#############################
# Ajuste Modelo Clásico
#############################
head(scotlip2)
lm_fit <- lm(TCLH1 ~ COORD_X + COORD_Y + AFF + AREA, data=scotlip@data)
summary(lm_fit)

lm_fit2 <- lm(TCLH1 ~ COORD_X + COORD_Y + PER + AREA, data=scotlip2@data)
summary(lm_fit2)
scotlip2$lm_res <- residuals(lm_fit2)
# Mapa residuales: Opción 1
spplot(scotlip2["lm_res"], col.regions = rev(terrain.colors(20)))
lm.morantest(lm_fit2, kn6,alternative = "two.sided")
LM_res <- lm.LMtests(lm_fit2, kn6, test = "all")
print(LM_res)
summary(LM_res)

                                                 
# Por lo tanto, los residuos estan autocorrelacionados

# Mapa residuales: Opción 2
pal2 <- colorRampPalette(c("red3", "wheat1", "blue3"))
spplot(columbus,"lm_res", col.regions = pal2(20))
lm_RME <- lm(RME ~ COORD_X + COORD_Y + AFF + AREA+POP, data=scotlip@data)
summary(lm_RME)

scotlip$lm_res <- residuals(lm_fit)
# Mapa residuales: Opción 1
spplot(scotlip["lm_res"], col.regions = rev(terrain.colors(20)))

# Mapa residuales: Opción 2
pal2 <- colorRampPalette(c("red3", "wheat1", "blue3"))
spplot(columbus,"lm_res", col.regions = pal2(20))

lm.morantest(lm_fit, scot.lknn6)

LM_res <- lm.LMtests(lm_fit, scot.lknn6, test = "all")
t(sapply(LM_res, function(x) unlist(x[1:3])))

# Validación supuestos
library(lmtest)
library(VIF)
bptest(lm_fit)
resettest(lm_fit)
raintest(lm_fit)
shapiro.test(residuals(lm_fit))
vif(lm_fit)

### MODELO SAR #####
solve(); ?lagsarlm
scot.lag.sm <- lagsarlm(TCLH1 ~ COORD_X + COORD_Y + PER + AREA, data=as.data.frame(scotlip2), 
                        listw=kn6, method="eigen",tol.solve = 1e-20)


summary(scot.lag.sm, Nagelkerke=T, correlation=TRUE)

predict.sarlm(scot.lag.sm)
AIC(scot.lag.sm)

scot.lag.mixto <- lagsarlm(TC1 ~ X+Y + AFF +AREA_KM2, data=as.data.frame(scotlip), 
                        listw=scot.lknn6, method="eigen",tol.solve = 1e-20,type = "mixed")
summary(scot.lag.mixto, Nagelkerke=T, correlation=TRUE)


scot.lag.mixto2 <- lagsarlm(TC1 ~ X+Y + AFF , data=as.data.frame(scotlip), 
                           listw=scot.lknn6, method="eigen",tol.solve = 1e-20,type = "mixed")
summary(scot.lag.mixto2, Nagelkerke=T, correlation=TRUE)

scot.error <- errorsarlm(TC1 ~ X+Y + AFF + AREA_KM2, data=as.data.frame(scotlip), 
                            listw=scot.lknn6, method="eigen",tol.solve = 1e-20)
summary(scot.error, Nagelkerke=T, correlation=TRUE)


library(lmtest)
library(VIF)
library(car)
library(EnvStats)
library(MASS)
# Transformación Box-Cox
bc <- boxcox(lm_fit2, data = scotlip2, lambda = seq(-3, 3, len = 20))
L <- bc$x[which.max(bc$y)]


## ESTADISTICOS VAR TRANSF
hist(log(scotlip2$TCLH1),col="moccasin",prob=TRUE,main="Histograma de Ln(TCLH1)",xlab = "Ln(TCLH1)")
lines(density(log(scotlip2$TCLH1)),lwd=2)

qqnorm(log(scotlip2$TCLH1))
qqline(log(scotlip2$TCLH1))
 X11()
ppt<-boxplot(log(scotlip2$TCLH1),col="orange",ylab="Ln(TCLH1)",outcol="midnightblue")
points(ppt$group, ppt$out, type = "p", pch=8)


lm_fit3 <- lm( ((((TCLH1)^L)-1)/L) ~ COORD_X + COORD_Y + AFF + AREA, data=scotlip@data)
summary(lm_fit3)
shapiro.test(residuals(lm_fit3))
x11()
par(mfrow=c(1,1))
plot(lm_fit3)
plot(lm_fit)

shapiro.test(residuals(lm_fit3)) ## NORMALIDAD
bptest(lm_fit3) ## HETEROCEDASTICIDAD
resettest(lm_fit3)
raintest(lm_fit3)
shapiro.test(residuals(lm_fit3))
vif(lm_fit3)  ## MULTICOLINEALIDAD

hist(scotlip2$TCLH1)
plot(density(scotlip2$TCLH1))
hist(scotlip2$TCLH1,col="moccasin",prob=TRUE,main="Histograma de TCLH1",xlab = "TCLH1")
lines(density(scotlip2$TCLH1),lwd=2)

qqnorm(scotlip2$TCLH1)
qqline(scotlip2$TCLH1)
X11()
pp<-boxplot(TCLH1 ~ REGION, data= scotlip,col=c("orange","darkolivegreen3","indianred"),las=2,
           ylab="TCLH1",outcol="midnightblue",par(cex.axis=0.7),par(cex.lab=1))

points(p$group, p$out, type = "p", pch=3)
pp<-boxplot(scotlip2$TCLH1,col="orange",ylab="TCLH1",outcol="midnightblue")
points(pp$group, pp$out, type = "p", pch=8)

?vif
scotlip2$Tvol <- (trees$Volume^L-1)/(L*geometric.mean(trees$Volume)^(L-1))
trees$T1vol <- (trees$Volume^L-1)/L
trees$Tpvol <- trees$Volume^L
shapiro.test(residuals(lm_fit2))





predict.sarlm(scot.lag.sm)
AIC(scot.lag.sm)


col.error.sm <- errorsarlm(TC1 ~ X + Y+ AFF, data=as.data.frame(scotlip), listw=scotlq2)

col.sarar <- sacsarlm(TC1 ~ X + Y+ AFF, data=as.data.frame(scotlip), listw=scotlq2, method="eigen") 


MBCI1 <- moranbi.test(col.poly@data$CRIME,col.poly@data$INC,a.lwq1,N=999,zero.policy =T,graph=T)


moran

glm(formula = )


# I Moran bivariate

set.seed(123)
MBCH1 <- moranbi.test(scotlip3$TCLH1,scotlip$AFF,kn6,N=999,zero.policy =T,graph=T) 
# Rechazo Ho
moranbi.plot(scotlip3$TCLH1,scotlip$AFF,quiet =F,zero.policy =F,listw=kn6)

x11()
source("D:/GEOESTADISTICA/Datos de area/Datos Área/Funciones/moran.cluster.R") 
# LISA Cluster Map
moran.cluster(scotlip3$TCLH1, kn6, zero.policy = T, scotlip3, significant=T)

source("D:/GEOESTADISTICA/Datos de area/Datos Área/Funciones/getis.cluster.R") 
# Getis Cluster Map
getis.cluster(scotlip3$TCLH1, kn6, zero.policy = T, scotlip, significant=T)

###############################################
#  Mapeando outliers e influyentes locales   ##
###############################################

x1<-(scotlip3$TCLH1 -mean(scotlip3$TCLH1))/sd(scotlip3$TCLH1)
mp<-moran.plot(x1, kn6, main="Gráfico de Dispersión de Moran")      
# Cambiar por "x1", para la estandarización

infl1 <- apply(mp$is.inf, 1, any)
lhx1 <- cut(col.poly@data$CRIME, breaks=c(min(col.poly@data$CRIME), mean(col.poly@data$CRIME), 
                                          max(col.poly@data$CRIME)), labels=c("L", "H"), include.lowest=TRUE)
wx1 <- lag(a.lwq1, columbus$CRIME)
lhwx1 <- cut(wx1, breaks=c(min(wx1), mean(wx1), max(wx1)), labels=c("L", "H"), include.lowest=TRUE)
lhlh1 <- interaction(lhx1, lhwx1, infl1, drop=TRUE)
cols1 <- rep(1, length(lhlh1))
cols1[lhlh1 == "L.L.TRUE"] <- 2
cols1[lhlh1 == "L.H.TRUE"] <- 3
cols1[lhlh1 == "H.L.TRUE"] <- 4
cols1[lhlh1 == "H.H.TRUE"] <- 5
plot(columbus, col=colcode[cols1])          # gray.colors(4, 0.95, 0.55, 2.2)[cols])
legend(locator(1), legend=c("None", "LL", "LH", "HL", "HH"), fill=attr(colcode, "palette"), bty="n", cex=1, y.intersp=0.8)


localmoran(scotlip3$TCLH1,kn6)
printCoefmat(I.local[DISTRICT,],row.names(scotlip3$DISTRICT[DISTRICT]),check.names=F)
haed
nci<-moran.plot(ken$Tasa_hurto,ken.nbq_2,labels=as.character(ken$id_polig))




### ANALISIS CONFIRMATORIO DE DATOS ESPACIALES - MODELOS ESPACIALES ####
library(spdep) ## LIB DATOS AREA
library(sp)
library(rgdal)
library(car)
library(spgwr)
library(maptools)
library(ggplot2)
library(tmap)
library(RColorBrewer)
library(lmtest)
library(lattice)

set.seed(123)
set.ZeroPolicyOption(T)
scot_reg<-readOGR("D:/Loca/cancer_reg3.shp")
scotlip<-readOGR("D:/GEOESTADISTICA/Trabajo Analisis Datos Espaciales/scot_final.shp")
coords_reg <- coordinates(scotlip)

scot.kn6<-knn2nb(knearneigh(coords_reg,k=6)) , row.names=scotlip$DISTRICT)
knn6<-nb2listw(scot.kn6,style="W")

W <- as.matrix(as_dgRMatrix_listw(knn6))
sp.correlogram(scot.kn6, scot_reg$TCLH1 , order=7, method="I", style="W", zero.policy=T)
source("D:/GEOESTADISTICA/Datos de area/Datos Área/Funciones/test.W.R")
summary(test.W(scot_reg$TCLH1,scot.kn6))
kn6

## ANAMORFOSIS GAUSSIANA ##
library(RGeostats)
demo(RGeostats.start)
??anam
names(scotlip)
### MODELO LINEAL ###
lineal1 <- lm(LnTCLH1 ~ COORD_X + COORD_Y + PER +AREA+CLH_OBS+CLH_ESP, data=scotlip@data)
summary(lineal1)

lineal11 <- lm(LnTCLH1 ~  COORD_Y + PER + CLH_OBS + CLH_ESP, data=scotlip@data)
summary(lineal11)
lm.morantest(lineal11,knn6,alternative = "two.sided")
summary(lm.LMtests(lineal11, knn6, test = "all"))
shapiro.test(residuals(lineal1))
vif(lineal1)
jarque.bera.test(lineal1)

anova(scot.lag.sm,scot.error.sm,scotlagsd,mod.final)

summary.sarlm(mod.final)

sar1 <- lm(LnTCLH1 ~ W.COORD_Y + W.LnTCLH1 + PER, data=scotlip@data)
moran.test(sar1$residuals,knn6,alternative = "two.sided")
summary(sar1)
shapiro.test(residuals(sar1))
vif(sar1)

x11()
my.palette <- brewer.pal(n = 8, name = "OrRd")
scotlip$reslineal <- residuals(lineal11)
spplot(scotlip["reslineal"],col.regions = my.palette,col="transparent",cuts=7)



## MODELO SPATIAL LAG ###
scot.lag.sm <- lagsarlm(LnTCLH1 ~  W.COORD_Y + PER+CLH_OBS+CLH_ESP, data=as.data.frame(scotlip), 
                        listw=knn6, method="eigen",tol.solve = 1e-30)
summary(scot.lag.sm, Nagelkerke=T, correlation=TRUE)

scotlip$W.PER <- lag.listw(knn6, scotlip$PER)
scotlip$W.AREA <- lag.listw(knn6, scotlip$AREA)
scotlip$W.LnTCLH1 <- lag.listw(knn6, scotlip$LnTCLH1)

SAR1<-lagsarlm(LnTCLH1 ~ PER+CLH_OBS+CLH_ESP, data=as.data.frame(scotlip), 
               listw=knn6, method="eigen",tol.solve = 1e-30)

summary(SAR1, Nagelkerke=T, correlation=TRUE)
moran.test(SAR1$residuals,knn6,alternative = "two.sided")
bptest.sarlm(SAR1)
shapiro.test(residuals(SAR1))


SAR1A<-lagsarlm(LnTCLH1 ~ PER+COORD_Y, data=as.data.frame(scotlip), 
                listw=kn6_orden2, method="eigen",tol.solve = 1e-30)
summary(SAR1A, Nagelkerke=T, correlation=TRUE)
moran.test(SAR1A$residuals,kn6_orden2,alternative = "two.sided")
bptest.sarlm(SAR1A)
shapiro.test(residuals(SAR1A))

SAR1B<-lagsarlm(LnTCLH1 ~ PER+COORD_Y, data=as.data.frame(scotlip), 
                listw=kn6_orden3, method="eigen",tol.solve = 1e-30)
summary(SAR1B, Nagelkerke=T, correlation=TRUE)
moran.test(SAR1A$residuals,kn6_orden2,alternative = "two.sided")
bptest.sarlm(SAR1A)
shapiro.test(residuals(SAR1A))


SAR1C<-lagsarlm(LnTCLH1 ~ PER+CLH_OBS+CLH_ESP+W3.COORD_Y, data=as.data.frame(scotlip), 
                listw=kn6_orden2, method="eigen",tol.solve = 1e-30)
summary(SAR1C, Nagelkerke=T, correlation=TRUE)
moran.test(SAR1A$residuals,kn6_orden2,alternative = "two.sided")
bptest.sarlm(SAR1A)
shapiro.test(residuals(SAR1A))


SAR4<-lagsarlm(LnTCLH1 ~   W.PER+W3.COORD_Y+CLH_OBS, data=as.data.frame(scotlip), 
               listw=knn6, method="eigen",tol.solve = 1e-30)
summary(SAR4, Nagelkerke=T, correlation=TRUE)
moran.test(SAR4$residuals,knn6,alternative = "two.sided")
bptest.sarlm(SAR4)
shapiro.test(residuals(SAR4))
my.palette <- brewer.pal(n = 5, name = "Greens")
scotlip$resSAR4 <- residuals(SAR4)
spplot(scotlip["resSAR4"],col.regions = my.palette,col="transparent",cuts=4)

x11()
SAR5<- lagsarlm(LnTCLH1 ~W.PER+CLH_OBS+CLH_ESP, data=as.data.frame(scotlip), 
                listw=kn6_orden2, method="eigen",tol.solve = 1e-30)
summary(SAR5, Nagelkerke=T, correlation=TRUE)
moran.test(SAR5$residuals,kn6_orden2,alternative = "two.sided")
bptest.sarlm(SAR5)
shapiro.test(residuals(SAR5))
my.palette <- brewer.pal(n = 6, name = "OrRd")
scotlip$resSAR5 <- residuals(SAR5)
spplot(scotlip["resSAR5"],col.regions = my.palette,col="transparent",cuts=5)


moran.mc(scot.lag.sm$residuals,knn6,999)


my.palette <- brewer.pal(n = 12, name = "RdBu")
scotlip$resSAR1 <- residuals(SAR1)
spplot(scotlip["resSAR1"],col.regions = my.palette,col="transparent",cuts=11)


spplot(scotlip,"resSAR1",findInterval(scotlip$resSAR1, breaks,all.inside =T) , length=8, 
       col.regions=rev(brewer.pal(7,"RdBu")) )
x11()
my.palette <- brewer.pal(n = 5, name = "BuPu")
spplot(scotlip, "resSAR1", col.regions = my.palette, cuts = 4, col = "transparent")


class(SAR1$residuals)

SAR2<-lagsarlm(LnTCLH1 ~ PER+W3.COORD_Y, data=as.data.frame(scotlip), 
               listw=knn6 , method="eigen",tol.solve = 1e-30)
summary(SAR2, Nagelkerke=T, correlation=TRUE)
moran.test(SAR1$residuals,knn6,alternative = "two.sided")
bptest.sarlm(SAR1)
shapiro.test(residuals(SAR1))

moran.mc(scot.lag.sm$residuals,knn6,999)
spplot(scot.lag.sm, "residuals", at=brks, col.regions=rev(brewer.pal(grps, "RdBu")),
       col="transparent")

scot.lag.sm2 <- lagsarlm(LnTCLH1 ~ PER, data=as.data.frame(scot_reg), 
                         listw=knn6, method="eigen",tol.solve = 1e-20)
summary(scot.lag.sm2, Nagelkerke=T, correlation=TRUE)

# Validación supuestos
library(lmtest)
library(VIF)
library(car)
library(normtest)
AIC(scot.lag.sm)
bptest.sarlm(scot.lag.sm)
shapiro.test(residuals(scot.lag.sm))
remove.packages(car)

library(usdm)
vif(scot.lag.sm)

hetero.plot(mod.final)

### MODELO SPATIAL ERROR:
ERROR2 <- errorsarlm(LnTCLH1 ~  PER+CLH_OBS+ W.CLH_OBS+CLH_ESP+W.COORD_Y, data=as.data.frame(scotlip), 
                     listw=knn6, method="eigen",tol.solve = 1e-20)
summary(ERROR2, correlation=TRUE, Nagelkerke=T)
moran.test(ERROR2$residuals,knn6,alternative = "two.sided")
bptest.sarlm(ERROR2)
shapiro.test(residuals(ERROR2))
my.palette <- brewer.pal(n = 6, name = "YlGnBu")
scotlip$resERROR2 <- residuals(ERROR2)
spplot(scotlip["resERROR2"],col.regions = my.palette,col="transparent",cuts=5)

X11()
NK(ERROR2,LnTCLH1)
predict.sarlm(scot.error.sm)
AIC(ERROR2)
deviance.sarlm(scot.error.sm)
residuals.sarlm(scot.error.sm)
coef.sarlm(scot.error.sm)

fitted.sarlm(ERROR)
bptest.sarlm(ERROR)
hetero.plot(ERROR2)

# SARAR
scot.sarar <- sacsarlm(LnTCLH1 ~  PER+ CLH_OBS +CLH_ESP, data=as.data.frame(scotlip), 
                       listw=knn6, method="eigen")
summary(scot.sarar, correlation=TRUE, Nagelkerke=T)
predict.sarlm(col.lag.sm)
AIC(col.sarar)
deviance.sarlm(col.sarar)
residuals.sarlm(col.sarar)
coef.sarlm(col.sarar)
fitted.sarlm(col.sarar)
bptest.sarlm(col.sarar)
hetero.plot(col.lag.sm)
shapiro.test(col.sarar$residuals)

library(olsrr)
ols_step_all_possible(lineal1)

# Durbin spatial 
DURBIN<-lagsarlm(LnTCLH1 ~ PER +COORD_Y, data=as.data.frame(scotlip), 
                 listw=knn6, type="mixed",tol.solve = 1e-20) 
summary(DURBIN, correlation=TRUE,Nagelkerke=T)
predict.sarlm(columlagsd)
AIC(columlagsd)
deviance.sarlm(columlagsd)
residuals.sarlm(columlagsd)
coef.sarlm(columlagsd)
fitted.sarlm(columlagsd)
bptest.sarlm(columlagsd)
shapiro.test(columlagsd$residuals)

lmSLX(LnTCLH1 ~ COORD_X + COORD_Y + PER + AREA, data=as.data.frame(scot_reg), 
      listw=knn6)


glm - poisson

?glm
# Modelo glm RME
scotlip$W.RME <- lag.listw(knn6, scotlip$RME)
scot.glm<-glm(RME ~ PER + W.COORD_Y, 
              data=as.data.frame(scotlip),
              family=poisson(link="log")) 
summary(scot.glm, correlation=TRUE, Nagelkerke=T)

scot.glm$effects

residuals(scot.glm)
hetero.plot(col.glm)
?glm
lm.morantest(scot.glm,knn6,alternative = "two.sided")

moran.plot(residuals(col.glm), nb2listw(o.nb, style="W"), xlim=c(-6.5,6.5),ylim=c(-2,2))
# Calculo pseudoR2:
pseudoR2.glm <- cor(exp(predict(scot.glm)),scotlip$RME)^2

bptest(scot.glm)
# Aquí el supuesto de normalidad no aplicaría dado que se trabaja como una Poisson
shapiro.test(col.glm$residuals)
# Comportamiento residuales x11()
par(mfrow = c(2, 2))
plot(scot.glm)


### MODELO SARMA ###
scotlip.autSMA <- spautolm(LnTCLH1 ~  COORD_Y + PER+RME, data=as.data.frame(scot_reg), 
                           listw=knn6, family="SMA") 
summary(scotlip.autSMA, Nagelkerke=T)

summary(spautolm(LnTCLH1 ~ COORD_X + COORD_Y + PER + AREA, data=as.data.frame(scot_reg), 
                 listw=knn6, family="CAR") )

## MEJOR MODELO ###
## MEJOR MODELO SPATIAL LAG ###
scotlip$W.CLH_OBS <- lag.listw(knn6, scotlip$CLH_OBS)
scotlip$W.PER <- lag.listw(knn6, scotlip$PER)
scotlip$W2.PER <- lag.listw(kn6_orden2 , scotlip$PER)
scotlip$W.CLH_ESP <- lag.listw(knn6, scotlip$CLH_ESP)
scotlip$W.COORD_Y <- lag.listw(knn6, scotlip$COORD_Y)
scotlip$W2.COORD_Y <- lag.listw(kn6_orden2, scotlip$COORD_Y)
scotlip$W3.COORD_Y <- lag.listw(kn6_orden3, scotlip$COORD_Y)

mod.final <- lagsarlm(LnTCLH1 ~ PER + W.COORD_Y, data=as.data.frame(scotlip), 
                      listw=knn6, method="eigen",tol.solve = 1e-30)
summary(mod.final, Nagelkerke=T, correlation=TRUE)
bptest.sarlm(mod.final,studentize = T)
?bptest.sarlm

shapiro.test(residuals(mod.final))

jarque.bera.test(mod.final$residuals)
jb.norm.test(mod.final$residuals,nrepl = 10000)
moran.test(mod.final$residuals,knn6,alternative = "two.sided")
plot(mod.final$residuals,mod.final$fitted.values )
median(scotlip$LnTCLH1)
library(moments)
mean(scotlip$TCLH1)

scotlip@data$res<-mod.final$residuals
x11()
spplot(scotlip,"res",at=seq(min(scotlip@data$res,na.rm=TRUE),
                            max(scotlip@data$res,na.rm=TRUE), length=8), 
       col.regions=rev(brewer.pal(7,"RdBu")) )

scotlip@data$ajust<-exp(mod.final$fitted.values)
scotlip@data$ajustERROR2<-exp(ERROR2$fitted.values)
sort(scotlip@data$ajustERROR2)

x11()

lm.null <- lm(LnTCLH1 ~ PER + W.COORD_Y, data=as.data.frame(scotlip))                      
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

x11()

## GRAFICA FINAL ##
my.palette <- brewer.pal(n = 6, name = "YlOrRd")
spplot(scotlip, "TCLH1", col.regions = my.palette, cuts = 5)

spplot(scotlip, "ajustERROR2", col.regions = my.palette, cuts = 5)

spplot(scotlip,"TCLH1",at=seq(min(scotlip@data$TCLH1 ,na.rm=TRUE),
                              max(scotlip@data$TCLH1 ,na.rm=TRUE)),length=6, 
       col.regions=brewer.pal(9,"BuGn") )


spplot(scotlip,"TCLH1",col.regions=brewer.pal(5,"BuGn"),length=6 )

?spplot
sort(scotlip$LnTCLH1)

mean(scotlip$TCLH1)

sort(mod.final$fitted.values)
sort(exp(mod.final$fitted.values))

scotlip@data

ols_step

impacts(SAR1,listw =  knn6)

library(lmtest)
library(VIF)
library(car)
library(tseries)
AIC(scot.lag.sm)
bptest.sarlm(mod.final)
shapiro.test(residuals(mod.final))
vif(mod.final)
jarque.bera.test(mod.final$residuals)
jb.norm.test(mod.final$residuals,nrepl = 10000)
plot(lineal1)

########################################################################
##########       Regresión Geográficamente Ponderada       #############
########################################################################

# https://gis.stackexchange.com/questions/241127/how-to-plot-output-from-gwr-in-r
# http://geokitchen.blogspot.com.co/2012/09/r-geographically-weighted-regression.html

scotlip4 <- as.data.frame(scotlip)
adapt <- gwr.sel(LnTCLH1 ~ COORD_Y + PER +  CLH_OBS +CLH_ESP, data=scotlip4, coords=cbind(scotlip4$COORD_X,scotlip4$COORD_Y)) 
gwr_fit <- gwr(LnTCLH1 ~ COORD_Y + PER +  CLH_OBS +CLH_ESP, data=scotlip4, coords=cbind(scotlip4$COORD_X,scotlip4$COORD_Y), adapt = 1, hatmatrix = TRUE)
gwr_fit
gwr1 <- gwr.sel(LnTCLH1 ~ COORD_Y + PER +  CLH_OBS +CLH_ESP, data=scotlip4, coords=cbind(scotlip4$COORD_X,scotlip4$COORD_Y), gweight=gwr.bisquare)
# Creando un objecto con el valor de Quasi-global R2
globalR2 <- (1 - (gwr_fit$results$rss/gwr_fit$gTSS))

results<-as.data.frame(gwr_fit$SDF)
head(results)

head(gwr_fit$SDF)
# Create a colour palette 
lm.palette <- colorRampPalette(c("red","orange", "white"), space = "rgb")
colours = c("dark blue", "blue", "red", "dark red")
scotlip$intercept_coef <- gwr_fit$SDF@data[,2]
scotlip$coordy_coef <- gwr_fit$SDF@data$COORD_Y
scotlip$per_coef <- gwr_fit$SDF@data$PER
scotlip$clh_obs_coef <- gwr_fit$SDF@data$CLH_OBS
scotlip$clh_esp_coef <- gwr_fit$SDF@data$CLH_ESP
scotlip$localR2 <- gwr_fit$SDF@data$localR2
c$pred <- gwr_fit$SDF@data$pred

x11()
my.palette <- brewer.pal(n = 9, name = "OrRd")
spplot(gwr_fit$SDF) 
spplot(scotlip, "coordy_coef", col.regions=my.palette,cuts=8, cex=0.6) 
spplot(scotlip, "per_coef", col.regions=my.palette,cuts=8, cex=0.6) 
spplot(scotlip, "clh_obs_coef", col.regions=my.palette,cuts=8, cex=0.6) 
spplot(scotlip, "clh_esp_coef", col.regions=my.palette,cuts=8, cex=0.6) 
spplot(scotlip, "localR2", col.regions=my.palette,cuts=8, cex=0.6) 

spplot(gwr_fit$SDF@data[,2], "localR2", key.space = "right", col.regions = lm.palette(20), cuts = 7, 
       main = "CRIME_localR2")
spplot(col.poly, "INC_coef", key.space = "right", col.regions = lm.palette(20), cuts = 7, 
       main = "INC_coef")
spplot(scotlip3$intercept_coef, "Intercept_coef", key.space = "right", col.regions = lm.palette(20), cuts = 7, 
       main = "Intercept_coef")

names(gwr_fit$SDF)

#Not unsurprisingly, the relationship is significant at a very high confidence level.
# But, what about the residuals ? Is there a geographical patterning to where the model over- or under-predicts? Let's look! 

# Makes the data into a "mappable" format 
map = scotlip3@data

# Lists the "mappable" data 
names(map) 

#calculate t-value
t = gwr_fit$SDF$CLH_OBS / gwr_fit$SDF$CLH_OBS_se 
t2 = gwr_fit$SDF$COORD_Y / gwr_fit$SDF$COORD_Y_se
sort(t)
sort(t2)
map@data$t = t 
scotlip$CLH_OBS.t <- t

names(map)

colours2=c("green","red") 
breaks=c(min(t),-1.96,max(t))  

#estimated GWR t- values, red indicates a relationship that is not significant
spplot(scotlip, "CLH_OBS.t",cuts=breaks, col.regions=colours2, main = "t - Value") 

length(breaks)

x11()
colours2=c("green3","red3","green3") 
breaks=c(min(t),-1.96,1.96,max(t)) 
range(t)
breaks <- c(min(scotlip$CLH_OBS.t), -1.96, 1.96, max(scotlip$CLH_OBS.t))
np <- findInterval(scotlip$CLH_OBS.t, breaks,all.inside =T) 

?findInterval

colors <- c("green3","red3","green3") 
plot(scotlip, col = colors[np])

#estimated GWR t- values, red indicates a relationship that is not significant
spplot(scotlip, "CLH_OBS.t", cuts=breaks, col.regions=colours2,  main = "t - Value") 


spplot(scotlip, "CLH_OBS.t", col.regions=colours2, cex=c(0.6,1,0.6), main = "t - Value")


findInterval()

## predicted and standard error
map@data$pred <- gwr_fit$SDF@data$pred
map@data$pred.se <- gwr_fit$SDF@data$pred.se
map@data$localR2 <- gwr_fit$SDF@data$localR2
spplot(map, "pred",col.regions=lm.palette(20), main = "Predicted Value")
spplot(map, "pred.se", col.regions=lm.palette(20), main = "Standard Error")

## bubble plot 
bubble(gwr_fit$SDF, "localR2",fill = F, main = "local R2 Coef")


#######Global tests of geographical weighted regressions ###########
##Four related test statistics for comparing OLS and GWR models based on papers by Brunsdon, 
#Fotheringham and Charlton (1999) and Leung et al (2000), and a development from the GWR book (2002).


##Brunsdon, Fotheringham & Charlton (1999) ANOVA
BFC99.gwr.test(gwr_fit)
#Brunsdon, Fotheringham & Charlton (2002, pp. 91-2)
BFC02.gwr.test(gwr_fit)

# anova
anova(gwr_fit)

col.poly2 <- col.poly1
coordinates(col.poly2) <- c("X", "Y")
xx <- gwr(CRIME~ INC + HOVAL + X + Y + WX + WY, col.poly2, bandwidth = 9.454624, hatmatrix=TRUE)
xx

# local R2
map@data$localR2b <- xx$SDF@data$localR2
spplot(map, "localR2b", key.space="right", col.regions = lm.palette(20), cuts=7, main = "Local R2")
?glm

t = gwr_fit$SDF$PER / gwr_fit$SDF$PER_se 
map@data$t = t 
colours2=c("green","red","green") 

summary(mod.final, correlation=TRUE,Nagelkerke=T)

(-2.78416) + (0.048539* scotlip$PER )+ (0.45237 * W %*% scotlip$LnTCLH1)+ (0.00000282633 * scotlip$W.COORD_Y) 

mod.final$fitted.values

mod.final$parameters

############################
##  Análisis de Impactos  ##
############################

# Cambiando la tasa de criminalidad en u
scot.new2 <- scotlip
scot.new2@data <- scotlip@data
scot.new2@data[scot.new2@data$DISTRICT == "2","PER"] <- 7

# Los valores de las predicciones originales
SAR1<-lagsarlm(LnTCLH1 ~ PER+CLH_OBS+CLH_ESP, data=as.data.frame(scotlip), 
               listw=knn6, method="eigen",tol.solve = 1e-30)
sort(exp(SAR1$fitted.values))
orig.pred<-exp(SAR1$fitted.values)



SAR1_impactos<-lagsarlm(LnTCLH1 ~ PER+CLH_OBS+CLH_ESP, data=as.data.frame(scot.new2), 
                        listw=knn6, method="eigen",tol.solve = 1e-30)
sort(exp(SAR1_impactos$fitted.values))
new.pred<-exp(SAR1_impactos$fitted.values)



# Las diferencias entre las predicciones
effect.10 <- orig.pred - new.pred  
sort(effect.10)
range(effect.10)

el <- data.frame(name = scot.new2$DISTRICT , dif_pred_TCLH1 = effect.10)
scot.new2$ef10 <- el$dif_pred_TCLH1

scot.new@data
range(scot.new$ef10)
mean(scot.new$ef10)
sort(scot.new$ef10)
# Ordenando los barrios por el valor absoluto del cambio en la predicción del CRIME
el <- el[rev(order(abs(el$dif_pred_TCHL1))), ]
el[1:10, ]  #muestra los 10 primeros barrios

# Mapear estos cambios es tambien importante:

breaks <- c(min(scot.new2$ef10), -0.03, 0.03, max(scot.new2$ef10))
labels <- c("Efecto negativo (< -0.03)", "Sin efecto (-0.03 a 0.03)", 
            "Efecto positivo (> 0.03)")

# faltaba all.inside =T para evitar un poligono en blanco asociado al valor máximo  
np <- findInterval(scot.new2$ef10, sort(breaks),all.inside =T) 
colors <- c("brown2", "cornsilk2", "forestgreen")

# Dibujando el mapa
x11()
scot.new2@data
plot(scot.new2, col = colors[np])
mtext("Efectos de un cambio en el distrito 2\n sobre los valores predichos\n en el modelo SAR 1", 
      side = 3, line = 1)
legend("topleft", legend = labels, fill = colors, bty = "n")

# También podríamos mapear la magnitud de los cambios causados a raíz de la disminución de la criminalidad en el barrio 10.

pal5 <- brewer.pal(6, "Spectral")
cats5 <- classIntervals(col.new$ef10, n = 5, style = "jenks")
colors5 <- findColours(cats5, pal5)
plot(col.new, col = colors5)
legend("topleft", legend = round(cats5$brks, 2), fill = pal5, bty = "n")
mtext("Efectos de un cambio en el barrio 10, (fijando CRIME = 10)\n sobre los valores de las 
      predicciones en un modelo Durbin Espacial", side = 3, line = 1)

# While these maps and the effects for individual counties may be useful they only tell us about a 
# change in one place. The impacts() function gives us something like OLS regression coefficients for
# a spatial lag model. The logic of the impacts() function is similar to the code above, it tells you
# the direct (local), indirect (spill-over), and total effect of a unit change in each of the 
# predictor variables. The changes reported by impacts are the global average impact:


W <- as(as_dgRMatrix_listw(knn6), "CsparseMatrix")
trMatc <- trW(W, type="mult")
trMC <- trW(W, type="MC")

impacts(mod.final, listw=knn6, tr = tr(W),)

# Impacto Direct PER:
SrW.I <- solve(diag(56)-as.numeric(mod.final["rho"])*as.matrix(W))%*%diag(56)*as.numeric(coefficients(mod.final)["PER"])
sum(diag(SrW.I))/56

# Impacto Direct HOVAL:
SrW.TCHL1 <- solve(diag(56)-as.numeric(mod.final["rho"])*as.matrix(W))%*%diag(56)*as.numeric(coefficients(mod.final)["LnTCLH1"])
sum(diag(SrW.TCHL1))/56

# Impacto Direct X:
SrW.X <- solve(diag(49)-as.numeric(col.lagx["rho"])*as.matrix(W))%*%diag(49)*as.numeric(coefficients(col.lagx)["X"])
sum(diag(SrW.X))/49

# Impacto Direct Y:
SrW.Y <- solve(diag(49)-as.numeric(col.lagx["rho"])*as.matrix(W))%*%diag(49)*as.numeric(coefficients(col.lagx)["Y"])
sum(diag(SrW.Y))/49

# Impacto Direct WX:
SrW.WX <- solve(diag(49)-as.numeric(col.lagx["rho"])*as.matrix(W))%*%diag(49)*as.numeric(coefficients(col.lagx)["WX"])
sum(diag(SrW.WX))/49

# Impacto Direct Y:
SrW.WY <- solve(diag(49)-as.numeric(col.lagx["rho"])*as.matrix(W))%*%diag(49)*as.numeric(coefficients(col.lagx)["WY"])
sum(diag(SrW.WY))/49

# Impacto Indirect INC:
(sum(SrW.I) - sum(diag(SrW.I)))/49

# Impacto Indirect INC:
(sum(SrW.H) - sum(diag(SrW.H)))/49

# Impacto Total INC:
rep(1,49)%*%(SrW.I)%*%matrix(rep(1,49),ncol=1)/49              # otra forma mas simple es: sum(SrW.I)/49

# Impacto Total HOVAL:
rep(1,49)%*%(SrW.H)%*%matrix(rep(1,49),ncol=1)/49              # otra forma mas simple es: sum(SrW.H)/49

# Nuestro modelo que trata las tasas de criminalidad CRIME (Robo de viviendas y vehículos por cada 1000 hogares): en terminos de los impactos 
# nos dice que un aumento del 100% en HOVAL (valor de la vivienda en miles de dolares) conduce a una caida del  34.9% (en promedio) en la tasa 
# de CRIME (Robo de viviendas y vehículos por cada 1000 hogares).

impacts(mod.final, tr=trMatc)
impacts(mod.final, tr=trMC)

summary(lagsarlm(LnTCLH1 ~ W.COORD_Y + PER, data=as.data.frame(scotlip), 
                 listw=knn6, method="eigen",tol.solve = 1e-30))

# Procedimiento para evaluar los niveles de significancia
lobj <- (lagsarlm(LnTCLH1 ~ W.COORD_Y + PER, data=as.data.frame(scotlip), 
                  listw=knn6, method="eigen",tol.solve = 1e-30))
summary(lobj)

lobj1 <- stsls(LnTCLH1 ~ W.COORD_Y + PER, data=as.data.frame(scotlip), 
               listw=knn6)

loobj1 <- impacts(lobj1, tr=trMatc, R=200)

summary(loobj1, zstats=TRUE, short=TRUE)

lobj1r <- stsls(LnTCLH1 ~ W.COORD_Y + PER, data=as.data.frame(scotlip), 
                listw=knn6, robust=TRUE)

loobj1r <- impacts(lobj1r, tr=trMatc, R=200)

summary(loobj1r, zstats=TRUE, short=TRUE)

lobjIQ5 <- impacts(lobj, tr=trMatc, R=200, Q=5)

summary(lobjIQ5, zstats=TRUE, short=TRUE)

summary(lobjIQ5, zstats=TRUE, short=TRUE, reportQ=TRUE)
impacts(mobj, listw=a.lw)
impacts(mobj, tr=trMatc)
impacts(mobj, tr=trMC)
summary(impacts(mobj, tr=trMatc, R=200), zstats=TRUE)
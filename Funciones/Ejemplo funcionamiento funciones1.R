############################################################
## Funcionamiento moran.bi, moranbi.test y moranbi.plot:  ##
############################################################

library(spdep)
example(columbus)
col_nbq <- poly2nb(columbus)
par.lags1 <- nblag(col_nbq, 6)                  # Orden 2
e.lw2 <- nb2listw(par.lags1[[6]], style="W",zero.policy=T)
a.lw <- nb2listw(col_nbq, style="W")
CRIME <- columbus$CRIME

####################################
########     moran.bi      #########
####################################

source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/moran.bi.R")
W <- as.matrix(as_dgRMatrix_listw(a.lw))
moran.bi(columbus$CRIME,columbus$INC,a.lw,zero.policy =T)

########################################
########     moranbi.test      #########
########################################

source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/moranbi.test.R")
source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/randomize_vector.R")
set.seed(123)
MBCrime <- moranbi.test(columbus$CRIME,columbus$INC,a.lw,999,graph=T,zero.policy =T,N=1000)
moranbi.test(columbus$INC,columbus$HOVAL,a.lw,999,graph=T,zero.policy =T,N=1000)
localmoran.bi(columbus$INC,columbus$HOVAL,a.lw,zero.policy =T,alternative = "two.sided")
########################################
########     moranbi.plot      #########
########################################

source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/moranbi.plot.R")
# Editando las etiquetas de los ejes
CRIME <- as.vector(scale(columbus$CRIME))
INCOME <- as.vector(scale(columbus$INC))
moranbi.plot(CRIME,INCOME,quiet =F,zero.policy =F,listw=a.lw)
# Sin editar la etiqueta de los ejes
moranbi.plot(as.vector(scale(columbus$CRIME)),as.vector(scale(columbus$INC)),quiet =F,zero.policy =F,listw=a.lw)

#########################################
########     moran.cluster      #########
#########################################

source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/moran.cluster.R")
# LISA Cluster Map: COLUMBUS
# Si desean el mapa de significancia deben escribir T, de lo contrario F
x11()
moran.cluster(CRIME, a.lw, zero.policy = FALSE, columbus, significant=T)
moran.cluster(CRIME, e.lw2, zero.policy = FALSE, columbus, significant=T)

#########################################
########     getis.cluster      #########
#########################################

source("D:/CARLOS/Estadistica Espacial/Funciones/Lattices/getis.cluster.R")
# Getis Cluster Map: COLUMBUS
# Si desean el mapa de significancia deben escribir T, de lo contrario F
x11()
getis.cluster(CRIME, a.lw, zero.policy = FALSE, columbus, significant=T)
getis.cluster(CRIME, e.lw2, zero.policy = FALSE, columbus, significant=T)


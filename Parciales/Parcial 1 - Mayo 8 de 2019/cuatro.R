#Número de zonas
n <- 9

#Ingresar las variables

Z <- c(45,23, 86, 64,38, 7,3,9, 56)

Y <- c(5, 3, 6,2,7,4,0,1,7)


#Definir vecino Orden 2 efectoReina

A1<- c(0,1,0,1,1,1,0,0,0)
A2 <- c(1,0,1,0,1,0,0,1,0)
A3 <- c(0,1,0,1,0,1,0,1,1)
A4 <- c(1,0,1,0,1,0,0,0,1)
A5 <- c(1,1,0,1,0,0,0,0,0)
A6 <- c(1,0,1,0,0,0,0,1,0)
A7 <- c(0,0,0,0,0,0,0,1,1)
A8 <- c(0,1,1,0,0,1,1,0,0)
A9 <- c(0,0,1,1,0,0,1,0,0)

#Estandarizar la matriz 
sumavecinos <- function(vecinos)
{
  suma <- 0
  for(vecino in vecinos){
    suma <- suma + vecino
  }
  return(suma)
}

svA1 <- sumavecinos(A1)
svA2 <- sumavecinos(A2)
svA3 <- sumavecinos(A3)
svA4 <- sumavecinos(A4)
svA5 <- sumavecinos(A5)
svA6 <- sumavecinos(A6)
svA7 <- sumavecinos(A7)
svA8 <- sumavecinos(A8)
svA9 <- sumavecinos(A9)

A1/svA1 -> A1
A2/svA2 -> A2
A3/svA3 -> A3
A4/svA4 -> A4
A5/svA5 -> A5
A6/svA6 -> A6
A7/svA7 -> A7
A8/svA8 -> A8
A9/svA9 -> A9


#Concatenamos los datos
A <- c(A1,A2,A3,A4,A5,A6,A7,A8,A9)
# Matriz W
W <-matrix(A, nrow=9, ncol=9, byrow=T)

#I de Moran

moran <- function(matrizw, z, n) {
  zm <- mean(z, trim=0, na.rm=TRUE)
  suma <- 0
  sumb <- 0
  sumc <- 0
  a <- {
    for(i in seq_len(nrow(matrizw))) {
      for(j in seq_len(ncol(matrizw))) {
        suma <- suma + matrizw[i, j] * (z[i] - zm) *  (z[j] - zm)
      }
    }
  }
  b <- for (i in seq_len(nrow(matrizw))) {
    for (j in seq_len(ncol(matrizw))) {
      sumb <- sumb + matrizw[i, j]
    }
  }
  c <-for (valorz in z) {
    sumc <- sumc + (valorz-zm)^2
  }
  
  I <- (n*a)/(b*c)
  return(I)
}

lol <- moran(W,Z, n)

I <- n









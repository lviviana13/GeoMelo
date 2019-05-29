#Número de zonas
n <- 9

#Definir vecino Orden 2 efectoReina

A1 <- c(0,1,0,1,1,1,0,0,0)
A2 <- c(1,0,1,0,1,0,0,1,0)
A3 <- c(0,1,0,1,0,1,0,1,1)
A4 <- c(1,0,1,0,1,0,0,0,1)
A5 <- c(1,1,0,1,0,0,0,0,0)
A6 <- c(1,0,1,0,0,0,0,1,0)
A7 <- c(0,0,0,0,0,0,0,1,1)
A8 <- c(0,1,1,0,0,1,1,0,0)
A9 <- c(0,0,1,1,0,0,1,0,0)
#Concatenamos los datos
A <- c(A1,A2,A3,A4,A5,A6,A7,A8,A9)
#Realizar Matriz de Contigüiddad Orden 2 efecto reina
MC <-matrix(A, nrow=9, ncol=9, byrow=T)





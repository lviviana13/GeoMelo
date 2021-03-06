assign("moran.bi1",
  function(x,y,listw,zero.policy = NULL,adjust.n = TRUE, NAOK=FALSE,...){
    if (!inherits(listw, "listw")) 
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (!is.numeric(x)) 
        stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    if (!is.numeric(y)) 
        stop(paste(deparse(substitute(y)), "is not a numeric vector"))  
    if (is.null(zero.policy)) 
        zero.policy <- get("zeroPolicy", envir = .spdepOptions)
    stopifnot(is.logical(zero.policy))  
   wc <- spweights.constants(listw, zero.policy = zero.policy, adjust.n = adjust.n)
   n <- wc$n
   morans<-(n/(wc$S0))%*%((t(scale(x))%*%as.matrix(as_dgRMatrix_listw(listw))%*%scale(y))/(t(scale(x))%*%scale(x)))
 	xx <- mean(x, na.rm=NAOK)
	z <- x - xx
	zz <- sum(z^2, na.rm=NAOK)
	K <- (length(x)*sum(z^4, na.rm=NAOK))/(zz^2)
 	res <- list(I=as.vector(morans), K=K)
	res
})

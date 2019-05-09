assign("moranbi1.test",
  function(x, y, listw, randomisation=TRUE, zero.policy=NULL,
	alternative="greater", rank = FALSE, na.action=na.fail, spChk=NULL, 
	adjust.n=TRUE) {
	alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
	if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if (!is.numeric(y)) stop(paste(deparse(substitute(y)),
		"is not a numeric vector"))   
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
 	if (spChk && !chkIDs(y, listw))
		stop("Check of data and weights ID integrity failed")

#	if (any(is.na(x))) stop("NA in X")
	xname <- deparse(substitute(x))
 	yname <- deparse(substitute(y))
	wname <- deparse(substitute(listw))
	NAOK <- deparse(substitute(na.action)) == "na.pass"
	x <- na.action(x)
	xna.act <- attr(x, "na.action")
	if (!is.null(xna.act)) {
	    subset <- !(1:length(listw$neighbours) %in% xna.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
 	y <- na.action(y)
	yna.act <- attr(y, "na.action")
	if (!is.null(yna.act)) {
	    subset <- !(1:length(listw$neighbours) %in% yna.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}

	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	if (n != length(y)) stop("objects of different length")	
	wc <- spweights.constants(listw, zero.policy=zero.policy, 
		adjust.n=adjust.n)
	S02 <- wc$S0*wc$S0
	res <- moran.bi1(x, y, listw, wc$n, wc$S0, zero.policy=zero.policy,NAOK=NAOK)
	I <- res$I
	K <- res$K
	if (rank) K <- (3*(3*wc$n^2 -7))/(5*(wc$n^2 - 1))
	EI <- (-1) / wc$n1
	if(randomisation) {
		VI <- wc$n*(wc$S1*(wc$nn - 3*wc$n + 3) - wc$n*wc$S2 + 3*S02)
		tmp <- K*(wc$S1*(wc$nn - wc$n) - 2*wc$n*wc$S2 + 6*S02)
                if (tmp > VI) warning("Kurtosis overflow,\ndistribution of variable does not meet test assumptions")
		VI <- (VI - tmp) / (wc$n1*wc$n2*wc$n3*S02)
                tmp <- (VI - EI^2)
                if (tmp < 0) warning("Negative variance,\ndistribution of variable does not meet test assumptions")
		VI <- tmp
	} else {
		VI <- (wc$nn*wc$S1 - wc$n*wc$S2 + 3*S02) / (S02*(wc$nn - 1))
                tmp <- (VI - EI^2)
                if (tmp < 0) warning("Negative variance,\ndistribution of variable does not meet test assumptions")
		VI <- tmp
	}
	ZI <- (I - EI) / sqrt(VI)
	statistic <- ZI
	names(statistic) <- "Bivariate Moran Z(I) statistic"
        if (alternative == "two.sided") 
		PrI <- 2 * pnorm(abs(ZI), lower.tail=FALSE)
        else if (alternative == "greater")
            PrI <- pnorm(ZI, lower.tail=FALSE)
        else PrI <- pnorm(ZI)
	if (!is.finite(PrI) || PrI < 0 || PrI > 1) 
		warning("Out-of-range p-value: reconsider test arguments")
	vec <- c(I, EI, VI)
	names(vec) <- c("Bivariate Moran I statistic", "Expectation", "Variance")
	method <- paste("Bivariate Moran I test under", ifelse(randomisation,
	    "randomisation", "normality"))
	data.name <- paste(xname, ifelse(rank,
		"using rank correction",""), "\nweights:",
		wname, ifelse(is.null(xna.act), "", paste("\nomitted:", 
	    paste(xna.act, collapse=", "))), ifelse(is.null(yna.act), "", paste("\nomitted:", 
	    paste(yna.act, collapse=", "))),"\n")
	res <- list(statistic=statistic, p.value=PrI, estimate=vec, 
	    alternative=alternative, method=method, data.name=data.name)
	if (!is.null(xna.act)) attr(res, "na.action") <- xna.act
 	if (!is.null(yna.act)) attr(res, "na.action") <- yna.act
	class(res) <- "htest"
	res
})


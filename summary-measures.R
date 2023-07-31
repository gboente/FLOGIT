###############################################################
# COMPUTES PREDITED VALUES ON A TESTING SAMPLE
# IF THE CLASS (0-1) IS AVAILABLE RETURN ALSO DEVIANCE AND 
#                                        PEARSON RESIDUALS
###############################################################

## yy      	: response (length nn) if available
## X       	: functional data (matrix nn * lt)
## alfa 	: intercept
## beta		: Slope function evaluated on the grid ttt 
## ttt  	: grid of points where the covariates are observed (length lt)

FunF <- function (t)
{
	FunF <- 1/(1+exp(-t))
	FunF
}


Flogit.predichos <- function(yy,X,ttt, alfa, beta){
	nn <- dim(X)[1]
	lt=length(ttt)
	deltat <- 1 / (lt - 1)

	predicho <- as.vector(X %*% beta  * dt + alfa )
	proba  <- FunF(predicho) 

	########################################
	# TO AVOID NUMERICAL PROBLEMS AND NA 
	# DEFINE proba when equals to 0 or 1
	# As 10^(-300) or 1-10^(-300)
	#########################################
	for (i in 1:nn){
		if(proba[i]==0){proba[i]<-10^(-300)} 
		if(proba[i]==1){proba[1]<-1-10^(-300)}
	}

	if(is.null(yy)){
		residuo <- rep(NA,nn)
		residuo.dev <- rep(NA,nn)
		}
	else{
	#################################
	## Pearson residual
	#####################################

	residuo  <- (yy-proba)/sqrt(proba*(1-proba))
	
	#################################
	## Deviance residual
	##################################

	residuo.dev <- sign(yy-proba)*sqrt(-2*((1-yy)*log(1-proba)+ yy * log(proba)))

	#################################
	## TO AVOID NA
	#####################################
	for (ires in 1:nn){
		if(yy[ires]==1 & 1-10^(-300)<=proba[ires]){
			residuo[ires] <- 0
			residuo.dev[ires] <- 0
			}	
		if(yy[ires]==0 & proba[ires]<10^(-300)){
			residuo[ires] <- 0
			residuo.dev[ires] <- 0
			}
		if(yy[ires]==1 & proba[ires]< 10^(-300)){
			print(c("Pearson Residual infinite, redifined as 10^50, for the observation", ires))
			residuo[ires] <-  10^(10)
			print(c("Deviance Residual infinite, redifined as 10^50, for the observation", ires))
			residuo.dev[ires] <-  10^(10)
			}
		if(yy[ires]==0 & 1-10^(-300)<=proba[ires]){
			print(c("Pearson Residual -infinite, redifined as -10^10, for the observation", ires))
			residuo[ires] <- -10^(10)
			print(c("Deviance Residual -infinite, redifined as -10^10, for the observation", ires))
			residuo.dev[ires] <-  -10^(10)
			}
		}
	}
	
	return(list(residuo.Pearson=residuo, residuo.deviance=residuo.dev, predicted = predicho, prob.predict= proba))
}

#########################################
####  Functions to compute qqplots 
####         for the deviance       #####
#########################################

#######################################################
# Deviance distribution as in  Garcia Ben - yohai 2004
#######################################################

FhatD <- function(d, phat){
  	nn <- length(phat)
	for (i in 1:nn){
		if(phat[i]==0){phat[i]<-10^(-50)} 
		if(phat[i]==1){phat[1]<-1-10^(-50)}
	}
  	iA <- which( sqrt(-2*log(phat)) <=d )
  	iB <- which(-sqrt(-2*log(1-phat)) <= d )
  	efe= (sum(phat[iA]) + sum(1-phat[iB]))/nn
	return(efe)
}
 

FhatDinv <- function(ydes, phat, lower = -10, upper =10,espacio=0.001){
  	xs <- seq(lower,upper, by=espacio)
  	Fdex <- sapply(xs, FhatD, phat = phat)
  	efeinv <- min(xs[Fdex>=ydes])
	return(efeinv)
}

QQDEV.plot <- function(prob.hat, resid.dev, plt = FALSE, nombre){
	nn <- length(resid.dev)
	########################
	# COMPUTE THE QUANTILES
	########################
	ll1 <- FhatDinv(0.9999,phat = prob.hat)
	ll2 <- FhatDinv(0.0001,phat = prob.hat)
	xs <- seq(ll1,ll2,length = 1000)
	Fdex <- sapply(xs, FhatD, phat = prob.hat)
	qi<-rep(NA, length=nn)
	for( i in 1:nn){
  		qi[i] <- min(xs[Fdex>=(i-0.5)/nn])
	}

	di <- resid.dev

	l2 <- FhatDinv(0.995,phat = prob.hat)
	l1 <- FhatDinv(0.005,phat = prob.hat)
 	limite <- c(l1,l2) 
	iout1 <- which(resid.dev<l1)
	iout2 <- which(resid.dev>l2)
	iout <- c(iout1, iout2)
	if (plt) {
		pdf(nombre,bg='transparent')
		qqplot(sort(qi), sort(di), xlab = "Quantiles of estimated distribution",
       			ylab ="Deviance residuals", col="blue")
		abline(0,1,col="gray30",lwd=2)
		abline(h=l1, lty = 2,col="gray60",lwd=2)
		abline(h=l2, lty = 2,col="gray60",lwd=2)
		dev.off()
	}
	return(list(limites=limite, index.outliers=iout))
}

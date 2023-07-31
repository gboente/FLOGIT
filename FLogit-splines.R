
## -----------------------------------------------------------
## Estimation procedures
## -----------------------------------------------------------
##  BETA \in L^2(0,1)
##
## INPUT
##
## yy      : response
## X       : functional data
## type_est: Estimation procedure 
##          	'ML': Maximum likelihood
##	      	'WML':  Weighted version of ML with weights given through type_weight
##          	'M' : Bianco Yohai estimator with croux function with constant cterho
##	      	'WM': Weighted version of 'M' with weights computed through type_weight
## 
## type_weight: Weights used
##		'NONE': no weighting method
##		'HDS': hard rejection weights computed from the Donoho-Stahel estimator
##			 with cutoff qchisq(percentil,p) 
##		'BQDS': bisquare weitghs computed from the Donoho-Stahel estimator
##			  with cutoff qchisq(percentil,p) 
##	      	'FBOX': weights hard-rejection (0-1) according to  the outlyingness of 
##                  the functional boxplot as depth "MBD"
##	      	'FBOXB': weights hard-rejection according to  the outlyingness of 
##                   the functional boxplot using as depth "Both"

##    
## percentil: quantile used in the weight function to control the Mahalanobis distance
## cterho   : tuning constant for the robust estimators
##  freq    : number of splines for beta
## norder   : Order for the splines (for cubic splines norder=4)
##  ttt     : grid of points where the slope and X_i are evaluated
## -----------------------------------------------------------
## OUTPUT
##
## ordenada :intercept
## beta_coef: Coeficients of beta on the spline basis
## slope  : function beta evaluated at the grid 
## value  : minimization value 
##          mean of the deviances for 'ML', 
##          mean of rho(deviances)*weights for the others
##
## conv    : TRUE if the algorithm converge, FALSE otherwise
##


Flogitfunctional.fit <- function (yy, X, type_est,type_weight, percentil=0.975, cterho = 0.5, range_freq=NULL, norder=4, ttt) { 
      	nn <- length(yy)

      	if(is.null(range_freq)){
		freqmin= floor(max(nn^(1/5),4))
		freqmax= floor(8+2* nn^(1/5))
		range_freq <- freqmin:freqmax 
	}
      
      	lfreq=length(range_freq)
	lt=length(ttt)

  	beta.est<- matrix(rep(NA,lt*lfreq),nrow=lfreq,ncol=lt )
	alfa.est <- criterio <- objetivo <- rep(NA, lfreq)

      	for (ii in 1:lfreq){
		freq=range_freq[ii]
		est <-  Flogitfunctional(yy=yy, X=X, type_est=type_est,type_weight=type_weight, 
				percentil=percentil, cterho = cterho, freq=freq, norder=norder, ttt=ttt) 
		if (est$conv==TRUE) {
			######################################################
        		## Slope function estimate
        		######################################################
			beta.est[ii,] <- est$beta
			alfa.est[ii] <- est$ordenada
			objetivo[ii] <- est$value
			criterio[ii]   =  objetivo[ii]    +   freq  * log(nn) /   nn
       		}
		else{print(c("Estimator did not converge for basis dimension = ",freq))}
      	}
     
      	mejor <- which.min(criterio)
      	value.mejor <- objetivo[mejor]
      	freq.mejor <- (freqmin:freqmax)[mejor]
      	beta.mejor  <- beta.est[mejor,] 
      	alfa.mejor  <- alfa.est[mejor]

      	return(list(ordenada = alfa.mejor, beta=beta.mejor, value.opt = value.mejor, freq.opt=freq.mejor))
}


Flogitfunctional <- function (yy, X, type_est,type_weight, percentil=0.975, cterho = 0.5, freq, norder=4, ttt) {
	lt=length(ttt)
	deltat <- 1 / (lt - 1)
   
	 
       
   	##################################
	# COMPUTES THE SPLINE BASIS 
	##################################
	grilla.tes   <- ttt
	nodos.spl   <- seq(min(grilla.tes), max(grilla.tes), length = freq - norder + 2)
	base.beta   <- create.bspline.basis(rangeval = c(min(grilla.tes), max(grilla.tes)), 
		norder = norder, breaks = nodos.spl)
	spl.beta <- getbasismatrix(grilla.tes, base.beta)
	  
	cov_dec <- spl.beta[, 1:freq]

	###############################################################
	# Coefficients of the covariates in the B-spline basis (by row)
	###############################################################

	xx_coef <- X %*% cov_dec * deltat

	###########################################################
      	## Define the weights 
	## according to the estimation method
      	###########################################################
      	if(type_est=='ML' & type_weight!='NONE'){
		type_weight<- 'NONE'
		print("With type_est= ML Weigths equal 1, we redefine type_weight as NONE")
	}

	if(type_est=='M' & type_weight!='NONE'){
		type_weight<- 'NONE'
		print("With type_est= M Weigths equal 1, we redefine type_weight as NONE")
	}

	pesos= rep(1, length=length(yy))

	if(type_weight=='NONE'){
		pesos= rep(1, length=length(yy))
		HARDREJ = TRUE
	}

	if(type_weight=='FBOX'){
		aaa=fbplot(fit= t(X), x=ttt,  method=  "MBD", xlim=c(0,1), xlab="t",ylab="X", plot=FALSE)
		indices= aaa$out
		pesos[indices]=0  
		HARDREJ = TRUE        
	}
	 
  

	if(type_weight=='FBOXB'){
		pesos= rep(1, length=length(yy))
		aaa=fbplot(fit= t(X), x=ttt,  method=  "Both", xlim=c(0,1), xlab="t",ylab="X", plot=FALSE)
		indices= aaa$out
		pesos[indices]=0  
		HARDREJ = TRUE         
          }

 
      if(type_weight=='BQDS'){ 
    		pesos <- w.DS(xx_coef,rbst=1, percentil=percentil,cstcorte= sqrt(qchisq(percentil,1))) 
		HARDREJ = FALSE       
	}
     
	if(type_weight=='HDS'){ 
    		pesos <- wHD.DS(xx_coef,rbst=1, percentil=percentil) 
		HARDREJ = TRUE       
	}

	est <- estimar(yy, xx_coef, type_est, HARDREJ , percentil=percentil,cterho=cterho,pesos)

	if(est$conv==TRUE){
		est_slope_fun <- cov_dec %*% est$coef[-1]	
	}else{
		est$coef<- rep(NA, (freq+1))
  		est_slope_fun <- rep(NA, length(ttt))
	}
  
  return(list(ordenada = est$coef[1], beta_coef= est$coef[-1], beta=est_slope_fun, value = est$value, conv =  est$conv))
}

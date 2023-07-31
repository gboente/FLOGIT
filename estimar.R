 
## -----------------------------------------------------------
## Estimation with euclidean covariates 
## -----------------------------------------------------------
## INPUT
##
## yy      : response
## xx_coef : coefficients of X on the considered basis
## type_est: Estimation procedure 
##          'ML': Maximum likelihood
##	      'WML':  Weighted version of ML with weights given through type_weight
##          'M' : Bianco Yohai estimator with croux function with constant cterho
##	      'WM': Weighted version of 'M' with weights computed through type_weight
## 
## HARDREJ : Are the weights Hard rejection? The WML only admits Hard rejection weights

##    
## percentil: quantile used in the weight function to control the Mahalanobis
## -----------------------------------------------------------
## OUTPUT
##
## coef   : regression coeficients including intercept
## value  : minimization value 
##          mean of the deviances for 'ML', 
##          mean of rho(deviances)*weights for the others
##
## conv    : TRUE if the algorithm converge, FALSE otherwise
##-----------------------------------------------------------

 
 


estimar <- function (yy, xx_coef, type_est,HARDREJ, percentil=0.975, cterho = 0.5,pesos=1) {
  ## Initialize
 
  vv <- convergio <- NA
  n <- length(yy)
  ########################################################### 
  ## MAXIMUM LIKELIHOOD
  ########################################################### 

  if (type_est == 'ML') {
    	fit <- glm(yy ~ xx_coef, family = "binomial") ##  with intercept
    	cf  <- fit$coef
    	vv  <- fit$deviance/(2*length(yy)) ###To ensure that the result is an average
    	convergio  <- fit$converged
  }
  
  if (type_est == 'WML' & HARDREJ==TRUE) {
     	subset.obs <- (1:n)[pesos!=0] 
	fit<-glm(yy ~ xx_coef, family = "binomial",subset=subset.obs) 
	cf  <- fit$coef
	vv  <- fit$deviance/(2*length(yy)) ###To ensure that the result is an average
	convergio  <- fit$converged
  }
  
  if (type_est == 'WML' & HARDREJ==FALSE) {
     	print("ERROR: ONLY 0-1 WEIGHTS ARE ADMITTED")
  }

  ########################################################### 
  ## Robust methods
  ########################################################### 
   
  if (type_est == 'WM') {
	fit <-WBYlogistic(x0=xx_coef,y=yy,pesos=pesos,const=cterho,kmax=1000,maxhalf=10)   
    	cf  <- fit$coef
    	vv  <- fit$objective
    	convergio  <- fit$convergence
  }

  if (type_est == 'M') {
	fit <- BYlogistic(xx_coef,y=yy,const=cterho,kmax=1000,maxhalf=10)
    	cf  <- fit$coef
    	vv  <- fit$objective
    	convergio  <- fit$convergence
  }

  ########################################################### 
  ## Estimated parameters
  ########################################################### 
  
  
  return(list(coef = cf, value = vv, conv =  convergio))
}

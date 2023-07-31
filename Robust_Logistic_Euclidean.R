library(MASS)

#  Computation of the estimator of Bianco and Yohai (1996) in logistic regression
#  -------------
#  Christophe Croux, Gentiane Haesbroeck
#  Version November 2002
#  This program computes the estimator of Bianco and Yohai in 
#  logistic regression. By default, an intercept term is included 
#  and p parameters are estimated. 
#
#  For more details we refer to 
#     Croux, C., and Haesbroeck, G. (2002), ``Implementing the Bianco and Yohai estimator for 
#     Logistic Regression'' 
#
#Inputs:
#------- 
# x0= n x (p-1) matrix of explanatory variables; 
# y=n-vector of binomial response (0 or 1);
#
# initwml= logical value for selecting one of the two possible methods for computing
#          the initial value of the optimization process. 
#	   If initwml=T (default), a weighted ML estimator is computed with weights 
#	   derived from the Stahel Donoho estimator computed on the explanatory variables. 
#          If initwml=F, a classical ML fit is perfomed.
#          When the explanatory variables contain binary observations, it is recommended
#          to set initwml to F or to modify the code of the algorithm to compute the weights
#          only on the continuous variables. 
# const= tuning constant used in the computation of the estimator (default=0.5);
# kmax= maximum number of iterations before convergence (default=1000);
# maxhalf= max number of step-halving (default=10).
  
#
#Outputs:
#--------
# list with
# 1st component: T or F if convergence achieved or not
# 2nd component: value of the objective function at the minimum
# p next components: estimates for the parameters.




BYlogistic<-function(x0,y,const=0.5,kmax=1000,maxhalf=10){
	n<-dim(x0)[[1]]
	p<-dim(x0)[[2]]+1  #TO INCLUDE INTERCEPT

	#################################################################
	#Smallest value of the scale parameter before implosion
	#################################################################

	sigmamin<-0.0001



	x<-cbind(rep(1,n),x0) #TO INCLUDE INTERCEPT
	x<-as.matrix(x)
	y<-as.numeric(y)

	#################################################################
	# Computation of the initial value of the optimization process
	#################################################################
	
	wrd<- wHD.DS(x0,rbst=1,percentil=0.975)
	subset.obs <- (1:n)[wrd==1]
	gstart<-glm(y ~ x0,family=binomial,subset=subset.obs)$coef


	sigmastart<-1/sqrt(sum(gstart^2))
	xistart<-as.vector(gstart*sigmastart)
	stscores<-x %*% xistart


	#Initial value for the objective function

	oldobj<-mean(phiBY3(stscores/sigmastart,y,const))
	kstep<-1
	jhalf<-1

	while (( kstep < kmax) & (jhalf<maxhalf)){	
		optimsig<-optimize(sigmaBY3,lower=0,upper=10^3,y=y,s=stscores,c3=const)
		sigma1<-optimsig$minimum
		if (sigma1<sigmamin){
		 	print("Explosion")
			kstep<-kmax
			resultat<-list(convergence=FALSE,objective=0,coef=t(rep(NA,p)))
			return(resultat)
		}
		else {
			gamma1<-xistart/sigma1
			scores<-stscores/sigma1
			newobj<-mean(phiBY3(scores,y,const))
			oldobj<-newobj
			gradBY3<-as.vector(apply(((derphiBY3(scores,y,const))%*%t(rep(1,p)))*x,2,mean))
			h<- -gradBY3+(as.vector(gradBY3 %*% xistart) *xistart)
			finalstep<-h/sqrt(sum(h^2))
			xi1<-xistart+finalstep
			xi1<-xi1/(sum(xi1^2))
			scores1<-(x%*%xi1)/sigma1
			newobj<-mean(phiBY3(scores1,y,const))

 			hstep<-1
			jhalf<-1
			while  ((jhalf <=maxhalf) & (newobj>oldobj)){
				hstep<-hstep/2
				xi1<-xistart+finalstep*hstep
   				xi1<-xi1/sqrt(sum(xi1^2))
   				scores1<-x%*%xi1/sigma1
  				newobj<-mean(phiBY3(scores1,y,const))
   				jhalf<-jhalf+1
			}
			if ((jhalf==maxhalf+1) & (newobj>oldobj)){
				print("Convergence Achieved")
			}
			else{
				jhalf<-1
				xistart<-xi1
				oldobj<-newobj
				stscores<-x%*% xi1
				kstep<-kstep+1
			}
		}
	}

	if (kstep == kmax){
		print("No convergence")
		resultat<-list(convergence=FALSE,objective=0,coef=t(rep(NA,p))) 
		return(resultat)
	}	
	else{
		gammaest<- xistart/sigma1
		resultat<-  list(convergence=TRUE,objective=oldobj,coef=t(gammaest)) 
		return(resultat)
	}
}

 

#################################################################################
#  Computation of the weighted M-estimator of Bianco and Yohai in logistic regression
# The same weights are used to compute the initial estimator and the final weighted M-estimator.
################################################################################


WBYlogistic<-function(x0,y,pesos,const=0.5,kmax=1000,maxhalf=10){
	n<-dim(x0)[[1]]
	p<-dim(x0)[[2]]+1  #TO INCLUDE INTERCEPT

	#################################################################
	#Smallest value of the scale parameter before implosion
	#################################################################

	sigmamin<-0.0001

	x<-cbind(rep(1,n),x0) #TO INCLUDE INTERCEPT
	x<-as.matrix(x)
	y<-as.numeric(y)

	#################################################################
	# Computation of the initial value of the optimization process
	#################################################################

	subset.obs <-(1:n)[pesos!=0]
	gstart<-glm(y~x0,family=binomial,subset=subset.obs)$coef 
 
	sigmastart<-1/sqrt(sum(gstart^2))
	xistart<-as.vector(gstart*sigmastart)
	stscores<-x %*% xistart

	#################################################################
	#Initial value for the objective function
	#################################################################

	oldobj<-mean(pesos*phiBY3(stscores/sigmastart,y,const))
	kstep<-1
	jhalf<-1

	while (( kstep < kmax) & (jhalf<maxhalf)){
		optimsig<-optimize(sigmaBY3wei,lower=0,upper=10^3,y=y,s=stscores,c3=const,wei=pesos)
		sigma1<-optimsig$minimum
		if (sigma1<sigmamin){
			print("Explosion")
			kstep<-kmax
			resultat<-list(convergence=FALSE,objective=0,coef=t(rep(NA,p)))
			return(resultat)
		}
		else{
			gamma1<-xistart/sigma1
			scores<-stscores/sigma1
			newobj<-mean(pesos*phiBY3(scores,y,const))
			oldobj<-newobj
			gradBY3<-as.vector(apply(((derphiBY3(scores,y,const)*pesos)%*%t(rep(1,p)))*x,2,mean))
			h<- -gradBY3+(as.vector(gradBY3 %*% xistart)*xistart)
			finalstep<-h/sqrt(sum(h^2))
			xi1<-xistart+finalstep
			xi1<-xi1/(sum(xi1^2))
			scores1<-(x%*%xi1)/sigma1
			newobj<-mean(pesos*phiBY3(scores1,y,const))

			####stephalving

			hstep<-1
			jhalf<-1
			while((jhalf <=maxhalf) & (newobj>oldobj)){
				hstep<-hstep/2
				xi1<-xistart+finalstep*hstep
   				xi1<-xi1/sqrt(sum(xi1^2))
   				scores1<-x%*%xi1/sigma1
  				newobj<-mean(pesos*phiBY3(scores1,y,const))
   				jhalf<-jhalf+1
			}
			if ((jhalf==maxhalf+1) & (newobj>oldobj)){
				print("Convergence Achieved")
			}
			else{
				jhalf<-1
				xistart<-xi1
				oldobj<-newobj
				stscores<-x%*% xi1
				kstep<-kstep+1
			}

			}
	}
	if (kstep == kmax){
		print("No convergence")
		resultat<-list(convergence=FALSE,objective=0,coef=t(rep(NA,p))) 
		return(resultat)
	}
	else{
		gammaest<- xistart/sigma1
 		resultat<-  list(convergence=TRUE,objective=oldobj,coef=t(gammaest))
 		return(resultat)
	}
}

############################################################################
#Functions needed for the computation of estimator of Bianco and Yohai
############################################################################

phiBY3<-function(s,y,c3){
	s<-as.double(s)
	dev<-log(1+exp(-abs(s)))+abs(s)*((y-0.5)*s<0)
	res<-rhoBY3(dev,c3)+GBY3Fs(s,c3)+GBY3Fsm(s,c3)
	res
}

rhoBY3 <- function(t,c3){
	(t*exp(-sqrt(c3))*as.numeric(t <= c3))+
	(((exp(-sqrt(c3))*(2+(2*sqrt(c3))+c3))-(2*exp(-sqrt(t))*(1+sqrt(t))))*as.numeric(t >c3))
}


GBY3Fs <- function(s,c3){
	Fs<-	exp(-(log(1+exp(-abs(s)))+abs(s)*(s<0)))
	resGinf<-exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fs))))-1)
	resGinf<-(resGinf+(Fs*exp(-sqrt(-log(Fs)))))*as.numeric(s <= -log(exp(c3)-1))
	resGsup<-((Fs*exp(-sqrt(c3)))+(exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1)))*as.numeric(s > -log(exp(c3)-1))
	resG<-resGinf+resGsup
	resG
}


GBY3Fsm <- function(s,c3){
	Fsm<-exp(-(log(1+exp(-abs(s)))+abs(s)*(s>0)))
	resGinf<-exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fsm))))-1)
	resGinf<-(resGinf+(Fsm*exp(-sqrt(-log(Fsm)))))*as.numeric(s >= log(exp(c3)-1))
	resGsup<-((Fsm*exp(-sqrt(c3)))+(exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1)))*as.numeric(s < log(exp(c3)-1))
	resG<-resGinf+resGsup
	resG
}

sigmaBY3wei<-function(sigma,s,y,c3,wei){
	mean(phiBY3(s/sigma,y,c3)*wei)
}


sigmaBY3<-function(sigma,s,y,c3){
	mean(phiBY3(s/sigma,y,c3))
}

 

derphiBY3<-function(s,y,c3){
	Fs<-	exp(-(log(1+exp(-abs(s)))+abs(s)*(s<0)))
	ds<-Fs*(1-Fs)
	dev<-log(1+exp(-abs(s)))+abs(s)*((y-0.5)*s<0)
	Gprim1<-log(1+exp(-abs(s)))+abs(s)*(s<0)
	Gprim2<-log(1+exp(-abs(s)))+abs(s)*(s>0)
	derphi<- -psiBY3(dev,c3)*(y-Fs)+((psiBY3(Gprim1,c3)-psiBY3(Gprim2,c3))*ds) 
	return(derphi)
}



psiBY3 <- function(t,c3){
	psi<- (exp(-sqrt(c3))*as.numeric(t <= c3))+(exp(-sqrt(t))*as.numeric(t >c3))
	return(psi)
}
 
 
  

 
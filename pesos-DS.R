#################################
# funcion de peso
#################################
wbi<-function(u,d){
    w = 1 * (abs(u)<=d)
    v = w * (1-(u/d)^2)^2
    return(v)
}

###########################
# Weight function for x
# only when dim(x)= 1
###########################

pesosx1 <- function(x, cstcorte = 4.685){
    mediana <-median(x) 
    z<-x-mediana 
    z<-z/mad(z) 
    w<-wbi(z, cstcorte)
    weight<-w 
    weight
}


 
###########################################################
# The Donoho-stahel estimator uses the library rrcov
###########################################################

pesos.dimp.DS <- function(x,percentil){
    p        <- ncol(x)
    DSestim  <- CovSde(x)
    media    <- DSestim@center
    cova     <- DSestim@cov
    norma2.x <- sqrt(mahalanobis(x, media, cova))
    
    w1 = wbi(norma2.x,sqrt(qchisq(percentil,p)))
    return(w1)
}


pesosHD.dimp.DS <- function(x,percentil){
    p        <- ncol(x)
    DSestim  <- CovSde(x)
    media    <- DSestim@center
    cova     <- DSestim@cov
    norma2.x <- sqrt(mahalanobis(x, media, cova))
    
    w1 = 1*(norma2.x<=sqrt(qchisq(percentil,p)))
    return(w1)
}


pesosHDx1 <- function(x, percentil){
    mediana <-median(x) 
    z<-x-mediana 
    z<-z/mad(z) 
    w<-1*(z<=sqrt(qchisq(percentil,1)))
    weight<-w 
    weight
}


###########################################################
# Define the weigths
# when rbst==0, the weigths equal 1
###########################################################
##########################
# WITH BISQUARE FUNCTION
##########################
w.DS <- function(x,rbst,percentil=0.975,cstcorte=4.685){
    p<- dim(x)[[2]]
    if(rbst==0){
        w=1
    }
    if((rbst>0) & (p==1)){
        w = pesosx1(x,cstcorte=cstcorte)
    }
    if((rbst>0) & (p!=1)){
        w = pesos.dimp.DS(x,percentil=percentil)
    }
    return(w)
}

################################
#HARD REJECTION WEIGHTS
################################

wHD.DS <- function(x,rbst,percentil=0.975){
    p<- dim(x)[[2]]
    if(rbst==0){
        w=1
    }
    if((rbst>0) & (p==1)){
        w = pesosHDx1(x,percentil=percentil)
    }
    if((rbst>0) & (p!=1)){
        w = pesosHD.dimp.DS(x,percentil=percentil)
    }
    return(w)
}


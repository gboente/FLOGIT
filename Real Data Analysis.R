
## type_est: Estimation procedure 
##          'ML': Maximum likelihood
##	      'WML':  Weighted version of ML with weights given through type_weight
##          'M' : Bianco Yohai estimator with croux function with constant cterho
##	      'WM': Weighted version of 'M' with weights computed through type_weight
## 
## type_weight: Weights used
##		'NONE': no weighting method
##		'HDS': hard rejection weights computed from the Donoho-Stahel estimator
##			 with cutoff qchisq(percentil,p) 
##		'BQDS': bisquare weitghs computed from the Donoho-Stahel estimator
##			  with cutoff qchisq(percentil,p) CANNOT BE USED WITH ML
##	      'FBOX': weights hard-rejection (0-1) according to  the outlyingness of 
##                  the functional boxplot as depth "MBD"
##	      'FBOXB': weights hard-rejection according to  the outlyingness of 
##                   the functional boxplot using as depth "Both"

rm(list=ls())
 
library(dplyr)
library(fda.usc)
library(rrcov)
require('fda')          # splines
require('robust') 	# WML
require('robustbase')   # WBY
require('MASS')
library(rrcov)

source("Robust_Logistic_Euclidean.R") 
source('pesos-DS.R')		# computes the weights on the Mahalanobis distance
source("estimar.R")  		# computes the estimators based on multivariate covariates
source('FLogit-splines.R')   	# computes the functional estimator
source("summary-measures.R") 	# computes residuals and QQ plot as in Garcia Ben and Yohai (2004)

data.all <- read.csv("Data_Jan_2006_Sep_2008.csv", sep = ";")
Sys.setlocale("LC_ALL", "English")

# Holidays
feier <- as.Date(c(
  "01.01.2006",# Neujahr
  "02.01.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "03.01.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "04.01.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "05.01.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "06.01.2006",# Heilige Drei Könige
  "14.04.2006",# Karfreitag
  "17.04.2006",# Ostern
  "30.04.2006",# Tanz in den Mai
  "01.05.2006",# Maifeiertag
  "25.05.2006",# Christi Himmelfahrt
  "26.05.2006",# Brückentag
  "05.06.2006",# Pfingsten
  "15.06.2006",# Fronleichnam
  "16.06.2006",# Brückentag
  "08.08.2006",# Friedensfest (nur in Augsburg)
  "15.08.2006",# Maria Himmelfahrt
  "02.10.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien (Herbstferien-Beginn)
  "03.10.2006",# Tag der deutschen Einheit
  "31.10.2006",# Reformationstag
  "01.11.2006",# Allerheiligen
  "22.11.2006",# Buß- und Bettag
  "21.12.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "22.12.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "23.12.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "24.12.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "25.12.2006",# 1. Weihnachtstag
  "26.12.2006",# 2. Weihnachtstag
  "27.12.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "28.12.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "29.12.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "30.12.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "31.12.2006",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "01.01.2007",# Neujahr
  "02.01.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "03.01.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "04.01.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "05.01.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "06.01.2007",# Heilige Drei Könige
  "06.04.2007",# Karfreitag
  "09.04.2007",# Ostern
  "30.04.2007",# Tanz-in den Mai
  "01.05.2007",# Maifeiertag
  "17.05.2007",# Christi Himmelfahrt
  "18.05.2007",# Brückentag
  "28.05.2007",# Pfingsten
  "07.06.2007",# Fronleichnam
  "08.06.2007",# Friedensfest
  "15.06.2007",# Maria Himmelfahrt
  "03.10.2007",# Tag der deutschen Einheit
  "31.10.2007",# Reformationstag
  "01.11.2007",# Allerheiligen
  "02.11.2007",# Brückentag
  "05.11.2007",# Brückentag
  "21.11.2007",# Buß- und Bettag
  "21.12.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "22.12.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "23.12.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "24.12.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "25.12.2007",# 1. Weihnachtstag
  "26.12.2007",# 2. Weihnachtstag
  "27.12.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "28.12.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "29.12.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "30.12.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "31.12.2007",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "01.01.2008",# Neujahr
  "02.01.2008",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "02.01.2008",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "03.01.2008",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "04.01.2008",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "05.01.2008",# KEIN-FEIERTAG-ABER-RAUS: Schulferien
  "06.01.2008",# Heilige drei Könige
  "07.01.2008",# KEIN-FEIERTAG-ABER-RAUS: Brückentag
  "04.02.2008",# Rosenmontag
  "05.02.2008",# Fastnachtsdienstag
  "06.02.2008",# Aschermittwoch
  "14.02.2008",# Valentinstag
  "21.03.2008",# Karfreitag
  "23.03.2008",# Ostern
  "24.03.2008",# Ostermontag
  "30.04.2008",# Tanz-in den Mai
  "01.05.2008",# Christi Himmelfahrt
  "02.05.2008",# Brückentag
  "11.05.2008",# Pfingstsonntag
  "11.05.2008",# Muttertag
  "12.05.2008",# Pfingstmontag
  "22.05.2008",# Fronleichnam
  "23.05.2008",# Brückentag
  "21.07.2008",# Ferienbeginn in 2 Bundesländern (Meckl.-Vorp u. Schles.-Holstein)
  "08.08.2008",# Friedensfest (Augsburg)
  "04.08.2008",# Beginn Schulferien in Bayern
  "05.08.2008",# Beginn Schulferien in Bayern
  "13.08.2008",# 	Brückentag zu Mariä Himmelfahrt
  "14.08.2008",# 	Brückentag zu Mariä Himmelfahrt
  "15.08.2008" # 	Mariä Himmelfahrt
), "%d.%m.%Y")


Da.vec <- matrix(data.all$Date, nrow=24)[1,]
Da.vec <- as.Date(Da.vec, "%d.%m.%Y")

# Non-working days
DUMMYSU    <- as.numeric(weekdays(Da.vec)=="Sunday")
DUMMYSA    <- as.numeric(weekdays(Da.vec)=="Saturday")
DUMMYMO    <- as.numeric(weekdays(Da.vec)=="Monday")
DUMMYFEIER    <- as.numeric(!is.na(match(Da.vec,feier)))
DUMMYPSEUDOMO <- as.numeric(!is.na(match(Da.vec,feier+1)))
non.workdays  <- c(1:length((Da.vec)))[c(as.logical(DUMMYSU)|as.logical(DUMMYSA)|as.logical(DUMMYFEIER))]

first.day   <- min(data.all$DayIndex)
last.day    <- max(data.all$DayIndex)
days        <- c(first.day:last.day)

## Remove non-workdays from the sample
tmp  <- match(non.workdays, days)
days <- days[-tmp]
bb <- ( data.all$DayIndex %in% days )
x <- data.all[bb, ]
# x$Date <- as.Date(x$Date, "%d.%m.%Y")

## arrange data in a 638 x 24 matrix
## keep Date (so, 638 x 25)

a <- unique(x$Hour)
xx <- data.frame(matrix(NA, 638, 25))
colnames(xx) <- c('Date', paste('Hour', 1:24, sep=''))
class(xx$Date) <- "Date"
for(j in a) {
  tmp <- x[ x$Hour == j, ]
  xx[, j+1] <- tmp$Price
  xx[, 1] <- as.Date(tmp$Date, "%d.%m.%Y")
}




# Sanity check
# extract curves in a list, as for Elliptical FDA
X <- vector('list', 3)
names(X) <- c('x', 'pp', 'Load')
X$x <- split(x$Price, x$DayIndex)
X$pp <- split(x$Hour, x$DayIndex)
X$Load <- split(x$Load, x$DayIndex)

# Then, the j-th curve should be the j-th row of xx
# except for the first entry in xx that is the Date
j <- sample(nrow(xx), 1)
all.equal(as.numeric(X$x[[j]]), as.numeric(xx[j, -1]) )


# Now extract response and non-param explanatory variables
# response: daily mean Load
# non-param covariate: log( daily mean Wind )

wind.mean <- load.mean <- vector('numeric', nrow(xx))
days <- unique(x$DayIndex)
for(j in 1:length(days)) {
  wind.mean[j] <- mean( x$Windinfeed[ x$DayIndex == days[j] ] )
  load.mean[j] <- mean( x$ResLoad[ x$DayIndex == days[j] ] )
  # load.max[j] <- max( x$Load[ x$DayIndex == days[j] ] )
}

y <- load.mean
u <- log(wind.mean)
tt <- 1:24
x <- as.matrix( xx[, -1] )
colnames(x) <- row.names(x) <- NULL

rng <- 3:12
par(mfrow=c(1,1))
matplot(t(x), type = "l")


X <- x
Demanda <- y
ttt<-tt
yy <- 1*(Demanda>median(Demanda))

#######################################
# PLOT DATA
########################################

nombre="datos-electricidad.pdf"

pdf(nombre, bg='transparent')
matplot(ttt,t(X), type = "l", col = y + 2, lty = 1,lwd=2, xlab="Time",ylab="Price")
dev.off()


nombre="datos-electricidad-grupos-rayas.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X[yy==1,]), type = "l", col = "maroon", 
          lty = 1,lwd=2,ylim=c(0,2500), xlab="Time",ylab="Price")
matplot(ttt,t(X[yy==0,]), type = "l", col = "chartreuse4", 
          lty = 2,lwd=1,ylim=c(0,2500), xlab="Time",ylab="Price",add=T)
dev.off()



nombre="datos-electricidad-grupos-lineas.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X[yy==1,]), type = "l", col = "maroon", 
          lty = 1,lwd=2,ylim=c(0,2500), xlab="Time",ylab="Price")
matplot(ttt,t(X[yy==0,]), type = "l", col = "chartreuse4", 
          lty = 1,lwd=1,ylim=c(0,2500), xlab="Time",ylab="Price",add=T)
dev.off()

aaa=fbplot(fit= t(X), x=ttt,  method=  "Both", 
	xlim=c(min(ttt),max(ttt)), ylim=c(0,2500),xlab="t",ylab="X")
 
indices= aaa$out
Fechas=xx$Date[indices]


nombre="datos-electricidad-grupos-Fechas.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X[yy==1,]), type = "l", col = "maroon", 
          lty = 1,lwd=2,ylim=c(0,2500), xlab="Time",ylab="Price")
matplot(ttt,t(X[yy==0,]), type = "l", col = "chartreuse4", 
          lty = 1,lwd=1,ylim=c(0,2500), xlab="Time",ylab="Price",add=T)

text(12, X[indices[6],12]+50, Fechas[6], cex=0.8, col="chartreuse4")

text(12, X[indices[8],12]+50, Fechas[8], cex=0.8, col="chartreuse4")
 
text(19, X[indices[10],19]+50, Fechas[10], cex=0.8, col = "maroon")
text(18, X[indices[18],18]+50, Fechas[18], cex=0.8, col = "maroon")
dev.off()

X.aaa <- X[-indices,]
yy.aaa <- yy[-indices]


nombre="datos-electricidad-grupos-rayas-sinoutfbplot.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X.aaa[yy.aaa==1,]), type = "l", col = "maroon", 
          lty = 1,lwd=2,ylim=c(0,350), xlab="Time",ylab="Price")
matplot(ttt,t(X.aaa[yy.aaa==0,]), type = "l", col = "chartreuse4", 
          lty = 2,lwd=2,ylim=c(0,350), xlab="Time",ylab="Price",add=T)
dev.off()

###########################################
#### COMPUTE ESTIMATORS
###########################################
nsamp <- dim(X)[1]
freqmin= floor(max(nsamp^(1/5),4))
freqmax= floor(8+2* nsamp^(1/5))
rango.freq=freqmin:freqmax
lfreq=length(rango.freq)
 

lt=length(ttt)
deltat <- 1 / (lt - 1)


grilla_spl <- ttt 
dt=(grilla_spl[2]-grilla_spl[1])/(max(grilla_spl)-min(grilla_spl))


norder <- 4
percentil <- 0.975
cterho = 0.5


#######################
# ML
#######################

type_est = 'ML'	

type_weight='NONE'

 
ajuste.ML <- Flogitfunctional.fit(yy=yy, X=X, type_est=type_est,type_weight=type_weight, range_freq=rango.freq, norder=4, ttt=ttt)  

alfa.ML <- ajuste.ML$ordenada
beta.ML <- ajuste.ML$beta
freq.ML <-ajuste.ML$freq.opt
mejor.ML <- ajuste.ML$value.opt  
 
result.ML <- Flogit.predichos(yy=yy,X=X,ttt, alfa=alfa.ML, beta=beta.ML) 
predicho.ML <- result.ML$predicted
prob.ML <- result.ML$prob.predict

residuo.ML <- result.ML$residuo.Pearson
residuo.dev.ML <- result.ML$residuo.deviance

#####################################
#M ESTIMATOR
######################################

type_est = 'M'	

type_weight='NONE'
 
ajuste.M  <- Flogitfunctional.fit(yy=yy, X=X, type_est=type_est,type_weight=type_weight, range_freq=rango.freq, norder=4, ttt=ttt)  

alfa.M  <- ajuste.M$ordenada
beta.M  <- ajuste.M$beta
freq.M <-ajuste.M$freq.opt
mejor.M  <- ajuste.M$value.opt   


result.M <- Flogit.predichos(yy=yy,X=X,ttt, alfa=alfa.M, beta=beta.M) 
predicho.M <- result.M$predicted
prob.M <- result.M$prob.predict

residuo.M <- result.M$residuo.Pearson
residuo.dev.M <- result.M$residuo.deviance


#################################################################
# WML with hard rejection weights based on Mahalanobis distance
#################################################################

type_est = 'WML'	
type_weight='HDS' 

ajuste.WML.HDS <- Flogitfunctional.fit(yy=yy, X=X, type_est=type_est,type_weight=type_weight, range_freq=rango.freq, norder=4, ttt=ttt)  

alfa.WML.HDS <- ajuste.WML.HDS$ordenada
beta.WML.HDS <- ajuste.WML.HDS$beta
freq.WML.HDS <-ajuste.WML.HDS$freq.opt
mejor.WML.HDS <- ajuste.WML.HDS$value.opt  
 
 
result.WML.HDS <- Flogit.predichos(yy=yy,X=X,ttt, alfa=alfa.WML.HDS, beta=beta.WML.HDS) 
predicho.WML.HDS <- result.WML.HDS$predicted
prob.WML.HDS <- result.WML.HDS$prob.predict

residuo.WML.HDS <- result.WML.HDS$residuo.Pearson
residuo.dev.WML.HDS <- result.WML.HDS$residuo.deviance
 

#################################################################
# WM  with hard rejection weights based on Mahalanobis distance
#################################################################

type_est = 'WM'	
type_weight='HDS' 

ajuste.WM.HDS <- Flogitfunctional.fit(yy=yy, X=X, type_est=type_est,type_weight=type_weight, range_freq=rango.freq, norder=4, ttt=ttt)  

alfa.WM.HDS <- ajuste.WM.HDS$ordenada
beta.WM.HDS <- ajuste.WM.HDS$beta
freq.WM.HDS <-ajuste.WM.HDS$freq.opt
mejor.WM.HDS <- ajuste.WM.HDS$value.opt 


result.WM.HDS <- Flogit.predichos(yy=yy,X=X,ttt, alfa=alfa.WM.HDS, beta=beta.WM.HDS) 
predicho.WM.HDS <- result.WM.HDS$predicted
prob.WM.HDS <- result.WM.HDS$prob.predict

residuo.WM.HDS <- result.WM.HDS$residuo.Pearson
residuo.dev.WM.HDS <- result.WM.HDS$residuo.deviance

  
#########################################################
# RESULTS INTERCEPT and DIMENSION B-SPLINES BASIS
#########################################################

c(  alfa.WM.HDS, alfa.WML.HDS,  alfa.M, alfa.ML) 
 

c(  freq.WM.HDS, freq.WML.HDS,  freq.M, freq.ML) 
 

################################
# PLOTS OF BETA
################################

 
nombre="beta-ML-M-WMHDS.pdf"
pdf(nombre, bg='transparent')
par(mar=c(4,5,3,3))

maximo=max(c(beta.WM.HDS,beta.M,beta.ML))
minimo=min(c(beta.WM.HDS,beta.M,beta.ML))


plot(ttt,beta.WM.HDS, type="l", col="blue",lwd=2, xlab= "Time (in hours)", 
      ylim=c(minimo,maximo), ylab=expression(hat(beta)))

lines(ttt,beta.M, type="l", col="cyan",lwd=2)
lines(ttt,beta.ML, type="l", col="tomato",lwd=2)
 
dev.off()
 


##################################
# BOXPLOTS PREDICTED PROBABILITIES
###################################
 
nombre="prob-predicha-ML.pdf"
pdf(nombre,bg='transparent')

boxplot(prob.ML[yy==0], prob.ML[yy==1], col=c("gold","violet"), names=c("Low Demand", "High Demand"))
   abline(h=0.5, col="gray30",lty=2,lwd=2)

dev.off()


nombre="prob-predicha-M.pdf"
pdf(nombre,bg='transparent')

boxplot(prob.M[yy==0], prob.M[yy==1], col=c("gold","violet"), names=c("Low Demand", "High Demand"))
   abline(h=0.5, col="gray30",lty=2,lwd=2)

dev.off()


nombre="prob-predicha-WML-HDS.pdf"
pdf(nombre,bg='transparent')

boxplot(prob.WML.HDS[yy==0], prob.WML.HDS[yy==1], col=c("gold","violet"), names=c("Low Demand", "High Demand"))
   abline(h=0.5, col="gray30",lty=2,lwd=2)

dev.off()

 

nombre="prob-predicha-WM-HDS.pdf"
pdf(nombre,bg='transparent')

boxplot(prob.WM.HDS[yy==0], prob.WM.HDS[yy==1], col=c("gold","violet"), names=c("Low Demand", "High Demand"))
   abline(h=0.5, col="gray30",lty=2,lwd=2)

dev.off()

#################################
## Boxplots Pearson residuals
#####################################
 

names(residuo.ML) <- names(residuo.M)<-   names(residuo.WM.HDS)<- names(residuo.WML.HDS)<-1:length(yy)

res <- abs(residuo.ML)
out.pearson<-as.numeric(names(adjbox(res)$out))
inliers <- res[out.pearson] < median(res,na.rm =TRUE )
out.res.pearson.ML <- out.pearson[!inliers] 
 

res <- abs(residuo.WML.HDS )
out.pearson<-as.numeric(names(adjbox(res)$out))
inliers <- res[out.pearson] < median(res,na.rm =TRUE)
out.res.pearson.WML.HDS <- out.pearson[!inliers] 
 


 
res <- abs(residuo.M )
out.pearson<-as.numeric(names(adjbox(res)$out))
inliers <- res[out.pearson] < median(res,na.rm =TRUE)
out.res.pearson.M <- out.pearson[!inliers] 
 

 

res <- abs(residuo.WM.HDS )
out.pearson<-as.numeric(names(adjbox(res)$out))
inliers <- res[out.pearson] < median(res,na.rm =TRUE)
out.res.pearson.WM.HDS <- out.pearson[!inliers] 
 

######################################################
## QQ PLOT OF Deviance RESIDUALS COMPUTED WITH WM-HDS
######################################################
 

names(residuo.dev.ML) <- names(residuo.dev.M)<- 
names(residuo.dev.WM.HDS)<- names(residuo.dev.WML.HDS)<- 
	names(prob.WM.HDS)<- names(prob.WML.HDS)<-1:length(yy)
 

nombre="QQPLOT-DEV-WM-HDS.pdf" 

grafico <- QQDEV.plot(prob.hat=prob.WM.HDS, resid.dev=residuo.dev.WM.HDS, plt = TRUE, nombre)   

bounds.WM.HDS<- grafico$limites
out.WM.HDS<- as.numeric(names(grafico$index.outliers))   
 
out.WM.HDS

##################################################
# CLASSICAL WITHOUT OUTLIERS Detected by WM-HDS
##################################################

outliers <- out.WM.HDS
X.sout<- X[-outliers,]
yy.sout <- yy[-outliers]

nsamp.sout=length(yy.sout)
freqmin.sout= floor(max(nsamp.sout^(1/5),4))
freqmax.sout= floor(8+2* nsamp.sout^(1/5))
rango.freq.sout <- freqmin.sout:freqmax.sout

criterio='bic1'
norder= 4
percentil=0.975
cterho = 0.5


#######################
# ML
#######################

type_est = 'ML'	
type_weight='NONE'

 
ajuste.ML.sout <- Flogitfunctional.fit(yy=yy.sout, X=X.sout, type_est=type_est,type_weight=type_weight, range_freq=rango.freq.sout, norder=4, ttt=ttt)  

alfa.ML.sout <- ajuste.ML.sout$ordenada
beta.ML.sout <- ajuste.ML.sout$beta
freq.ML.sout <-ajuste.ML.sout$freq.opt
mejor.ML.sout <- ajuste.ML.sout$value.opt  
 
result.ML.sout <- Flogit.predichos(yy=yy,X=X,ttt, alfa=alfa.ML.sout, beta=beta.ML.sout) 
predicho.ML.sout <- result.ML.sout$predicted
prob.ML.sout <- result.ML.sout$prob.predict
 

################################################
# RESULT
##################################################

c(  alfa.WM.HDS, alfa.WML.HDS,   alfa.M, alfa.ML,alfa.ML.sout) 
 
c( freq.WM.HDS, freq.WML.HDS,  freq.M, freq.ML,freq.ML.sout) 
 

################################
# PLOTS OF BETA
################################

 

nombre="beta-ML-WMHDS-MLSOUT-OUTHDS.pdf"
pdf(nombre, bg='transparent')
par(mar=c(4,5,3,3))

maximo=max(c(beta.WM.HDS,beta.ML.sout,beta.ML))
minimo=min(c(beta.WM.HDS,beta.ML.sout,beta.ML))


plot(ttt,beta.WM.HDS, type="l", col="blue",lwd=2, xlab= "Time (in hours)", 
      ylim=c(minimo,maximo), ylab=expression(hat(beta)))
lines(ttt,beta.ML, type="l", col="tomato",lwd=2)
lines(ttt,beta.ML.sout, type="l", lty=2, col="maroon",lwd=2)
 
dev.off()



nombre="beta-ML-MLSOUT-black.pdf"
pdf(nombre, bg='transparent')
par(mar=c(4,5,3,3))



maximo=max(c(beta.WM.HDS,beta.ML.sout,beta.ML))
minimo=min(c(beta.WM.HDS,beta.ML.sout,beta.ML))


plot(ttt,beta.WM.HDS, type="n", col="black",lwd=3,  xlab= "Time (in hours)", 
      ylim=c(minimo,maximo), ylab=expression(hat(beta)), cex.lab=1.3,cex=1.2)

lines(ttt,beta.ML, type="l", col="red",lwd=3)

lines(ttt,beta.ML.sout, type="l", lty=2, col="maroon",lwd=3)
 
dev.off()


nombre="beta-ML-M-WMLHDS-black2.pdf"
pdf(nombre, bg='transparent')
par(mar=c(4,5,3,3))

maximo=max(c(beta.WM.HDS,beta.ML.sout,beta.ML))
minimo=min(c(beta.WM.HDS,beta.ML.sout,beta.ML))



plot(ttt,beta.WM.HDS, type="l", col="black",lwd=3, xlab= "Time (in hours)", 
      ylim=c(minimo,maximo), ylab=expression(hat(beta)), cex.lab=1.3,cex=1.2)

lines(ttt,beta.M, type="l", col="gray70",lwd=3)
lines(ttt,beta.ML, type="l", col="red",lwd=3)

lines(ttt,beta.WML.HDS, type="l", col="red",lty=2,lwd=3)


 


dev.off()

nombre="beta-ML-WMHDS-MLSOUT-OUTHDS-black3.pdf"
pdf(nombre, bg='transparent')
par(mar=c(4,5,3,3))



maximo=max(c(beta.WM.HDS,beta.ML.sout,beta.ML))
minimo=min(c(beta.WM.HDS,beta.ML.sout,beta.ML))


plot(ttt,beta.WM.HDS, type="l", col="black",lwd=3,  xlab= "Time (in hours)", 
      ylim=c(minimo,maximo), ylab=expression(hat(beta)), cex.lab=1.3,cex=1.2)

lines(ttt,beta.ML, type="l", col="red",lwd=3)

lines(ttt,beta.ML.sout, type="l", lty=3, col="magenta",lwd=3)
 
dev.off()

 

################################################
# PREDICTED PROBABILTIES OVER THE WHOLE SAMPLE
################################################

nombre="prob-predicha-ML-SOUT-OUTHDS.pdf"
pdf(nombre, bg='transparent')

boxplot(prob.ML.sout[yy==0], prob.ML.sout[yy==1], col=c("gold","violet"), names=c("No metabolic syndrome", "Metabolic syndrome"))
   abline(h=0.5, col="gray30",lty=2,lwd=2)

dev.off()



########################################################### 
#  PREDICTED PROBABILTIES OVER THE SAMPLE WITHOUT OUTLIERS
###########################################################

prob.ML.sout.X.sout <- prob.ML.sout[-outliers]
 
 
prob.ML.X.sout <- prob.ML[-outliers]
 
 
prob.M.X.sout <- prob.M[-outliers]
 
prob.WML.HDS.X.sout <- prob.WML.HDS[-outliers]
 

prob.WM.HDS.X.sout <- prob.WM.HDS[-outliers]
 
 

nombre="prob-predicha-ML-SOUT-XSOUT-OUTHDS.pdf"
pdf(nombre, bg='transparent')

boxplot(prob.ML.sout.X.sout[yy.sout ==0], prob.ML.sout.X.sout[yy.sout ==1], col=c("gold","violet"), names=c("Low Demand", "High Demand"))
   abline(h=0.5, col="gray30",lty=2,lwd=2)

dev.off()

nombre="prob-predicha-ML-XSOUT-OUTHDS.pdf"
pdf(nombre, bg='transparent')

boxplot(prob.ML.X.sout[yy.sout==0], prob.ML.X.sout[yy.sout==1], col=c("gold","violet"), names=c("Low Demand", "High Demand"))
   abline(h=0.5, col="gray30",lty=2,lwd=2)

dev.off()

nombre="prob-predicha-WML-HDS-XSOUT-OUTHDS.pdf"
pdf(nombre, bg='transparent')

boxplot(prob.WML.HDS.X.sout[yy.sout==0], prob.WML.HDS.X.sout[yy.sout==1], col=c("gold","violet"), names=c("Low Demand", "High Demand"))
   abline(h=0.5, col="gray30",lty=2,lwd=2)

dev.off()

 

nombre="prob-predicha-M-XSOUT-OUTHDS.pdf"
pdf(nombre, bg='transparent')

boxplot(prob.M.X.sout[yy.sout==0], prob.M.X.sout[yy.sout==1], col=c("gold","violet"),names=c("Low Demand", "High Demand"))
   abline(h=0.5, col="gray30",lty=2,lwd=2)

dev.off()

 

nombre="prob-predicha-WM-HDS-XSOUT-OUTHDS.pdf"
pdf(nombre, bg='transparent')

boxplot(prob.WM.HDS.X.sout[yy.sout==0], prob.WM.HDS.X.sout[yy.sout==1], col=c("gold","violet"), names=c("Low Demand", "High Demand"))
   abline(h=0.5, col="gray30",lty=2,lwd=2)

dev.off()

 
###############################################
# TABLE WITH OUTLIERS AND RESIDUALS of WM-HDS
###############################################


tabla=cbind(outliers, cbind.data.frame(xx$Date[outliers]), residuo.dev.WM.HDS[outliers], 
       yy[outliers], prob.WM.HDS[outliers])

tabla

tablita=data.frame(outliers, cbind.data.frame(xx$Date[outliers]), residuo.dev.WM.HDS[outliers], 
       yy[outliers], prob.WM.HDS[outliers])

dates <- sapply(x,FUN = function(x){class(x) == "Date"})
tablita[,dates] <- as.character(tablita[,dates])

############################################
# TABLE
########################################


library(xtable) 

xtable2 <- function(x, ...) {
   # get the names of variables that are dates by inheritance
   datevars <- colnames(x)[vapply(x, function(y) {
       inherits(y, c("Date", "POSIXt", "POSIXct"))
   }, logical(1))]
   for (i in datevars){
        x[ , i] <- as.character(x[, i])
   }
   xtable::xtable(x, ...)
}

colnames(tabla[,-1])   <- c('Date','Residual', 'y', expression(hat(p))) 


 
rownames(tabla[,-1]) <-outliers

print(xtable2(tabla[,-1],  digits=3 ,align='c|cccc|', 
caption='Outliers identified by WM with weights based on the Mahalanobis distance and hard rejection weight function'), type='latex')

 
#####################################################
# SOME PLOTS  
#####################################################


nombre="datos-electricidad-outliers-OUTHDS.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X), type = "l", col = "gray60", 
          lty = 1,lwd=2,ylim=c(0,2500), xlab="Time",ylab="Price")
matplot(ttt,t(X[outliers,]), type = "l", col = "blue4", 
          lty = 2,lwd=3,ylim=c(0,2500), xlab="Time",ylab="Price",add=T)
dev.off()


nombre="datos-electricidad-grupo1-outliers-OUTHDS.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X[yy==1,]), type = "l", col = "maroon", 
          lty = 1,lwd=2,ylim=c(0,2500), xlab="Time",ylab="Price")

for (i in outliers){
if(yy[i]==1){ lines(ttt, X[i,] , type = "l", col = "blue4", 
          lty = 2,lwd=3)
	}
} 

dev.off()


nombre="datos-electricidad-grupo0-outliers-OUTHDS.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X[yy==0,]), type = "l", col = "chartreuse4", 
          lty = 1,lwd=2,ylim=c(0,2500), xlab="Time",ylab="Price")

for (i in outliers){
if(yy[i]==0){ lines(ttt, X[i,] , type = "l", col = "blue4", 
          lty = 2,lwd=3)
	}
} 

dev.off()



nombre="datos-electricidad-grupo1-outliers-OUTHDS-Fechas-black.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X[yy==1,]), type = "l", col = "gray40", 
          lty = 1,lwd=2,ylim=c(0,2500), xlab="Time",ylab="Price", cex=1.2,cex.lab=1.3)

for (i in outliers){
if(yy[i]==1){ lines(ttt, X[i,] , type = "l", col = "black", 
          lty = 2,lwd=3)
	}
} 
 
text(19, X[indices[10],19]+50, Fechas[10], cex=1.2, col = "maroon")
text(18, X[indices[18],18]+50, Fechas[18], cex=1.2, col = "maroon")
  


dev.off()


nombre="datos-electricidad-grupo1-outliers-OUTHDS-Fechas-black-una.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X[yy==1,]), type = "l", col = "gray40", 
          lty = 1,lwd=2,ylim=c(0,2500), xlab="Time",ylab="Price", cex=1.2,cex.lab=1.3)

for (i in outliers){
if(yy[i]==1){ lines(ttt, X[i,] , type = "l", col = "black", 
          lty = 2,lwd=3)
	}
} 
 
text(19, X[indices[10],19]+50, Fechas[10], cex=1.2, col = "maroon")
 


dev.off()

nombre="datos-electricidad-grupo1-outliers-OUTHDS-Fechas-blackall.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X[yy==1,]), type = "l", col = "gray40", 
          lty = 1,lwd=2,ylim=c(0,2500), xlab="Time",ylab="Price", cex=1.2,cex.lab=1.3)

for (i in outliers){
if(yy[i]==1){ lines(ttt, X[i,] , type = "l", col = "black", 
          lty = 2,lwd=3)
	}
} 
 
text(19, X[indices[10],19]+50, Fechas[10], cex=1.2, col = "maroon")
text(18, X[indices[18],18]+50, Fechas[18], cex=1.2, col = "maroon")
 
text(22, X[indices[2],19]+49, Fechas[2], cex=1.2, col = "maroon")


text(15, X[indices[22],18]+50, Fechas[22], cex=1.2, col = "maroon")


dev.off()


nombre="datos-electricidad-grupo0-outliers-OUTHDS-Fechas-black.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X[yy==0,]), type = "l", col = "gray40", 
          lty = 1,lwd=2,ylim=c(0,2500), xlab="Time",ylab="Price", cex=1.2,cex.lab=1.3)


for (i in outliers){
if(yy[i]==0){ lines(ttt, X[i,] , type = "l", col = "black", 
          lty = 2,lwd=3)
	}
} 

text(12, X[indices[6],12]+50, Fechas[6], cex=1.2, col="maroon")

text(12, X[indices[8],12]+50, Fechas[8], cex=1.2, col="maroon")

dev.off()


nombre="datos-electricidad-grupo1-outliers-OUTHDS-black.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X[yy==1,]), type = "l", col = "gray40", 
          lty = 1,lwd=2,ylim=c(0,2500), xlab="Time",ylab="Price", cex=1.2,cex.lab=1.3)



for (i in outliers){
if(yy[i]==1){ lines(ttt, X[i,] , type = "l", col = "black", 
          lty = 2,lwd=3)
	}
} 
  
dev.off()
 

nombre="datos-electricidad-grupo0-outliers-OUTHDS-black.pdf"
pdf(nombre, bg='transparent')

matplot(ttt,t(X[yy==0,]), type = "l", col = "gray40", 
          lty = 1,lwd=2,ylim=c(0,2500), xlab="Time",ylab="Price", cex=1.2,cex.lab=1.3)



for (i in outliers){
if(yy[i]==0){ lines(ttt, X[i,] , type = "l", col = "black", 
          lty = 2,lwd=3)
	}
} 
  
dev.off()

# 



archivo <- "German-Electricity-2024-11-04-OUTHDS.RData"
save.image(archivo)

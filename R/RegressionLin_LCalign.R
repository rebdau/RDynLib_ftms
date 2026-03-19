#Regression between pos and neg chromatograms obtained for the same
#sample and run on the same instrument should be a straight curve. 
#Any peak in pos mode should have the same retention time in neg mode.
#By looking at the result from different experiments, robust regression 
#was shown to be more accurate than traditional least squares. 

RegressionLin_LCalign<-function(LCal){
	y=as.numeric(LCal[,4]) # pos retention times
	x=as.numeric(LCal[,2]) # neg retention times
	library(MASS)
	out<-rlm(y~x)
	rg.intercept<-summary(out)$coef[1,1]
	rg.slope<-summary(out)$coef[2,1]
	rg<-c(rg.intercept,rg.slope)
	return(rg)
}

PlotLin_LCalign<-function(LCal,rg){
	LCal.FT<-as.numeric(LCal[,2]) # neg retention times
	LCal.Syn<-as.numeric(LCal[,4]) # pos retention times
	plot(LCal.FT~LCal.Syn,xlab="FT pos",ylab="FT neg",mgp=c(2,0.5,0),col=5,cex.axis=0.8)
	min.FT=max(LCal.FT)/1000
	x.p=seq(min.FT,max(LCal.FT),min.FT)
	y.p=rg[1]+rg[2]*x.p
	lines(y.p,x.p) 
	pst<-paste("IWLS: Pos=",round(rg[1],digits=2),"+",round(rg[2],digits=2),"Neg")
	title(main="FTneg vs FTpos retention times (min)",sub=pst,cex.main=1,cex.sub=0.8)
}

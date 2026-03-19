Dynam_plotPie<-function(LCal,rg,z,v,pres1,pres2){
	PlotPie_LCalign(LCal,rg)
	iso.ft<-round(pres1[,1],digits=2)
	iso.syn<-round(pres2[,1],digits=2)
	for (i in 1:length(iso.ft)){
	 abline(h=iso.ft[i],lty=2,col=4)
	}
	for(i in 1:length(iso.syn)){
	 abline(v=iso.syn[i],lty=2,col=4)
	}
	points(pres2[v,1],pres1[z,1],col=2,cex=1.2,pch=21)
}
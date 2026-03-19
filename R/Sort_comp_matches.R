Sort_comp_matches<-function(pres,Assoc,time.ftng,COMPID.ftng,rg){
	if(dim(pres)[1]==1){
	 Assoc[COMPID.ftng,2]<-pres[1,3]
		# COMPID.ftng refers to COMPID entry, and is used here to
		# select the row (!) in the Assoc dataframe rather than the 
		# row corresponding with the COMPID itself.   
	}else{
	 pred_time=rg[1]+rg[2]*time.ftng
	 timediff<-abs(pres[,1]-pred_time)
	 pres.corrected<-cbind(pres,timediff)
	 o<-order(pres.corrected[,5])
	 pres.corr<-pres.corrected[o,]
	 Assoc[COMPID.ftng,2]<-pres.corr[1,3]
	}
	return(Assoc)
}

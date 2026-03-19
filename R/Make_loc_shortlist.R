Make_loc_shortlist<-function(COMPID,reg,ft.exp){
	#order on retention time
	ft.o<-ft.exp[order(ft.exp[,1]),] #reg is equal to nr of rows in ft.o or 1 less
	i<-which(ft.o[,3]==COMPID) #middle of chromatogram
	regdiv<-reg/2
	i.min<-i-regdiv
	if(i.min<0){
	 stop("Sorry,your compound is eluting very early. You have to decrease 'reg'")
	}
	i.max<-i+regdiv
	if(i.max>dim(ft.o)[1]){
	 stop("Sorry,your compound is eluting very late. You have to decrease 'reg'")
	}
	rng<-c(i.min:(i-1),(i+1):i.max) 
	ft.sh<-ft.o[rng,]	#all entries except peak in middle of chromatogram
	return(ft.sh)
}


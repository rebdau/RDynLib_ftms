CommonDotProd<-function(ac,bd){
	y<-intersect(ac[,1],bd[,1])
	if (length(y)==0){
	 common.1<-0
	 theta.1<-0
	}else{
	 zad<-length(y)+1
	 v<-rep(0,length(y))
	 w<-rep(0,length(y))
	 i=1
	 repeat{
	  # sometimes multiple product ions might be present in a MS/MS spectrum
	  # that have the same nominal mass: only the product ion with the highest
	  # relative intensity is considered
	  v[i]<-max(ac[which(ac[,1]==y[i]),2])
	  w[i]<-max(bd[which(bd[,1]==y[i]),2])
	  i=i+1
	  if (i==zad) break
	 }
	 common.1<-length(y)
	 theta.1<-sum(v*w)/(sqrt(sum(v*v))*sqrt(sum(w*w)))
	}
	return(list(common.1,theta.1))	
}
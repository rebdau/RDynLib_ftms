Find_isomers<-function(dbkey,SubDB,syn.exp,ft.exp,err=NULL){
	if(is.null(err)) err=0.02
	ft.o<-ft.exp[order(ft.exp[,2]),]
	syn.o<-syn.exp[order(syn.exp[,2]),]
	pres1<-array(dim=c(0,4)) # contains FT isomers
	pres2<-array(dim=c(0,4)) # contains QTOF isomers
	if((SubDB=="FTneg")|(SubDB=="FTpos")){
	 sl<-which(ft.exp[,3]==dbkey)
	 lb=ft.exp[sl,2]-err	# err is large!
	 ub=ft.exp[sl,2]+err
	 j=1
	 while (j<=dim(ft.o)[1]) {
	  if (ft.o[j,2]>lb&ft.o[j,2]<ub){
	   int<-ft.o[j,1:4]
	   pres1<-rbind(pres1,int)
	  }
	  if (ft.o[j,2]>=ub) break
	  j=j+1
	 }
	 j=1
	 while (j<=dim(syn.o)[1]) {
	  if (syn.o[j,2]>lb&syn.o[j,2]<ub){
	   int<-syn.o[j,1:4]
	   pres2<-rbind(pres2,int)
	  }
	  if (syn.o[j,2]>=ub) break
	  j=j+1
	 }
	}else if ((SubDB=="QTOFneg")|(SubDB=="QTOFpos")) {
	 sl<-which(syn.exp[,3]==dbkey)
	 lb=syn.exp[sl,2]-err	# err is large!
	 ub=syn.exp[sl,2]+err
	 j=1
	 while (j<=dim(syn.o)[1]) {
	  if (syn.o[j,2]>lb&syn.o[j,2]<ub){
	   int<-syn.o[j,1:4]
	   pres2<-rbind(pres2,int)
	  }
	  if (syn.o[j,2]>=ub) break
	  j=j+1
	 }
	 j=1
	 while (j<=dim(ft.o)[1]) {
	  if (ft.o[j,2]>lb&ft.o[j,2]<ub){
	   int<-ft.o[j,1:4]
	   pres1<-rbind(pres1,int)
	  }
	  if (ft.o[j,2]>=ub) break
	  j=j+1
	 }
	}else{
	 return("SubDB doesn't exist.")
	}
	return(list(pres1,pres2))
}

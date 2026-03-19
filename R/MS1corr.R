MS1corr<-function(fc,proc,namMS1){
	if ((length(fc)+4)!=dim(proc)[2]) print("fc.txt and nodes.txt are not compatible")
	sl<-which(proc[,3]%in%namMS1)
	procsl<-proc[sl,]
	nam.ind<-as.integer(rep(NA,length(namMS1)))
	for (i in 1:length(namMS1)){
	 sl<-which(procsl[,3]%in%namMS1[i])
	 nam.ind[i]<-sl
	}
	procsl<-t(procsl[nam.ind,5:dim(procsl)[2]])
	procslf<-data.frame(fc,procsl)
	corrcoef<-as.numeric(rep(NA,length(namMS1)))
	i=2
	repeat{
	 corrVar<-as.numeric(rep(NA,max(fc)))
	 j=1
	 while(j<=max(fc)){
	  procslfs<-procslf[procslf$fc==j,]
	  if((length(which(!is.na(procslfs[,2])))<3)|
			(length(which(!is.na(procslfs[,i])))<3)){
	   corrVar[j]<-NA
	  } else {
	   corrVar[j]<-cor(procslfs[,2],procslfs[,i],
					use="pairwise.complete.obs")
	  }
	  j=j+1
	 }
	 corrcoef[i-1]<-round(mean(corrVar,na.rm=T),2)	
	 if(i==dim(procslf)[2]) break
	 i=i+1
	}
	return(corrcoef)
}
File_CSPPname<-function(inp.x){
	msdet<-data.frame(inp.x[,1],inp.x[,3],inp.x[,11])
	msdet[,2]<-as.character(msdet[,2])
	msdet[,3]<-as.character(msdet[,3])
	colnames(msdet)<-c("COMPID","COMPNAME","CONVERSION")
	i=1
	while (i<=dim(msdet)[1]) {
	 if (is.na(msdet[i,2])) msdet[i,2]<-"NULL"
	 i=i+1
	}
	i=1
	while (i<=dim(msdet)[1]) {
	 if (is.na(msdet[i,3])) msdet[i,3]<-"NULL"
	 i=i+1
	}
	return(msdet)
}

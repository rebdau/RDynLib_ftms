GeneralizeSubgroupXCMS<-function(XCMS){
  #Ahlam: Here I modified the columns access
  CON.col <- which(colnames(XCMS) == "CON.new")
  SUBGRP.col <- which(colnames(XCMS) == "subgrp")
  
	o<-order(XCMS[,CON.col],XCMS[,SUBGRP.col])
	XCMS<-XCMS[o,]
	CON<-XCMS[,CON.col]
	SUBGRP<-XCMS[,SUBGRP.col]
	tim<-length(CON)-1
	CON.new<-c(2,rep(NA,tim))
	a=2
	i=2
	while(i<=length(CON)){
	 if(CON[i]==CON[i-1]){
	  if(SUBGRP[i]==SUBGRP[i-1]){
	   CON.new[i]<-a
	  }else{
	   a=a+1
	   CON.new[i]<-a
	  }
	 }else{
	  a=a+1
	  CON.new[i]<-a
	 }
	 i=i+1
	}
	XCMS<-data.frame(XCMS,CON.new)
	return(XCMS)
}

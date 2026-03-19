out.targMS2comp<-function(dbkey1,dbkey2,subDB,AnalMS){
	out<-targMS2comp(dbkey1,dbkey2,subDB,AnalMS)
	outcm<-((out$forward_ions==1)|(out$reverse_ions==1))&(out$prec_nr_ions>1)&(out$cmp_nr_ions>1)
	if (is.na(outcm)){
	 return(NA)
	}else{
	 if (outcm==TRUE){
	  if (is.nan(out$dot_ions)){
	   return(0)
	  }else{
	   if (out$dot_ions>0.9){
	    return(out$dot_ions)
	   }else{
	    return(NA)
	   }
	  }
	 }else{
	  return(NA)
	 }
	}
}

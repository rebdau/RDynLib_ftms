ProdIonMatch<-function(prod_ion.num,prdion,err){
	i=1
	repeat{
	 lowmz<-prdion[i,1]-err
	 highmz<-prdion[i,1]+err
	 j=1
	 repeat{
	  if((prod_ion.num[j]>=lowmz)&(prod_ion.num[j]<=highmz)){
	   infostr<-paste("m/z",prod_ion.num[j],"--",prdion[i,2],
				"--",prdion[i,3],sep=" ")
	   print(infostr)
	  }
	  if(j==length(prod_ion.num))break
	  j=j+1
	 }
	 if(i==dim(prdion)[1])break
	 i=i+1
	}
	writeLines("\n")
}


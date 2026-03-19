# take care that compound_add.txt is present in the working directory
# ave.cnt is the average of common neutrals and common ions. If half of the
# product ions are common and the other half yield neutral losses in common,
# than the ave.cnt is 0.5! Noticeably, 0.5 would also be obtained if that half
# of the product ions that are in common are also responsible for the common
# neutral losses.
# ave.dot refers to the average dot product obtained via either the common
# product ions or the common neutral losses.
# Comment: a better dot product might be based on all product ions of a parti-
# cular spectrum whether present or not in the other spectrum???? This would
# take the count into account in the dot product. Two possibilities are then a
# forward and reverse search.
# CSPPs are ordered first on decreasing number of ions in common, then on ave.cnt
# and finally on ave.dot. For each substrate COMPID, only one product COMPID is
# retained, i.e. the one that has the "highest similarity", namely that is the
# highest on the ordered CSPP.df list.
# Other CSPPs for the particular conversion type, e.g. methylation,
# that contain either the same substrate COMPID or the same product COMPID
# are removed from the cspp.df list before continuation. If you want to keep a
# complete record of all generated CSPPs, write the initial cspp.df file away to
# a .txt file.
# the conversion type allows to determine the column (indicated by conv.col) of
# comp_app (and finally of compound_add.txt) to which the results need to be written. 

rank.cspp<-function(cspp.df,conv.col,comp_add){
	ave.cnt<-(cspp.df[,12]+cspp.df[,14])/2
	ave.dot<-(cspp.df[,8]+cspp.df[,10])/2
	cspp.df<-data.frame(cbind(cspp.df,ave.cnt,ave.dot))
	cspp.df<-cspp.df[order(-cspp.df[,6],-cspp.df[,15],-cspp.df[,16]),]
	while(dim(cspp.df)[1]!=0){
	 compid.sub<-cspp.df[1,1]
	 compid.prod<-cspp.df[1,4]
	 cspp.name<-paste("!!",cspp.df[1,6],"!",
				cspp.df[1,15],"!",cspp.df[1,16],
				"!!",compid.prod,sep="")
	 comp_add[compid.sub,conv.col]<-cspp.name
	 res.prod<-which(cspp.df[,4]%in%compid.prod)
	 res.sub<-which(cspp.df[,1]%in%compid.sub)
	 cspp.df<-cspp.df[-sort(union(res.sub,res.prod)),]
	}
	return(comp_add)
}

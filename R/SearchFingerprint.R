# searching CID spectra in DynLib experiment for those COMPID that have the 
# MS spectral fingerprint represented by the output from CIDfingerprint()

SearchFingerprint<-function(base.dir,finlist,SubDB,exp,common_ion){
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB<-TMP[[1]]
	AnalMS<-TMP[[2]]
	subDB_sel<-subDB[[1]][subDB[[1]][,4]==exp,]
	COMPIDstart<-subDB_sel[1,3]
	COMPIDstart_row<-which(subDB[[1]][,3]%in%COMPIDstart)
	COMPIDend<-subDB_sel[dim(subDB_sel)[1],3]
	COMPIDend_row<-which(subDB[[1]][,3]%in%COMPIDend)
	MS2list<-subDB[[3]][COMPIDstart_row:COMPIDend_row]
	unk<-common_ion
	minIons<-length(unk)-2
	msnSect<-MultMatch(MS2list,unk,minIons)
	ENTRIES<-subMultMatch(msnSect,unk,MS2list,minIons)
	nrCOMPID_rows<-ENTRIES+COMPIDstart_row-1
	nrCOMPID<-subDB[[1]][nrCOMPID_rows,3]
	return(nrCOMPID)
}
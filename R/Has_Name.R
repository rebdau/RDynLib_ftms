Has_Name<-function(name.cha,No.name,sel){
	HasName<-as.integer()
	j=1
	while (j<=length(name.cha)) {
	 Ent.name<-as.character(c(substr(name.cha[j],1,1),name.cha[j]))
	 if (length(which(Ent.name%in%No.name))!=0){ # No.name is NULL or !
	 # Ent.name contains first character of the name to check for "!"
	 # and full name to check for "NULL"
	  HasName<-append(HasName,0)
	 }else{
	  HasName<-append(HasName,sel[j])
	 }
	 j=j+1
	}
	return(HasName)
}

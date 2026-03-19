SearchAdduct<-function(mtr,buf,dt.s,dt,s2,anno.1) {
	t=0	
	repeat{
		t=t+1
		nmp<-dt.s[mtr[mtr$z3==buf,]$rw[t],1]	# (A)
		anno.1[which(dt[,1]==nmp)]="adduct"
		if (t==s2) break
	}
	return(anno.1)
}

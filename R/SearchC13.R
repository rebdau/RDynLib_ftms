SearchC13<-function(mtr,dt.s,dt,s1,anno,findmn){
	t=0
	anno.t<-as.character(c())
	repeat{
		t=t+1
		nms<-dt.s[mtr[mtr$z3==1.0,]$cw[t],1]		
			# 'substrate' name, i.e. peak name of the candidate C12 (A)
		mngs<-dt.s[mtr[mtr$z3==1.0,]$cw[t],findmn]		
			# amount of candidate C12 peak				(A)
		nmp<-dt.s[mtr[mtr$z3==1.0,]$rw[t],1]			
			# 'product' name, i.e. peak name of candidate C13	(A)
		mngp<-dt.s[mtr[mtr$z3==1.0,]$rw[t],findmn]		
			# amount of candidate C13 peak				(A)
		if (mngs>mngp) {						
			# candidate C13 is only annotated as 'C13' if its amounts 
			#are lower than that of the candidate C12, otherwise it is 
			#annotated as 'C12'
			anno[which(dt[,1]==nmp)]="C13"
		}else{
			anno[which(dt[,1]==nmp)]="C12"
		}	
		anno.t<-anno[which(dt[,1]==nms)]
			# candidate C12 wordt geannoteerd als 'C12' enkel als er 
			#al geen annotatie als 'C13' eerder gebeurd was 
		if (anno.t=="C13"){
			if (t==s1) break
			next
		}else{
			anno[which(dt[,1]==nms)]="C12"
		}
		if (t==s1) break
	}
	return(anno)
}

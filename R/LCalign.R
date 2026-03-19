#LCalign will take 5 neighboring peaks (in the case t.ini=5) and compute all
#possible retention differences of each of the 5 FT individual m/z values with 
#all possible isomeric QTOF m/z values. For each combination of 5 FTm/z,QTOFm/z
#pairs, the standard deviation on the retention time differences is computed.
#The combination that yields the lowest standard deviation is selected and the
#5 FTm/z,QTOFm/z pairs added to a final list (res.F). Then this local sequence of 5 
#neighboring FT peaks is shifted with one peak towards the right of the
#chromatogram and the procedure repeated. For this "local retention time alignment"
#procedure, the LCalign function closely interacts with the Differ_of_retdif
#function. As each FT peak will enter in such a 'local neighboring peak sequence'
#5 times (first as the fifth peak, then as the fourth peak,....), each FT peak will
#yield 5 FTm/z,QTOFm/z pairs in the final list (res.F). Lets say that a particular
#FTm/z peak has multiple isomeric QTOFm/z peak options to which it might be paired
#and that always the same QTOFm/z isomer is taken in the 5 local sequences, this
#gives strong support that indeed the right FTm/z,QTOFm/z peak association was made.
#As this peak pair combination will be 5 times added to the res.F final list, they
#will contribute strongly to the subsequent regression following LCalign.
#In contrast, if the same FTm/z peak is not always associated with the same QTOFm/z
#peak, the regression will be slightly biased. This can happen if you have a
#sequence of e.g. 7 neighboring FT peaks and e.g. the fourth peak has two possible
#QTOFm/z isomers yielding a retention time difference of e.g. 0.5 min and 1 min. If
#the first three FTm/z peaks differ by 0.5 min with their respective QTOFm/z isomer,
#their combinations with the fourth FTm/z peak will lead to selection of the first
#QTOFm/z isomer for the latter FTm/z peak to obtain a low sd value. If the last
#three FTm/z peaks differ by 1 min with their respective QTOFm/z isomers, their
#combinations with the fourth FTm/z peak will lead to selection of the second
#QTOFm/z isomer for this fourth FTm/z peak. One could also consider the multiple
#contribution of each FTm/z peak to the final res.F list as a weighing effect:
#assume that a good association between the FTm/z peak and one of the QTOFm/z
#isomers is not possible, then adding multiple associations for the FTm/z peak will
#yield a less biased regression than if only one, false positive, association would
#have been added. Of course, the first and the last few
#FT chromatogram peaks cannot contribute 5 times to the final res.F list and, thus,
#will contribute less to the regression. If only one QTOFm/z isomer is present,
#this one will be selected to form a FTm/z,QTOFm/z peak pair. In case the latter
#association is wrong, this will yield a false positive association. Such false
#positive associations are expected to contribute prominently to the final res.F
#list, yet they are expected to deviate randomly with regard to the finally obtained
#regression curve following the LCalign procedure. Pay attention that at least one
#QTOFm/z isomer will be present due to the prior application of the Remov_empties
#function. 
  


LCalign<-function(err,t.ini,syn.exp,ft.sh){
	o<-order(syn.exp[,2]) # order on m/z value
	syn.o<-syn.exp[o,] 
	res.F<-array(dim=c(0,7))
	Lst<-list()
	a=1
		#do the following while loop only for the e.g. 4 (depends on
		#t.ini-1) lowest m/z values.	
	while (a<t.ini){  # one less than the number of adjacent peaks
	 A<-ft.sh[a,2]	# A is FTm/z value
	 Lst[[a]]<-array(dim=c(0,7))
	 j=1			# rows in QTOF DynLib experiment
	 repeat{
	  if (syn.o[j,2]>A-err&syn.o[j,2]<A+err){
		#if QTOF row-based m/z value falls into the FT-selected m/z window:
		#compute retention time difference between QTOF m/z and FT m/z
		#add all info, i.e. FT COMPID, FT RT, FT m/z, QTOF RT, QTOF m/z,
		#retention time difference, QTOF COMPID, to the Lst list at
		#position [[a]]. Lst[[a]] is a data.frame containing all isomeric
		#QTOF m/z values corresponding to the selected FT m/z value.
	   retdiff.A<-ft.sh[a,1]-syn.o[j,1]
	   crd.A<-c(as.character(ft.sh[a,3]),ft.sh[a,1],ft.sh[a,2],syn.o[j,1],
			syn.o[j,2],retdiff.A,syn.o[j,3])
	   Lst[[a]]<-rbind(Lst[[a]],crd.A)
        }
		#Stop searching the QTOF m/z values if you reach the row where
		#the QTOF m/z value is bigger than the upper threshold of the
		#selected FT m/z value.
	  if ((syn.o[j,2]>A+err)|(j==dim(syn.o)[1])) break
        j=j+1
	 }
	 a=a+1
	}
	i=t.ini
		#do the next repeat loop for all m/z values except for the
		#lowest e.g. 4 (depends on t.ini-1) m/z values
	repeat{
	 C<-ft.sh[i,2]	#C is FTm/z value starting from the 't.ini' peak onwards
	 Lst[[t.ini]]<-array(dim=c(0,7))
	 j=1
	 repeat{
		#see while loop above
	  if (syn.o[j,2]>C-err&syn.o[j,2]<C+err){
	   retdiff.C<-ft.sh[i,1]-syn.o[j,1]
	   crd.C<-c(as.character(ft.sh[i,3]),ft.sh[i,1],ft.sh[i,2],syn.o[j,1],syn.o[j,2],retdiff.C,syn.o[j,3])
	   Lst[[t.ini]]<-rbind(Lst[[t.ini]],crd.C)
	  }
	  if ((syn.o[j,2]>C+err)|(j==dim(syn.o)[1])) break
	  j=j+1
	 }
		#the next if clausule should never happen as Remov_empties resulted in
		#the elimination of those FTm/z peaks that do not have at least one 
		#QTOFm/z isomer and so, the else clausule should be entered.
	 if (is.null(Lst[[t.ini]])){
	  le<-length(Lst)-1
	  Lst<-Lst[1:le]
#	  i=i-1 #
	 }else{
		#compute the best local retention time alignment among all possible
		#pairs for these 5 FTm/z peaks
	  res.F<-Differ_of_retdif(Lst,res.F)
		#for t.ini=5, the length of Lst remains always 5. Now, the second Lst
		#entry becomes the first entry, the third Lst entry becomes the second,
		#the fourth becomes the third and the fifth entry becomes the fourth.
		#this is a for loop with 4 iterations:
	  le<-length(Lst)-1
	  for (b in 1:le){ 
	   Lst[[b]]<-Lst[[b+1]]
	  }
		#remove the fifth entry of Lst so that a new Lst[[t.ini]] entry can be
		#added and a next round of local retention alignment can be performed:
	  Lst<-Lst[1:le]
	 }
	 if (i==dim(ft.sh)[1]) break
	 i=i+1	#next peak in the FT chromatogram
	}
	return(res.F)
}

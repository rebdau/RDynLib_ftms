#' @title Breaking featureDefinitions() rows.
#'
#' @description
#' creates list of breaking rows to enable the subsequent peakgrouping of a
#' xcmsExperiment object. 
#'
#' @param tR (e.g. 1 sec) defines the retention time difference between two 
#'       subsequent rows that is needed to define a breaking row
#'
#' @param dat the output of the featureDefinitions() of the xcmsExperiment 
#'        object
#' 
#' @param rtmed the median retention times extracted from the 
#'        featureDefinitions() function.
#'        
#'@return a list of indices of each breaking row.
#'
#'@author Ahlam Mentag
#'
#'@export
BreakRowsXCMS_SQL <- function(dat, tR, rtmed){
  
  o<-order(rtmed)
  dat.o<-dat[o,]
  mx<-dim(dat.o)[1]
  if(mx<=1000){
    stp<-c(1,mx)
  }else{
    stp<-c(0)
    i=500
    while(i+500<mx){
      repeat{
        j=i+1
        ch.1<-dat.o[i,8]
        ch.2<-dat.o[j,8]
        ch=ch.2-ch.1
        if(ch>=tR){					
          stpm<-i
          stp<-append(stp,stpm)
          i=i+500
          break
        }else{
          i=i+1
          if(i==mx)break
        }
      }
    }
  }
  stp<-append(stp,mx)
  return(stp)
}
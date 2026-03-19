#' @title Grouping of the features.
#'
#' @description
#'  the `AutomFeatureGroup_SQL()` function sort the feature information by 
#'  retention time and then calls the featureGroupsXcms() function to perform 
#'  the grouping instructions inside it.
#'  
#'
#' @param tR (e.g. 1 sec) defines the retention time difference between two 
#'       subsequent rows that is needed to define a breaking row
#'
#' @param dat the output of the featureDefinitions() of the xcmsExperiment 
#'        object.
#'        
#' @param ints the intensities of the features in each sample.
#' 
#' @param rtmed the median retention times extracted from the 
#'        featureDefinitions() function.
#'        
#' @param stp set within the XCMSgrouping() function as the output from the 
#'        BreakRowsXCMS() function: represents the indices of features with a 
#'        particular retention time on which to split.
#'    
#' @param y10 maximum retention time deviation for feature alignment across 
#'  chromatograms, 1 sec by default.
#' 
#' @param y20 minimum feature abundance-based Pearson correlation for feature  
#'  alignment across chromatograms, 0.8 by default.
#'  
#' @param fc set within the XCMSgrouping() function as the output from the 
#'  LoadProc() function: represents the genotype/treatment class information.
#'
#'@return an XCMSExperiment object with the new grouping columns.
#'
#'@author Ahlam Mentag
#'
#'@export
AutomFeatureGroup_SQL <-function(dat,rtmed,ints,stp,y10,y20,fc){
  stp<-as.integer(stp)

  o<-order(rtmed)
  dat.o<-dat[o,]
  XCMS<-data.frame()
  i=1
  repeat{
    j=i+1
    strt<-stp[i]+1
    idx <- o[strt:stp[j]]      # row indices in original dat
    dt  <- dat[idx, , drop = FALSE]
    ints.sub <- ints[idx, , drop = FALSE]
    
    xcms2 <- FeatureGroupXCMS_SQL(dt, ints.sub, y10, y20, fc)
    

    
    XCMS<-rbind(XCMS,xcms2)	
    if(j==length(stp))break
    i=i+1
  }
  return(XCMS)
}

#' @title Further grouping of the features.
#'
#' @description
#'  the `FeatureSubgroupXCMS_SQL()` function performs further grouping based 
#'  on the median retention time of each current group
#'  
#'
#' @param XCMS an `XCMSExperiment` object with the previous grouping columns
#'    
#' @param y10 maximum retention time deviation for feature alignment across 
#'  chromatograms, 1 sec by default.
#'  
#'@return an XCMSExperiment object with the new grouping columns.
#'
#'@author Ahlam Mentag
#'
#'@export
FeatureSubgroupXCMS_SQL <- function(XCMS, y10) {
  
  group_col <- colnames(XCMS)[ncol(XCMS)]
  grp_vals  <- XCMS[[group_col]]
  
  ## Preallocate result
  XCMS$subgrp <- integer(nrow(XCMS))
  
  ## Loop over groups ONCE
  for (i in unique(grp_vals)) {
    
    idx <- which(grp_vals == i)
    if (length(idx) == 0) next
    
    rt <- XCMS$rtmed[idx]
    tRmed <- median(rt)
    
    low  <- tRmed - y10
    high <- tRmed + y10
    
    ## Vectorized assignment 
    subgrp <- ifelse(
      rt < low,  1L,
      ifelse(rt > high, 3L, 2L)
    )
    
    XCMS$subgrp[idx] <- subgrp
  }
  
  return(XCMS)
}




#' @title featureDefinitions Extraction and sample numbers to a list.
#'
#' @description 
#'      'LoadProc_SQL()' function creates list of featureDefinitions and sample
#'       numbers. 
#'
#' @param x `XcmsExperiment` object with the results from xcms preprocessing.
#'        
#' @return a list of featureDefinitions() and fc (sample number).
#'
#' @author Ahlam Mentag
#'
#' @export
LoadProc_SQL <- function(x) {
  
  ## proc replaces nodes.txt, it is now the feature definitions table
  proc <- featureDefinitions(x)
  
  ## fc replaces fc.txt, it is extracted from sample metadata
  fc <- sampleData(x)$fc
  
  return(list(proc, fc))
}

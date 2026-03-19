#' @title Annotate peaks in XCMS feature groups
#'
#' @description
#'  the `FeatureAnnotXCMS_SQL()` performs the annotation of the features within 
#'  each FG as 13C isotopes ‘anno’ column , buffer adducts ‘anno.1’ column
#'  or ion/neutral complex between the molecule and its ion ‘anno.2’ column.
#'                                                                         
#' @param ints the intensities of the features in each sample.
#'    
#' @param XCMS an `XCMSExperiment` object with the previous grouping columns
#' 
#' @param y11 monoisotopic mass of the buffer adduct, 46 Da (formic acid) 
#'  by default.
#'  
#' @param fc set within the XCMSgrouping() function as the output from the 
#'  LoadProc() function: represents the genotype/treatment class information.
#'
#'@return an XCMSExperiment object with the new annotation columns.
#'
#'@author Ahlam Mentag
#'
#'@export
FeatureAnnotXCMS_SQL <- function(XCMS, y11, fc, ints) {
  
  buf <- y11
  dt <- XCMS
  anno <- anno.1 <- anno.2 <- rep("ann", nrow(dt))
  
  # Determine the range of feature groups
  if(!"CON.new" %in% colnames(dt)) stop("dt must have CON.new column")
  mx <- max(dt$CON.new)
  f <- min(dt$CON.new) - 1  # start at first group
  
  repeat {
    f <- f + 1
    if (f > mx) break
    
    dt.s <- dt[dt$CON.new == f, ]
    if (nrow(dt.s) < 2) next  # skip singleton groups
    
    # Ensure features exist in ints
    feature_rows <- rownames(dt.s) %in% rownames(ints)
    if (!any(feature_rows)) next
    dt.s <- dt.s[feature_rows, , drop = FALSE]
    
    # Compute mean intensity per feature safely
    mng <- apply(ints[rownames(dt.s), , drop = FALSE], 1, mean, na.rm = TRUE)
    if (length(mng) != nrow(dt.s)) next
    dt.s$mng <- mng
    
    # Mass difference table (nominal differences)
    k <- nrow(dt.s)
    z3 <- outer(dt.s$mzmed, dt.s$mzmed, "-")
    z3 <- z3[lower.tri(z3)]
    z3 <- round(z3, 1)
    
    # Pairwise indices
    rw <- unlist(lapply(1:(k-1), function(i) seq(i+1, k)))
    cw <- unlist(lapply(1:(k-1), function(i) rep(i, k-i)))
    mtr <- data.frame(cw = cw, rw = rw, z3 = z3)
    if (nrow(mtr) == 0) next
    
    # Dimer annotation
    for (u in seq_len(nrow(mtr))) {
      cw_mean <- dt.s[mtr$cw[u], "mng"]
      if (length(cw_mean) == 0) next
      cw_mean <- cw_mean[1]  # force scalar
      if (!is.na(cw_mean) && round(mtr$z3[u], 1) == round(cw_mean + 1, 1)) {
        nmp <- rownames(dt.s)[mtr$rw[u]]
        anno.2[which(rownames(dt) == nmp)] <- "dimer"
      }
    }
    
    # C13 / adduct annotation 
    s1 <- sum(mtr$z3 == 1.0)
    s2 <- sum(mtr$z3 == buf)
    
    if (s1 > 0) {
      anno <- tryCatch(
        SearchC13(mtr, dt.s, dt, s1, anno, length(fc) + max(fc) + 17),
        error = function(e) {
          message("[FeatureAnnotXCMS_SQL] Warning: SearchC13 skipped for group ", f)
          return(anno)
        }
      )
    }
    if (s2 > 0) {
      anno.1 <- tryCatch(
        SearchAdduct(mtr, buf, dt.s, dt, s2, anno.1),
        error = function(e) {
          message("[FeatureAnnotXCMS_SQL] Warning: SearchAdduct skipped for group ", f)
          return(anno.1)
        }
      )
    }
  }
  
  dt$anno <- anno
  dt$anno.1 <- anno.1
  dt$anno.2 <- anno.2
  return(dt)
}

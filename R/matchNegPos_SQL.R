#' @title Match  MS2 spectra between the same instruments type 
#'        using SQL databases. 
#'
#' @description
#' 
#' `matchNegPos_SQL()` Removes false positive negFTm/z,posFTm/z peak pairs by 
#                      checking whether a certain percentage of the neg MS2  
#                      spectral peaks can be traced in the pos MS2 spectrum. 
#                      By default, 20% of the neg MS2 peaks are traced in the
#                      pos MS2 spectrum.
#'
#' @param LCal 
#' `data.frame` returned by `Aligning_General_SQL()` function.
#' 
#' @param FTn_con `DBIConnection`  
#'   A `DBI connection` object to the reference SQLite database.
#'
#' @param FTp_con `DBIConnection`  
#'   A `DBI connection` object to the SQLite database to align with.
#'
#' @param polarity_ftn 
#'        `Integer (0/1)` Polarity of negative FTMS or QTOF experiment.
#'        
#' @param polarity_ftp 
#'        `Integer (0/1)` Polarity of positive FTMS or QTOF experiment.
#' 
#' @param minPeaks `Numeric(1)` Minimum matching ratio (default 0.2).
#' 
#'
#' @return Filtered `LCal` with MS2-supported neg–pos pairs.
#'
#' @import DBI RSQLite
#' @export
matchNegPos_SQL <- function(
    LCal,
    FTn_con,
    FTp_con,
    polarity_ftn = 0,
    polarity_ftp = 1,
    minPeaks = 0.2
) {
  
  ## Load negative-mode compound table
  subdb_neg <- dbGetQuery(
    FTn_con,
    "SELECT retention_time, mass_measured, compound_id
     FROM ms_compound ORDER BY compound_id"
  )
  

  ## Load negative-mode MS2 peaks
  neg_cols <- dbListFields(FTn_con, "msms_spectrum")
  
  
  ms2_neg_df <- dbGetQuery(FTn_con, sprintf("
    SELECT s.compound_id, p.mz
    FROM msms_spectrum s
      JOIN msms_spectrum_peak p USING(spectrum_id)
    WHERE s.ms_level = 2
      AND s.polarity = %d
    ORDER BY s.compound_id, p.mz
  ", polarity_ftn))
  
  ms2_neg_list <- split(round(ms2_neg_df$mz), ms2_neg_df$compound_id)
  ms2_neg_list <- lapply(ms2_neg_list, unique)
  
  ## Load positive-mode MS2 peaks
  pos_cols <- dbListFields(FTp_con, "msms_spectrum")

  
  ms2_pos_df <- dbGetQuery(FTp_con, sprintf("
    SELECT s.compound_id, p.mz
    FROM msms_spectrum s
      JOIN msms_spectrum_peak p USING(spectrum_id)
    WHERE s.ms_level = 2
      AND s.polarity = %d
    ORDER BY s.compound_id, p.mz
  ", polarity_ftp))
  
  ms2_pos_list <- split(round(ms2_pos_df$mz), ms2_pos_df$compound_id)
  ms2_pos_list <- lapply(ms2_pos_list, unique)
  
  ## Align negative MS2 list to compound table
  ms2_neg <- ms2_neg_list[ as.character(subdb_neg$compound_id) ]
  ms2_neg[sapply(ms2_neg, is.null)] <- list(integer(0))
  
  ## Replace NULL in positive MS2 list with empty vectors
  ms2_pos_list[sapply(ms2_pos_list, is.null)] <- list(integer(0))
  

  ## Apply peak-overlap filtering
  if (is.null(LCal) || nrow(LCal) == 0) {
    message("No candidate neg–pos alignments found; returning empty LCal.")
    return(LCal)
  }
  
  i <- 1
  while (i <= nrow(LCal)) {
    
    neg_id <- LCal[i, 1]   # FT negative COMPID
    pos_id <- LCal[i, 7]   # FT positive COMPID
    
    neg_row <- which(subdb_neg$compound_id == neg_id)
    
    ## Remove pair if negative COMPID not found
    if (length(neg_row) == 0) {
      LCal <- LCal[-i, ]
      next
    }
    neg_peaks <- ms2_neg[[neg_row]]
    pos_peaks <- ms2_pos_list[[as.character(pos_id)]]
    if (length(neg_peaks) == 0 || length(pos_peaks) == 0) {
      LCal <- LCal[-i, ]
      next
    }
    
    pos_corrected <- pos_peaks - 2
    same_ion <- sum(neg_peaks %in% pos_corrected)
    
    if (same_ion / length(neg_peaks) < minPeaks) {
      LCal <- LCal[-i, ]
      next
    }
    
    i <- i + 1
  }
  
  return(unique(LCal))
}

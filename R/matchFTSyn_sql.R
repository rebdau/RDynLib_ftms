#' @title Match FT and QTOF MS2 spectra based on local alignments 
#'
#' @description
#' `matchFTSyn_SQL()` filters locally aligned FTMS–QTOF candidate pairs in
#' (`LCal`) by comparing MS2 peak lists stored in SQL databases.
#'
#' @details
#' The function:
#' - Loads FT and QTOF MS2 spectra from SQLite via DBI
#' - Applies polarity and spectrum-type filtering
#' - Computes MS2 peak overlap for each FT–QTOF candidate pair
#' - Removes pairs with less matched peaks based on a threshold
#'
#' @param LCal 
#' `data.frame` A data.frame returned by `Aligning_General_SQL()` function.
#' 
#' @param FT_con `character(1)` A DBI connection object 
#' to the FTMS SQLite database (the reference database).
#' 
#' @param QTOF_con `character(1)` A DBI connection object to the QTOF
#' SQLite database.
#' 
#' @param polarity_ft `numeric(1)` Integer (0 or 1), refers to the polarity 
#' of the ftms  experiment to align.
#' 
#' @param polarity_qtof `numeric(1)`Integer (0 or 1), refers to the polarity 
#' of the qtof experiment to align.
#' 
#' @param minPeaks `Numeric(1)` Minimum MS2 matching ratio required.
#' 
#' @param aggregated_Ft `logical(1)` Use aggregated FTMS spectra.
#' 
#' @param aggregated_QTOF `logical(1)` Use aggregated QTOFMS spectra.
#'
#' @return A filtered version of `LCal`, keeping only MS2 supported matches.
#'
#' @importFrom DBI dbConnect
#'
#' @importFrom DBI dbDisconnect
#'
#' @importFrom RSQLite SQLite
#' 
#' @importFrom DBI dbGetQuery
#' 
#' @importMethodsFrom DBI dbListFields
#' 
#' @author Ahlam Mentag
#'
#' @export
matchFTSyn_SQL <- function(LCal, FT_con, QTOF_con, polarity_ft = 0, 
                           polarity_qtof = 0,minPeaks = 0.6, 
                           aggregated_Ft = FALSE, aggregated_QTOF = FALSE) {
  
  
  # Load FT and QTOF compound tables
  subdb_ft  <- dbGetQuery(FT_con,
                          "SELECT retention_time, mass_measured, compound_id
     FROM ms_compound ORDER BY compound_id")
  
  subdb_syn <- dbGetQuery(QTOF_con,
                          "SELECT retention_time, mass_measured, compound_id
     FROM ms_compound ORDER BY compound_id")
  
  
  # Retrieve MS2 peak lists
  # FTMS MS2 peaks
  ft_cols <- dbListFields(FT_con, "msms_spectrum")
  
  if ("spectrum_type" %in% ft_cols) {
    ft_filter <- if (!aggregated_Ft) {
      "(s.spectrum_type IS NULL OR s.spectrum_type = '' OR s.spectrum_type = 'not_assembled')"
    } else {
      "s.spectrum_type = 'assembled'"
    }
  } else {
    # IF spectrum_type column then no filtering
    ft_filter <- "1=1"
  }
  
  
  ms2_ft_df <- dbGetQuery(FT_con, sprintf("
    SELECT s.compound_id, p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p USING(spectrum_id)
    WHERE s.ms_level = 2
      AND s.polarity = %d
      AND %s
    ORDER BY s.compound_id, p.mz", polarity_ft,ft_filter))
  
  
  # Round and split into lists per compound_id
  ms2_ft_list <- split(round(ms2_ft_df$mz), ms2_ft_df$compound_id)
  ms2_ft_list <- lapply(ms2_ft_list, unique)
  
  
  # QTOF MS2 peaks 
  qtof_cols <- dbListFields(QTOF_con, "msms_spectrum")
  
  if ("spectrum_type" %in% qtof_cols) {
    qtof_filter <- if (!aggregated_QTOF) {
      "(s.spectrum_type IS NULL OR s.spectrum_type = '' OR s.spectrum_type = 'intersect_single_energy')"
    } else {
      "s.spectrum_type = 'merged_MSMS_all_energies'"
    }
  } else {
    qtof_filter <- "1=1"
  }
  
  
  ms2_syn_df <- dbGetQuery(QTOF_con, sprintf("
    SELECT s.compound_id, p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p USING(spectrum_id)
    WHERE s.ms_level = 2
      AND s.polarity = %d
      AND %s
    ORDER BY s.compound_id, p.mz", polarity_qtof,qtof_filter ))
  
  
  ms2_syn_list <- split(round(ms2_syn_df$mz), ms2_syn_df$compound_id)
  ms2_syn_list <- lapply(ms2_syn_list, unique)
  
  # Align MS2 lists by compound ID order

  ms2.ft  <- ms2_ft_list[ as.character(subdb_ft$compound_id) ]
  ms2.syn <- ms2_syn_list[ as.character(subdb_syn$compound_id) ]
  
  # Replace NULL or NA lists with empty vectors
  ms2.ft[sapply(ms2.ft, is.null)]   <- list(integer(0))
  ms2.syn[sapply(ms2.syn, is.null)] <- list(integer(0))
  

  
  if (is.null(LCal) || nrow(LCal) == 0) {
    message("No candidate alignments found; returning empty LCal.")
    return(LCal)
  }
  
  i <- 1
  while (i <= nrow(LCal)) {
    
    ft.compid  <- LCal[i, 1]
    syn.compid <- LCal[i, 7]
    
    ft.row  <- which(subdb_ft$compound_id  == ft.compid)
    syn.row <- which(subdb_syn$compound_id == syn.compid)
    
    # Skip pair if compound_id not found
    if (length(ft.row) == 0 || length(syn.row) == 0) {
      LCal <- LCal[-i, ]
      next
    }
    
    ms2ion  <- ms2.ft[[ft.row]]
    msmsion <- ms2.syn[[syn.row]]
    
    # Skip if MS2 list is empty 
    if (length(ms2ion) == 0 || length(msmsion) == 0) {
      LCal <- LCal[-i, ]
      next
    }
    
    # Compute matching peak ratio
    same_ion <- sum(ms2ion %in% msmsion)
    
    if (same_ion / length(ms2ion) < minPeaks) {
      LCal <- LCal[-i, ]
      next
    }
    
    i <- i + 1
  }
  
  return(unique(LCal))
}

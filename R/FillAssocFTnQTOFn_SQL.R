#' @title Associate FT and QTOF features using SQL and MS2 matching
#'
#' @description
#' `FillAssocFTnQTOFn_SQL()` finds associations between FTMS and QTOFMS
#' compounds based on retention time alignment, MS2 peak similarity, and
#' alignment rules. It updates and returns the association table (`Assoc`).
#'
#' @details
#' The function:
#' - loads FT and QTOF compounds from SQL databases,
#' - computes expected RT alignment using regression parameters,
#' - finds candidate matches using SQL-assisted search,
#' - filters candidates using MS2 peak similarity,
#' - appends validated associations to `Assoc`.
#'
#' @param FT_con `DBIConnection`  
#'   A DBI connection object to the FT SQLite database.
#'
#' @param QTOF_con `DBIConnection`  
#'   A DBI connection object to the QTOF SQLite database.
#'
#' @param Assoc `data.frame` A table storing previously detected associations. 
#'   Will be filled with the new association if they are not already 
#'   in the assoc table.
#'
#' @param FT_expnr `numeric(1)`  
#'   The reference experiment number to load from the FTMS SQL database.
#'
#' @param QTOF_expnr `numeric(1)`  
#'   Experiment number to load from the QTOF SQL database.
#'
#' @param cutoff `numeric(1)`  
#'   Minimum retention time to consider for FT compounds.
#'
#' @param rg `numeric(6)`  
#'   Regression coefficients for retention-time alignment.
#'
#' @param lc.err `numeric(1)`  
#'   Allowed retention-time error window after regression transformation.
#'
#' @param err `numeric(1)`  
#'   Error threshold for candidate match search (`Find_cand_matches_SQL`).
#'
#' @param minIon `numeric(1)`  
#'    The minimum matching ions between spectra of peak pairs, which are found
#'    after applying the regression model, and which are the final matching 
#'    written in the Assoc table. (default `0.02`).
#'
#' @param polarity_ft `integer(1)`  
#'   Polarity filter for FTMS MS2 spectra (0 or 1).
#'
#' @param polarity_qtof `integer(1)`  
#'   Polarity filter for QTOF MS2 spectra (0 or 1).
#'
#' @param FT_path `character(1)`  
#'   File path to the FTMS database .
#'
#' @param QTOF_path `character(1)`  
#'   File path to the QTOF database.
#'
#' @param aggregated_Ft `logical(1)`  
#'   Whether to use aggregated FTMS spectra.
#'
#' @param aggregated_QTOF `logical(1)`  
#'   Whether to use aggregated QTOF MS2 spectra.
#'
#' @return `data.frame`  
#'   The updated association table with new FTMS–QTOF matches appended.
#'
#' @importFrom DBI dbConnect
#'
#' @importFrom DBI dbDisconnect
#'
#' @importFrom RSQLite SQLite
#'
#' @author Ahlam Mentag
#'
#' @export

FillAssocFTnQTOFn_SQL <- function(
    FT_con, QTOF_con, Assoc, FT_expnr, QTOF_expnr,
    cutoff, rg, lc.err, err, minIon = 0.02,
    polarity_ft = 0, polarity_qtof = 0,
    FT_path, QTOF_path,  aggregated_Ft = FALSE, aggregated_QTOF = FALSE
) {
  
  # Load FT and QTOF compounds
  ft.exp <- dbGetQuery(FT_con, sprintf(
    "SELECT retention_time, mass_measured, compound_id, expid 
     FROM ms_compound WHERE expid = %d", FT_expnr))
  
  syn.exp <- dbGetQuery(QTOF_con, sprintf(
    "SELECT retention_time , mass_measured, compound_id, expid 
     FROM ms_compound WHERE expid = %d", QTOF_expnr))
  
  # MS2 filtering based on aggregated flag
  ft_cols <- dbListFields(FT_con, "msms_spectrum")
  
  if ("spectrum_type" %in% ft_cols) {
    ft_filter <- if (!aggregated_Ft) {
      "(s.spectrum_type IS NULL OR s.spectrum_type = '' OR s.spectrum_type = 'not_assembled')"
    } else {
      "s.spectrum_type = 'assembled'"
    }
  } else {
    #if no spectrum_type column then no filtering
    ft_filter <- "1=1"
  }
  
  
  ms2_ft <- dbGetQuery(FT_con, sprintf("
    SELECT s.compound_id, p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p USING(spectrum_id)
    WHERE s.ms_level = 2
      AND s.polarity = %d
      AND %s
    ORDER BY s.compound_id, p.mz", polarity_ft, ft_filter))
  
  ms2_ft_list <- split(round(ms2_ft$mz), ms2_ft$compound_id)
  
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
  
  
  ms2_qtof <- dbGetQuery(QTOF_con, sprintf("
    SELECT s.compound_id, p.mz
    FROM msms_spectrum s
    JOIN msms_spectrum_peak p USING(spectrum_id)
    WHERE s.ms_level = 2
      AND s.polarity = %d
      AND %s
    ORDER BY s.compound_id, p.mz", polarity_qtof, qtof_filter))
  
  ms2_qtof_list <- split(round(ms2_qtof$mz), ms2_qtof$compound_id)
  
  # Loop over FT compounds
  for (i in seq_len(nrow(ft.exp))) {
    COMPID <- ft.exp$compound_id[i]
    x.tR   <- ft.exp$retention_time[i]
    if (is.na(x.tR) || x.tR < cutoff) next
    
    t1.tR <- ifelse(rg[3] != 0 & x.tR > rg[3], x.tR - rg[3], 0)
    t2.tR <- ifelse(rg[5] != 0 & x.tR > rg[5], x.tR - rg[5], 0)
    y.tR <- rg[1] + rg[2] * x.tR + rg[4] * t1.tR + rg[6] * t2.tR
    y.l  <- y.tR - lc.err
    y.h  <- y.tR + lc.err
    
    # Candidate matches
    pres <- Find_cand_matches_SQL(COMPID, err, syn.exp, ft.exp)
    if (nrow(pres) == 0) next
    
    ms2ion <- ms2_ft_list[[as.character(COMPID)]]
    if (is.null(ms2ion)) next
    
    pres1 <- data.frame()
    w.diff <- numeric()
    
    for (w in seq_len(nrow(pres))) {
      rt_w <- pres$retention_time[w]
      target_compid <- pres$compound_id[w]
      if (is.na(rt_w) || rt_w <= y.l || rt_w >= y.h) next
      
      msmsion <- ms2_qtof_list[[as.character(target_compid)]]
      if (is.null(msmsion)) next
      
      same_ion <- length(intersect(ms2ion, msmsion))
      if (same_ion / length(ms2ion) >= minIon) {
        pres1 <- rbind(pres1, pres[w, ])
        w.diff <- c(w.diff, abs(y.tR - rt_w))
      }
    }
    
    if (nrow(pres1) == 0) next
    pres1 <- pres1[order(w.diff), ]
    
    new_row <- data.frame(
      ref_compid      = COMPID,
      target_compid   = pres1$compound_id[1],
      ref_database    = basename(FT_path),
      target_database = basename(QTOF_path),
      stringsAsFactors = FALSE
    )
    
    # Only add if not already present
    is_new <- !any(
      Assoc$ref_compid      == new_row$ref_compid &
        Assoc$target_compid   == new_row$target_compid &
        Assoc$ref_database    == new_row$ref_database &
        Assoc$target_database == new_row$target_database
    )
    
    if (is_new) {
      Assoc <- rbind(Assoc, new_row)
    }
  }
  
  return(Assoc)
}
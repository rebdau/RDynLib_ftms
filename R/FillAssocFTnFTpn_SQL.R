#' @title Associate FT and QTOF features using SQL and MS2 matching
#'
#' @description
#' `FillAssocFTnFTpn_SQL()` finds associations between experiments from the same
#' instrument type, the one with negative polarity should be considered as
#' reference, based on retention time alignment, MS2 peak similarity, 
#' and alignment rules. It updates and returns the association table (`Assoc`).
#'
#' @details
#' The function:
#' - loads both databases compounds from SQL databases,
#' - computes expected RT alignment using regression parameters,
#' - finds candidate matches using SQL-assisted search,
#' - filters candidates using MS2 peak similarity,
#' - appends validated associations to `Assoc`.
#'
#' @param FTn_con `DBIConnection`  
#'   A DBI connection object to the reference SQLite database.
#'
#' @param FTp_con `DBIConnection`  
#'   A DBI connection object to the SQLite database to align with.
#'
#' @param Assoc `data.frame` A table storing previously detected associations. 
#'   Will be filled with the new association if they are not already 
#'   in the assoc table.
#'
#' @param FTn_expnr `numeric(1)` the reference experiment to align.
#'
#' @param FTp_expnr `numeric(1)` the experiment to align with the previous
#'        reference experiment.
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
#' @param polarity_ftn 
#'        `Integer (0/1)` Polarity of negative FTMS or QTOF experiment.
#'        
#' @param polarity_ftp 
#'        `Integer (0/1)` Polarity of positive FTMS or QTOF experiment.
#'
#' @param FTn_path `character(1)` the path to the reference database with 
#'        negative polarity.
#'        
#' @param FTp_path `character(1)` the path of the SQL database to align with.
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
FillAssocFTnFTpn_SQL <- function(
    FTn_con, FTp_con, Assoc, FTn_expnr, FTp_expnr, cutoff = 1, rg = rg,
    lc.err, err, minIon, polarity_ftn = 0, polarity_ftp = 1, FTn_path, FTp_path)
  {
  
  # Load FT and QTOF compounds
  ftn.exp <- dbGetQuery(FTn_con, sprintf(
    "SELECT retention_time, mass_measured, compound_id, expid 
     FROM ms_compound WHERE expid = %d", FTn_expnr))
  
  ftp.exp <- dbGetQuery(FTp_con, sprintf(
    "SELECT retention_time , mass_measured, compound_id, expid 
     FROM ms_compound WHERE expid = %d", FTp_expnr))
  
  o <- order(ftn.exp$mass_measured)
  ftng.o <- ftn.exp[o, ]
  o <- order(ftp.exp$mass_measured)
  ftps.o <- ftp.exp[o, ]
  
  i <- 1
  repeat {
    mass.ftng <- ftng.o$mass_measured[i]
    time.ftng <- ftng.o$retention_time[i]
    COMPID.ftng <- ftng.o$compound_id[i]
    
    pres <- Find_pos_compound(mass.ftng, time.ftng, ftps.o, err, lc.err, rg)
    
    if (nrow(pres) == 0) {
      if (i == nrow(ftng.o)) break
      i <- i + 1
      next
    } else {
      Assoc <- Sort_comp_matches(pres, Assoc, time.ftng, COMPID.ftng, rg)
    }
    
    if (i == nrow(ftng.o)) break
    i <- i + 1
  }
  
  return(Assoc)
}

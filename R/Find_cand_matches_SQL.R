#' @title Find Candidate Mass Matches Between FT and QTOF Compounds
#'
#'@description
#'
#' This function searches for candidate matching compounds in the QTOF dataset
#' (`syn.exp`) whose measured mass falls within a given mass tolerance (`err`)
#' around the FT compound mass (`ft.exp`).
#'
#' @param COMPID `integer(1)` The compound ID from the FT dataset for which
#'  candidate matches are requested. 
#' 
#' @param err `numeric(1)` The absolute mass tolerance window.
#' 
#' @param syn.exp `data.frame`  
#'   A data frame extracted from the QTOF SQLite database. Must contain columns:
#'   `"retention_time"`, `"mass_measured"`, `"compound_id"`, `"expid"`.
#'
#' @param ft.exp `data.frame`  
#'   A data frame extracted from the FT SQLite database. Must contain columns:
#'   `"retention_time"`, `"mass_measured"`, `"compound_id"`, `"expid"`.
#'
#' @return `data.frame`  
#'  A data frame containing QTOF compounds whose mass falls inside the tolerance
#'   window. Columns returned:
#'   - `retention_time` (`numeric`)
#'   - `mass_measured` (`numeric`)
#'   - `compound_id` (`integer`)
#'   - `expid` (`integer`)
#'
#' If no matches are found, an empty data.frame with the same columns is returned.
#' 
#'@noRd
Find_cand_matches_SQL <- function(COMPID, err, syn.exp, ft.exp) {
  # Find the FT compound
  i <- which(ft.exp$compound_id == COMPID)
  if (length(i) == 0) return(data.frame(
    retention_time = numeric(0),
    mass_measured  = numeric(0),
    compound_id    = integer(0),
    expid          = integer(0)
  ))
  
  ft_mass <- ft.exp$mass_measured[i]
  lb <- ft_mass - err
  ub <- ft_mass + err
  
  # Order syn.exp by mass_measured
  syn.o <- syn.exp[order(syn.exp$mass_measured), ]
  
  # Filter candidates within tolerance
  pres <- syn.o[syn.o$mass_measured > lb & syn.o$mass_measured < ub, ]
  
  # Ensure output is always a data.frame with the right columns
  if (nrow(pres) == 0) {
    pres <- data.frame(
      retention_time = numeric(0),
      mass_measured  = numeric(0),
      compound_id    = integer(0),
      expid          = integer(0)
    )
  }
  
  return(pres)
}

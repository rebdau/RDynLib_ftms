#' @title Aligning FTMS and QTOF experiments with the same/different polarities
#'
#' @description
#'
#' `Aligning_FT_QTOF_SQL()` aligns FTMS and QTOF experiments with the same
#' or different polarities, by calling a set of defined functions, and adds the
#' aligned compounds to a given text file.
#'
#' @details
#'
#' the function reads the both sql databases FTMS and QTOF, and extracts
#' their mode, create an output file to store the alignment results in a .txt
#' file if it's already given by the user otherwise we create a new one.
#' Then calls a set of defined functions :
#'
#' - [Aligning_General_SQL()]: performs the local alignment, and adjust the mz
#'   values if the alignment is between two different polarities.
#'
#' - `matchFTSyn_SQL()`: removes the matches found by the previous function
#'   if they contain fewer than `minIon` matched peaks.
#'
#' - `RemoveOutliers()`: calculate the standard deviation in the rt difference
#'   found between each matched feature and remove all rows outside a range.
#'
#' - `RegressionPie_LCalign()`: apply a regression model  to reduce the
#'   influence of false matches and outliers.
#'
#' - `FillAssocFTnQTOFn_SQL()`: store the final aligned compounds to a .txt
#'   file.
#'
#' @param FT_path `character(1)` the path to the FTMS SQL database.
#'
#' @param QTOF_path `character(1)` the path to the QTOF SQL database.
#'
#' @param FT_expnr `numeric(1)` the FTMS reference experiment to align.
#'
#' @param QTOF_expnr `numeric(1)` the QTOF experiment to align with the previous
#'        FTMS experiment.
#'
#' @param Assoc `character(1)` the file that contains some previous alignment
#'        results done by the user. If `NULL` we create it in the function.
#'
#' @param err `numeric(1)` is the m/z tolerance window used to decide which QTOF
#'        peaks are potential isomeric matches for a given FT peak.
#'
#' @param t.ini `numeric(1)` used in the `Aligning_General_SQL()` to define the
#'        number of neighboring peaks.
#'
#' @param lc.err `numeric(1)` the retention-time alignment tolerance.
#'
#' @param rng `numeric(1)` used in `RemoveOutliers()` refers to the number of
#'        standard  deviations used to define the outlier cutoff.
#'
#' @param minIon `numeric(1)`  
#'    The minimum matching ions between spectra of peak pairs, which are found
#'    after applying the regression model, and which are the final matching 
#'    written in the Assoc table. (default `0.01`).
#'
#' @param minPeaks Numeric. Minimum MS2 matching ratio required.
#'
#' @param startpoint `numeric(1)` used in `RegressionPie_LCalign()` refers to
#'        the minimum retention time for the first knot in piecewise regression.
#'
#' @param save_assoc `logical(1)` If true we add the alignment result to the
#'        Assoc text file.
#'
#' @param aggregated_Ft `logical(1)` If true the use of the aggregated ftms
#'        spectra of all levels is involved.
#'
#' @param aggregated_QTOF `logical(1)` If true the use of the aggregated QTOF
#'        spectra of all collision energies is involved.
#'
#' @return A `data.frame` containing the matched FT–QTOF compounds with columns:
#' - `ref_compid`: FTMS compound ID
#' - `target_compid`: QTOF compound ID
#' - `ref_database`: filename of the FTMS database
#' - `target_database`: filename of the QTOF database.
#' This table is updated with new associations found by the alignment procedure.
#'
#'
#' @importFrom DBI dbGetQuery
#'
#' @importFrom DBI dbConnect
#'
#' @importFrom DBI dbDisconnect
#'
#' @importFrom RSQLite SQLite
#'
#' @importFrom utils read.table
#'
#' @importFrom utils write.table
#'
#' @author Ahlam Mentag
#'
#' @export
Aligning_FT_QTOF_SQL <- function(
    FT_path, QTOF_path, FT_expnr, QTOF_expnr,
    Assoc = NULL, err = 0.02, t.ini = 5,
    lc.err = 1, rng = 2, minIon = 0.01, minPeaks = 0.6, startpoint = 1,
    save_assoc = FALSE, aggregated_Ft = FALSE, aggregated_QTOF = FALSE
) {
  
  if (is.character(Assoc)) {
    
    assoc_file <- Assoc
    
    if (file.exists(Assoc)) {
      Assoc_df <- read.table(Assoc, header = TRUE, sep = "\t")
    } else {
      warning("File ", assoc_file, " not found. Proceeding without.")
      Assoc_df <- data.frame(
        ref_compid = integer(0),
        target_compid = integer(0),
        ref_database = character(0),
        target_database = character(0),
        stringsAsFactors = FALSE
      )
    }
    
  } else if (is.data.frame(Assoc)) {
    
    Assoc_df <- Assoc
    assoc_file <- NULL
    
  } else {
    
    Assoc_df <- data.frame(
      ref_compid = integer(0),
      target_compid = integer(0),
      ref_database = character(0),
      target_database = character(0),
      stringsAsFactors = FALSE
    )
    assoc_file <- NULL
  }
  
  
  FT_con   <- dbConnect(SQLite(), FT_path)
  QTOF_con <- dbConnect(SQLite(), QTOF_path)
  on.exit(dbDisconnect(FT_con))
  on.exit(dbDisconnect(QTOF_con))
  
  # Detect polarity
  ft_mode  <- dbGetQuery(FT_con, sprintf("SELECT mode FROM experiment WHERE expid = %d", FT_expnr))$mode
  qtof_mode <- dbGetQuery(QTOF_con, sprintf("SELECT mode FROM experiment WHERE expid = %d", QTOF_expnr))$mode
  
  if (!length(ft_mode))
    stop("Experiment '", FT_expnr,"' not found in ", FT_path)
  if (!length(qtof_mode))
    stop("Experiment '", QTOF_expnr,"' not found in ", QTOF_path)
  
  ft_mode  <- tolower(ft_mode)
  qtof_mode <- tolower(qtof_mode)
  
  negative_values <- c("neg", "negative", "-", "negativ")
  
  polarity_ft  <- ifelse(ft_mode %in% negative_values, 0, 1)
  polarity_qtof <- ifelse(qtof_mode %in% negative_values, 0, 1)
  
  # Local LCal
  LCal <- Aligning_General_SQL(FT_con, QTOF_con, FT_expnr, QTOF_expnr, err, t.ini)
  
  LCal <- matchFTSyn_SQL(
    LCal, FT_con, QTOF_con, minPeaks = minPeaks,
    polarity_ft = polarity_ft, polarity_qtof = polarity_qtof,
    aggregated_Ft = aggregated_Ft, aggregated_QTOF = aggregated_QTOF
  )
  
  LCal <- RemoveOutliers(LCal, rng)
  
  # Regression
  rg <- RegressionPie_LCalign_SQL(LCal, startpoint)

  if (interactive() && capabilities("cairo")) {
    try(PlotPie_LCalign(LCal, rg), silent = TRUE)
  }
  
  # Fill assoc table
  Assoc_df <- FillAssocFTnQTOFn_SQL(
    FT_con = FT_con, QTOF_con = QTOF_con, Assoc = Assoc_df,
    FT_expnr = FT_expnr, QTOF_expnr = QTOF_expnr,
    cutoff = 1, rg = rg, lc.err = lc.err, err = err, minIon = minIon,
    polarity_ft = polarity_ft, polarity_qtof = polarity_qtof,
    FT_path = FT_path, QTOF_path = QTOF_path,
    aggregated_Ft = aggregated_Ft, aggregated_QTOF = aggregated_QTOF
  )
  
  # Save assoc if requested
  if (save_assoc && !is.null(assoc_file)) {
    write.table(Assoc_df, file = assoc_file, sep = "\t",
                row.names = FALSE, quote = FALSE)
  }
  
  return(Assoc_df)
}

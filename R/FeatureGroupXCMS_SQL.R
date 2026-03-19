#' @title Grouping of the features.
#'
#' @description
#'  the `FeatureGroupXCMS_SQL()` function performs the m/z feature grouping
#'   based on equal retention times and a high Pearson correlation across
#'   biological replicates.
#'  
#'
#' @param tR (e.g. 1 sec) defines the retention time difference between two 
#'       subsequent rows that is needed to define a breaking row
#'
#' @param dt the output of the featureDefinitions() of the xcmsExperiment 
#'        object ordred by retention time.
#'        
#' @param ints the intensities of the features in each sample.
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
FeatureGroupXCMS_SQL <- function(dt, ints, y10, y20, fc) {
  
  ## Order features by name
  xcms2i <- dt
  o <- order(xcms2i$name)
  xcms2 <- xcms2i[o, ]
  xo <- xcms2
  
  ## Ensure intensities are always a matrix
  if (is.null(dim(ints))) {
    ints <- matrix(
      ints,
      ncol = 1,
      dimnames = list(names(ints), "sample1")
    )
  }
  
  ## Reorder intensities using the same ordering as dt
  fc1 <- ints[o, , drop = FALSE]
  
  fc2 <- 1                # current sample group
  fc3 <- max(fc)          # number of sample groups
  fc4 <- 1                # first replicate column index
  
  repeat {
    
    fc5 <- length(fc[fc == fc2])   # number of replicates
    fc6 <- fc4 + fc5 - 1           
    
    x1 <- t(fc1[, fc4:fc6, drop = FALSE])
    k <- ncol(x1)
    nm <- as.character(xo$name)
    
    ## if fewer than 2 features then assign singleton groups
    if (k < 2) {
      grp_name <- paste0("nm.Grpo.Grp.", fc2)
      xcms2[[grp_name]] <- seq_along(nm)
      if (fc2 == fc3) break
      fc2 <- fc2 + 1
      fc4 <- fc4 + fc5
      next
    }
    
    colnames(x1) <- nm
    
    ## Compute correlations and RT distances
    z2 <- cor(x1, use = "p")
    z2 <- z2[lower.tri(z2)]
    
    z1 <- outer(xo$rtmed, xo$rtmed, "-")
    z1 <- abs(z1[lower.tri(z1)])
    
    jj <- (1:(k * (k - 1) / 2))[
      (z1 < y10) & !is.na(z1) & (z2 > y20) & !is.na(z2)
    ]
    
    ## Build index matrices
    rw <- cw <- integer(0)
    i <- 1
    while (i < k) {
      rw <- c(rw, seq(i + 1, k))
      cw <- c(cw, rep(i, k - i))
      i <- i + 1
    }
    
    res <- matrix(nrow = 0, ncol = 2)
    if (length(jj) > 0) {
      res <- cbind(rw[jj], cw[jj])
    }
    
    ## if no correlated pairs, assign singleton groups
    if (nrow(res) == 0) {
      grp_name <- paste0("nm.Grpo.Grp.", fc2)
      xcms2[[grp_name]] <- seq_along(nm)
      if (fc2 == fc3) break
      fc2 <- fc2 + 1
      fc4 <- fc4 + fc5
      next
    }
    
    ## Grouping correlated features
    grp <- 2
    Grp <- integer(length(nm))
    nm.Grp <- data.frame(nm = nm, Grp = Grp, stringsAsFactors = FALSE)
    
    repeat {
      set <- as.integer(res[1, ])
      l1 <- length(set)
      
      repeat {
        set <- unique(c(
          set,
          res[res[, 2] %in% set, 1],
          res[res[, 1] %in% set, 2]
        ))
        l2 <- length(set)
        if (l2 == l1) break
        l1 <- l2
      }
      
      nm.Grp$Grp[set] <- grp
      res <- res[!(res[, 1] %in% set), , drop = FALSE]
      
      if (nrow(res) < 1) break
      grp <- grp + 1
    }
    
    o <- order(nm.Grp$nm)
    grp_name <- paste0("nm.Grpo.Grp.", fc2)
    xcms2[[grp_name]] <- nm.Grp$Grp[o]
    
    if (fc2 == fc3) break
    fc2 <- fc2 + 1
    fc4 <- fc4 + fc5
  }
  
  ## Final consensus grouping across all sample groups
  nrcr <- ncol(xcms2) - fc3
  ii <- do.call(order, data.frame(xcms2[, -(1:nrcr), drop = FALSE]))
  xcms2 <- xcms2[ii, ]
  
  xcms2.cr <- xcms2[, -(1:nrcr), drop = FALSE]
  con <- rep(1L, nrow(xcms2.cr))
  
  i <- 2
  j <- 1
  l <- ncol(xcms2.cr)
  
  repeat {
    fin <- logical(l)
    for (k in seq_len(l)) {
      fin[k] <- (xcms2.cr[i, k] == xcms2.cr[j, k]) ||
        (xcms2.cr[j, k] == 1)
    }
    
    con[i] <- if (all(fin)) con[j] else con[j] + 1
    
    if (i == nrow(xcms2.cr)) break
    i <- i + 1
    j <- j + 1
  }
  
  xcms2 <- cbind(xcms2, con)
  return(xcms2)
}

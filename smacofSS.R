source("smacofAuxiliaries.R")
source("smacofDataUtilities.R")
source("smacofPlots.R")

smacofSS <- function(theData,
                     ndim = 2,
                     xinit = NULL,
                     weighted = FALSE,
                     ordinal = FALSE,
                     itmax = 1000,
                     eps = 1e-10,
                     digits = 10, 
                     width = 12,
                     verbose = FALSE) {
  if (!weighted && !ordinal) {
    source("smacofSSUR.R")
    return(
      smacofSSUR(
        theData = theData,
        ndim = ndim,
        xinit = xinit,
        itmax = itmax,
        eps = eps,
        digits = digits, 
        width = width,
        verbose = verbose
      )
    )
  }
  if (weighted && !ordinal) {
    source("smacofSSWR.R")
    return(
      smacofSSWR(
        theData = theData,
        ndim = ndim,
        xinit = xinit,
        itmax = itmax,
        eps = eps,
        digits = digits, 
        width = width,
        verbose = verbose
      )
    )
  }
  if (!weighted && ordinal) {
    source("smacofSSUO.R")
    return(
      smacofSSUO(
        theData = theData,
        ndim = ndim,
        xinit = xinit,
        ties = ordinal,
        itmax = itmax,
        eps = eps,
        digits = digits, 
        width = width,
        verbose = verbose
      )
    )
  }
  if (weighted && ordinal) {
    source("smacofSSWO.R")
    return(
      smacofSSWO(
        theData = theData,
        ndim = ndim,
        xinit = xinit,
        ties = ordinal,
        itmax = itmax,
        eps = eps,
        digits = digits, 
        width = width,
        verbose = verbose
      )
    )
  }
}


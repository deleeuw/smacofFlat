dyn.load("monotone.so")
dyn.load("smacofStep.so")
source("smacofUtils.R")
source("smacofSSUR.R")
source("smacofSSWR.R")
source("smacofSSUO.R")
source("smacofSSWO.R")

smacofSS <- function(theData,
                     ndim,
                     xold,
                     weighted,
                     ordinal,
                     ties,
                     itmax,
                     eps,
                     verbose) {
  if (!weighted && !ordinal) {
    return(
      smacofSSUR(
        theData = theData,
        ndim = ndim,
        xold = xold,
        itmax = itmax,
        eps = eps,
        verbose = verbose
      )
    )
  }
  if (weighted && !ordinal) {
    return(
      smacofSSWR(
        theData = theData,
        ndim = ndim,
        xold = xold,
        itmax = itmax,
        eps = eps,
        verbose = verbose
      )
    )
  }
  if (!weighted && ordinal) {
    return(
      smacofSSUO(
        theData = theData,
        ndim = ndim,
        xold = xold,
        ties = ties,
        itmax = itmax,
        eps = eps,
        verbose = verbose
      )
    )
  }
  if (weighted && ordinal) {
    return(
      smacofSSWO(
        theData = theData,
        ndim = ndim,
        xold = xold,
        ties = ties,
        itmax = itmax,
        eps = eps,
        verbose = verbose
      )
    )
  }
}





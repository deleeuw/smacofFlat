dyn.load("smacofSS.so")

source("smacofAuxiliaries.R")
source("smacofDataUtilities.R")
source("smacofPlots.R")
source("smacofTorgerson.R")

smacofSS <- function(theData,
                     ndim = 2,
                     xinit = NULL,
                     ties = 1,
                     itmax = 100,
                     eps = 1e-10,
                     digits = 10,
                     width = 15,
                     verbose = 1,
                     weighted = 1,
                     ordinal = 1) {
  if (is.null(xinit)) {
    xinit <- smacofTorgerson(theData, 2)$conf
  }
  xold <- xinit
  nobj <- theData$nobj
  ndat <- theData$ndat
  itel <- 1
  iind <- theData$iind
  jind <- theData$jind
  dhat <- theData$delta
  wght <- theData$weights
  blks <- theData$blocks
  edis <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    edis[k] <- sqrt(sum((xold[i, ] - xold[j, ])^2))
  }
  dhat <- dhat * sqrt(ndat / sum(dhat^2))
  sdd <- sum(edis^2)
  sde <- sum(dhat * edis)
  lbd <- sde / sdd
  edis <- lbd * edis
  xold <- lbd * xold
  sold <- sum((dhat - edis)^2) / ndat
  snew <- 0.0
  xold <- as.vector(xold)
  xnew <- xold
  h <- .C(
    "smacofSSEngine",
    nobj = as.integer(nobj),
    ndim = as.integer(ndim),
    ndat = as.integer(ndat),
    itel = as.integer(itel),
    ties = as.integer(ties),
    itmax = as.integer(itmax),
    digits = as.integer(digits),
    width = as.integer(width),
    verbose = as.integer(verbose),
    ordinal = as.integer(ordinal),
    weighted = as.integer(weighted),
    sold = as.double(sold),
    snew = as.double(snew),
    eps = as.double(eps),
    iind = as.integer(iind - 1),
    jind = as.integer(jind - 1),
    blks = as.integer(blks),
    wght = as.double(wght),
    edis = as.double(edis),
    dhat = as.double(dhat),
    xold = as.double(xold),
    xnew = as.double(xnew)
  )
  result <- list(
    delta = theData$delta,
    dhat = h$dhat,
    confdist = h$edis,
    conf = matrix(h$xnew, nobj, ndim),
    weightmat = theData$weights,
    stress = h$snew,
    ndim = ndim,
    init = xinit,
    niter = h$itel,
    nobj = nobj,
    iind = h$iind,
    jind = h$jind,
    weighted = weighted,
    ordinal = ordinal
  )
  class(result) <- c("smacofSSResult", "smacofSSUOResult")
  return(result)
}

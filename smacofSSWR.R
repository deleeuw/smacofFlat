dyn.load("smacofSSWREngine.so")

source("smacofAuxiliaries.R")
source("smacofDataUtilities.R")
source("smacofPlots.R")

smacofSSWR <- function(theData,
                       ndim = 2,
                       xinit = NULL,
                       itmax = 1000,
                       eps = 1e-10,
                       digits = 10,
                       width = 15,
                       verbose = TRUE) {
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
  wsum <- sum(wght)
  vinv <- makeMPInverseV(theData)
  edis <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    edis[k] <- sqrt(sum((xold[i, ] - xold[j, ])^2))
  }
  dhat <- dhat * sqrt(wsum / sum(wght * dhat^2))
  sdd <- sum(wght * edis^2)
  sde <- sum(wght * dhat * edis)
  lbd <- sde / sdd
  edis <- lbd * edis
  xold <- lbd * xold
  sold <- sum(wght * (dhat - edis)^2) / wsum
  snew <- 0.0
  xold <- as.vector(xold)
  xnew <- rep(0, nobj * ndim)
  h <- .C(
    "smacofSSWREngine",
    nobj = as.integer(nobj),
    ndim = as.integer(ndim),
    ndat = as.integer(ndat),
    itel = as.integer(itel),
    itmax = as.integer(itmax),
    digits = as.integer(digits),
    width = as.integer(width),
    verbose = as.integer(verbose),
    wsum = as.double(wsum),
    sold = as.double(sold),
    snew = as.double(snew),
    eps = as.double(eps),
    iind = as.integer(iind),
    jind = as.integer(jind),
    edis = as.double(edis),
    dhat = as.double(dhat),
    wght = as.double(wght),
    vinv = as.double(vinv),
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
    weighted = TRUE,
    ordinal = FALSE
  )
  class(result) <- c("smacofSSResult", "smacofSSWRResult")
  return(result)
}

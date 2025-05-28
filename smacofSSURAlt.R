dyn.load("smacofSSUREngine.so")

source("smacofUtils.R")
source("smacofSSData.R")

smacofSSURAlt <- function(theData,
                          ndim = 2,
                          xold = torgerson(theData, ndim),
                          itmax = 1000,
                          eps = 1e-10,
                          verbose = TRUE) {
  nobj <- nrow(xold)
  ndat <- nrow(theData)
  itel <- 1
  iind <- theData[, 1]
  jind <- theData[, 2]
  dhat <- theData[, 3]
  edis <- rep(0, ndat)
  h <- .C(
    "smacofDistances",
    nobj = as.integer(nobj),
    ndim = as.integer(ndim),
    ndat = as.integer(ndat),
    iind = as.integer(iind),
    jind = as.integer(jind),
    edis = as.double(edis),
    xold = as.double(xold)
  )
  edis <- h$edis
  dhat <- dhat / sqrt(sum(dhat^2))
  sdd <- sum(edis^2)
  sde <- sum(dhat * edis)
  lbd <- sde / sdd
  edis <- lbd * edis
  xold <- lbd * xold
  sold <- sum((dhat - edis)^2)
  snew <- 0.0
  xnew <- rep(0, nobj * ndim)
  h <- .C(
    "smacofSSUREngine",
    nobj = as.integer(nobj),
    ndim = as.integer(ndim),
    ndat = as.integer(ndat),
    itel = as.integer(itel),
    itmax = as.integer(itmax),
    verbose = as.integer(verbose),
    sold = as.double(sold),
    snew = as.double(snew),
    eps = as.double(eps),
    iind = as.double(iind),
    jind = as.double(jind),
    edis = as.double(edis),
    dhat = as.double(dhat),
    xold = as.double(xold),
    xnew = as.double(xnew)
  )
  return(
    list(
      data = theData,
      conf = h$xnew,
      loss = h$snew,
      edis = h$edis,
      dhat = h$dhat,
      itel = h$itel
    )
  )
}

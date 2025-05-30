dyn.load("smacofSSUREngine.so")

source("../smacofUtils.R")
source("../smacofSSData.R")

smacofSSUR <- function(theData,
                       ndim = 2,
                       xold = torgerson(theData, ndim),
                       itmax = 1000,
                       eps = 1e-10,
                       digits = 10,
                       width = 15,
                       verbose = TRUE) {
  nobj <- nrow(xold)
  ndat <- nrow(theData)
  itel <- 1
  iind <- theData[, 1]
  jind <- theData[, 2]
  dhat <- theData[, 3]
  edis <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    edis[k] <- sqrt(sum((xold[i, ] - xold[j, ])^2))
  }
  dhat <- dhat / sqrt(sum(dhat^2))
  sdd <- sum(edis^2)
  sde <- sum(dhat * edis)
  lbd <- sde / sdd
  edis <- lbd * edis
  xold <- lbd * xold
  sold <- sum((dhat - edis)^2)
  snew <- 0.0
  xold <- as.vector(xold)
  xnew <- rep(0, nobj * ndim)
  h <- .C(
    "smacofSSUREngine",
    nobj = as.integer(nobj),
    ndim = as.integer(ndim),
    ndat = as.integer(ndat),
    itel = as.integer(itel),
    itmax = as.integer(itmax),
    digits = as.integer(digits),
    width = as.integer(width),
    verbose = as.integer(verbose),
    sold = as.double(sold),
    snew = as.double(snew),
    eps = as.double(eps),
    iind = as.integer(iind),
    jind = as.integer(jind),
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

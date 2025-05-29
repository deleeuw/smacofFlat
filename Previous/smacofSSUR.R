dyn.load("smacofStep.so")

source("../smacofUtils.R")
source("../smacofSSData.R")

smacofSSUR <- function(theData,
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
  repeat {
    xnew <- matrix(0, nobj, ndim)
    h <- .C(
      "smacofStepUnweighted",
      nobj = as.integer(nobj),
      ndim = as.integer(ndim),
      ndat = as.integer(ndat),
      iind = as.integer(iind),
      jind = as.integer(jind),
      edis = as.double(edis),
      dhat = as.double(dhat),
      xold = as.double(as.vector(xold)),
      xnew = as.double(as.vector(xnew))
    )
    xnew <- matrix(h$xnew, nobj, ndim)
    edis <- h$edis
    snew <- sum((dhat - edis)^2)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, digits = 3, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    sold <- snew
    xold <- xnew
    itel <- itel + 1
  }
  return(list(
    data = theData,
    conf = xnew,
    loss = snew,
    edis = edis,
    dhat = dhat,
    itel = itel
  ))
}

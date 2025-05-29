smacofSSWR <- function(theData,
                       ndim,
                       xold,
                       itmax,
                       eps,
                       verbose) {
  nobj <- nrow(xold)
  ndat <- nrow(theData)
  itel <- 1
  iind <- theData[, 1]
  jind <- theData[, 2]
  dhat <- theData[, 3]
  wght <- theData[, 5]
  vinv <- makeMPInverseV(theData)
  itel <- 1
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
  dhat <- dhat / sqrt(sum(wght * dhat^2))
  sdd <- sum(wght * edis^2)
  sde <- sum(wght * dhat * edis)
  lbd <- sde / sdd
  edis <- lbd * edis
  xold <- lbd * xold
  sold <- sum(wght * (dhat - edis)^2)
  repeat {
    xtmp <- matrix(0, nobj, ndim)
    xnew <- matrix(0, nobj, ndim)
    h <- .C(
      "smacofStepWeighted",
      nobj = as.integer(nobj),
      ndim = as.integer(ndim),
      ndat = as.integer(ndat),
      iind = as.integer(iind),
      jind = as.integer(jind),
      edis = as.double(edis),
      dhat = as.double(dhat),
      wght = as.double(wght),
      vinv = as.double(vinv),
      xold = as.double(as.vector(xold)),
      xtmp = as.double(as.vector(xtmp)),
      xnew = as.double(as.vector(xnew))
    )
    xnew <- matrix(h$xnew, nobj, ndim)
    edis <- h$edis
    snew <- sum(wght * (dhat - edis)^2)
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
    dist = edis,
    dhat = dhat,
    itel = itel
  ))
}

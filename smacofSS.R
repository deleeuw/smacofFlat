dyn.load("monotone.so")
dyn.load("smacofStep.so")
source("smacofUtils.R")

small <- as.dist(matrix(c(0, 1, 3, 3, 1, 0, 1, 3, 3, 1, 0, 1, 3, 3, 1, 0), 4, 4))
smallData <- makeMDSData(small)

data(ekman, package = "smacof")
ekmanData <- makeMDSData((1 - ekman)^3)

data(wish, package = "smacof")
wishData <- makeMDSData(7 - wish, 1 / (7 - wish))

data(morse, package = "smacof")
morseData <- makeMDSData(morse)



smacofSS <- function(theData,
                     ndim = 2,
                     xold = torgerson(theData, ndim),
                     weighted = FALSE,
                     ordinal = FALSE,
                     ties = 2,
                     itmax = 1000,
                     eps = 1e-10,
                     verbose = TRUE) {
  if (!weighted && !ordinal) {
    return(
      smacofSSUR(
        theData,
        ndim = ndim,
        xold = torgerson(theData, ndim),
        itmax = itmax,
        eps = eps,
        verbose = verbose
      )
    )
  }
  if (weighted && !ordinal) {
    return(
      smacofSSWR(
        theData,
        ndim = ndim,
        xold = torgerson(theData, ndim),
        itmax = itmax,
        eps = eps,
        verbose = verbose
      )
    )
  }
  if (!weighted && ordinal) {
    return(
      smacofSSUO(
        theData,
        ndim = ndim,
        xold = torgerson(theData, ndim),
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
        theData,
        ndim = ndim,
        xold = torgerson(theData, ndim),
        ties = ties,
        itmax = itmax,
        eps = eps,
        verbose = verbose
      )
    )
  }
}

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

smacofSSWR <- function(theData,
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

smacofSSUO <- function(theData,
                       ndim = 2,
                       xold = torgerson(theData, ndim),
                       ties = 2,
                       itmax = 1000,
                       eps = 1e-10,
                       verbose = TRUE) {
  nobj <- nrow(xold)
  ndat <- nrow(theData)
  itel <- 1
  iind <- theData[, 1]
  jind <- theData[, 2]
  dhat <- theData[, 3]
  blks <- theData[, 4]
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
    smid <- sum((dhat - edis)^2)
    if (ties ==  1) {
      h <- .C(
        "primaryApproach",
        ndat = as.integer(ndat),
        blks = as.integer(theData[, 4]),
        dhat = as.double(edis),
        wvec = as.double(rep(1, ndat)),
        indi = 1:ndat
      )
      theData[, 1] <- theData[h$indi, 1]
      theData[, 2] <- theData[h$indi, 2]
      edis <- edis[h$indi]
      iind <- iind[h$indi]
      jind <- jind[h$indi]
    }
    if (ties == 2) {
      h <- .C(
        "secondaryApproach",
        ndat = as.integer(ndat),
        blks = as.integer(theData[, 4]),
        dhat = as.double(edis),
        wvec = as.double(rep(1, ndat))
      )
      
    }
    if (ties == 3) {
      h <- .C(
        "tertiaryApproach",
        ndat = as.integer(ndat),
        blks = as.integer(theData[, 4]),
        dhat = as.double(edis),
        wvec = as.double(rep(1, ndat))
      )
      
    }
    dhat <- h$dhat
    dhat <- dhat / sqrt(sum(dhat^2))
    snew <- sum((dhat - edis)^2)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, digits = 3, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "smid ",
        formatC(smid, digits = 10, format = "f"),
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

smacofSSWO <- function(theData,
                       ndim = 2,
                       xold = torgerson(theData, ndim),
                       ties = 2,
                       itmax = 1000,
                       eps = 1e-10,
                       verbose = TRUE) {
  nobj <- nrow(xold)
  ndat <- nrow(theData)
  itel <- 1
  iind <- theData[, 1]
  jind <- theData[, 2]
  dhat <- theData[, 3]
  blks <- theData[, 4]
  wght <- theData[, 5]
  vinv <- makeMPInverseV(theData)
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
    smid <- sum(wght * (dhat - edis)^2)
    if (ties ==  1) {
      h <- .C(
        "primaryApproach",
        ndat = as.integer(ndat),
        blks = as.integer(theData[, 4]),
        dhat = as.double(edis),
        wvec = as.double(rep(1, ndat)),
        indi = 1:ndat
      )
      theData[, 1] <- theData[h$indi, 1]
      theData[, 2] <- theData[h$indi, 2]
      edis <- edis[h$indi]
      iind <- iind[h$indi]
      jind <- jind[h$indi]
    }
    if (ties == 2) {
      h <- .C(
        "secondaryApproach",
        ndat = as.integer(ndat),
        blks = as.integer(theData[, 4]),
        dhat = as.double(edis),
        wvec = as.double(rep(1, ndat))
      )
      
    }
    if (ties == 3) {
      h <- .C(
        "tertiaryApproach",
        ndat = as.integer(ndat),
        blks = as.integer(theData[, 4]),
        dhat = as.double(edis),
        wvec = as.double(rep(1, ndat))
      )
      
    }
    dhat <- h$dhat
    dhat <- dhat / sqrt(sum(wght * dhat^2))
    snew <- sum(wght * (dhat - edis)^2)
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, digits = 3, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "smid ",
        formatC(smid, digits = 10, format = "f"),
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

library(RSpectra)
makeMPInverseV <- function(theData) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  wght <- matrix(0, nobj, nobj)
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    wght[i, j] <- wght[j, i] <- theData$weights[k]
  }
  vmat <- -wght
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  return(as.vector(as.dist(vinv)))
}

smacofTorgerson <- function(theData, ndim) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  dmat <- matrix(0, nobj, nobj)
  mdel <- mean(theData$delta)^2
  dmat <- mdel * (1 - diag(nobj))
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    dmat[i, j] <- dmat[j, i] <- theData$delta[k]
  }
  dmat <- dmat^2
  dr <- apply(dmat, 1, mean)
  dm <- mean(dmat)
  bmat <- -(dmat - outer(dr, dr, "+") + dm) / 2
  ev <- eigs_sym(bmat, ndim, which = "LA")
  evev <- diag(sqrt(pmax(0, ev$values)))
  x <- ev$vectors[, 1:ndim] %*% evev
  result <- list(
    delta = theData$delta,
    dhat = dhat,
    confdist = edis,
    conf = x,
    weightmat = theData$weights,
    stress = sstress,
    ndim = ndim,
    nobj = nobj,
    iind = theData$iind,
    jind = theData$jind
  )
  class(result) <- c("smacofSSResult", "smacofGuttmanResult")
  return(result)}

smacofGuttman <- function(theData, ndim) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  delta <- theData$delta
  wght <- theData$weights
  wsum <- sum(wght)
  dhat <- (theData$delta)^2
  dhat <- dhat * sqrt(wsum / sum(wght * dhat^2))
  bmat <- matrix(0, nobj, nobj)
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    bmat[i, j] <- bmat[j, i] <- wght[k] * dhat[k]
  }
  bmat <- -bmat
  diag(bmat) <- -rowSums(bmat)
  ev <- eigs_sym(bmat, ndim, which = "LA")
  evev <- diag(sqrt(pmax(0, ev$values)))
  x <- ev$vectors[, 1:ndim] %*% evev
  edis <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    edis[k] <- sum((x[i,] - x[j, ])^2)
  }
  sstress <- sum(wght * (dhat - edis)^2)
  result <- list(
    delta = theData$delta,
    dhat = dhat,
    confdist = edis,
    conf = x,
    weightmat = theData$weights,
    stress = sstress,
    ndim = ndim,
    nobj = nobj,
    iind = theData$iind,
    jind = theData$jind
  )
  class(result) <- c("smacofSSResult", "smacofGuttmanResult")
  return(result)
}

smacofRandomConfiguration <- function(theData, ndim = 2) {
  nobj <- theData$nobj
  x <- matrix(rnorm(nobj * ndim), nobj, ndim)
  return(columnCenter(x))
}

columnCenter <- function(x) {
  apply(x, 2, function(x) x - mean(x))
}

matrixPrint <- function(x,
                        digits = 6,
                        width = 8,
                        format = "f",
                        flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}

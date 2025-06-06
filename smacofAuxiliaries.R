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

torgerson <- function(theData, ndim) {
  nobj <- theData$nobj
  dmat <- matrix(0, nobj, nobj)
  ndat <- length(theData$iind)
  mdel <- mean(theData$delta)^2
  dmat <- mdel * (1 - diag(nobj))
  for (k in 1:ndat) {
    i <- theData$iind
    j <- theData$jind
    dmat[i, j] <- dmat[j, i] <- theData$delta[k]
  }
  dmat <- dmat^2
  dr <- apply(dmat, 1, mean)
  dm <- mean(dmat)
  cmat <- -(dmat - outer(dr, dr, "+") + dm) / 2
  ev <- eigs_sym(cmat, ndim, which = "LA")
  x <- ev$vectors[, 1:ndim] %*% diag(sqrt(drop(ev$values[1:ndim])))
  return(x)
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

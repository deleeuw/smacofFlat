library(RSpectra)

makeMPInverseV <- function(theData) {
  nobj <- max(theData[, 1:2])
  ndat <- nrow(theData)
  wght <- matrix(0, nobj, nobj)
  for (k in 1:ndat) {
    i <- theData[k, 1]
    j <- theData[k, 2]
    wght[i, j] <- wght[j, i] <- theData[k, 5]
  }
  vmat <- -wght
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  return(as.vector(as.dist(vinv)))
}

torgerson <- function(theData, ndim) {
  nobj <- max(theData[, 1])
  dmat <- matrix(0, nobj, nobj)
  ndat <- nrow(theData)
  mdel <- mean(theData[, 3])^2
  dmat <- mdel * (1 - diag(nobj))
  for (k in 1:ndat) {
    i <- theData[k, 1]
    j <- theData[k, 2]
    dmat[i, j] <- dmat[j, i] <- theData[k, 3]
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

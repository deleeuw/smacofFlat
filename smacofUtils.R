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

makeMDSData <- function(delta, weights = NULL) {
  nobj <- attr(delta, "Size")
  ndat <- nobj * (nobj - 1) / 2
  theData <- NULL
  k <- 1
  for (j in 1:(nobj - 1)) {
    for (i in (j + 1):nobj) {
      if (is.null(weights)) {
        theData <- rbind(theData, c(i, j, delta[k], 0))
      } else {
        theData <- rbind(theData, c(i, j, delta[k], 0, weights[k]))
      }
      k <- k + 1
    }
  }
  if (is.null(weights)) {
    colnames(theData) <- c("i", "j", "delta", "blocks")
  } else {
    colnames(theData) <- c("i", "j", "delta", "blocks", "weights")
  }
  delta <- theData[, 3]
  theData <- theData[order(delta), ]
  theData[, 4] <- makeTieBlocks(theData[, 3])
  return(theData)
}

makeTieBlocks <- function(x) {
  n <- length(x)
  y <- rep(0, n)
  k <- 1
  repeat {
    m <- length(which(x == x[k]))
    y[k] <- m
    k <- k + m
    if (k > n) {
      break
    }
  }
  return(y)
}

makeMPInverseV <- function(theData) {
  nobj <- max(theData[, 1])
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
    dmat[i, j] <- dmat[j, i] <- theData[k, 3]^2
  }
  dr <- apply(dmat, 1, mean)
  dm <- mean(dmat)
  cmat <- -(dmat - outer(dr, dr, "+") + dm) / 2
  ev <- eigen(cmat)
  x <- ev$vectors[, 1:ndim] %*% diag(sqrt(ev$values[1:ndim]), nrow = ndim, ncol = ndim)
  return(x)
}

fromMDSData <- function(theData) {
  ndat <- nrow(theData)
  nobj <- (1 + sqrt(1 + 8 * ndat)) / 2
  delta <- matrix(0, nobj, nobj)
  weights <- matrix(0, nobj, nobj)
  for (k in 1:ndat) {
    i <- theData[k, 1]
    j <- theData[k, 2]
    delta[i, j] <- delta[j, i] <- theData[k, 3]
    weights[i, j] <- weights[j, i] <- theData[k, 5]
  }
  return(list(delta = as.dist(delta), weights = as.dist(weights)))
}


library("numDeriv")

source("smacofDataUtilities.R")
source("ekmanData.R")

smacofDerivatives <- function(h) {
  iind <- h$iind + 1
  jind <- h$jind + 1
  dhat <- h$dhat
  wght <- h$weightmat
  ndat <- length(dhat)
  nobj <- h$nobj
  ndim <- h$ndim
  x <- h$conf
  v <- matrix(0, nobj, nobj)
  b <- matrix(0, nobj, nobj)
  d <- as.matrix(dist(x))
  gg <- matrix(0, nobj * ndim, nobj * ndim)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    v[i, j] <- -wght[k]
    b[i, j] <- -wght[k] * dhat[k] / d[i, j]
  }
  v <- v + t(v)
  b <- b + t(b)
  diag(v) <- -rowSums(v)
  diag(b) <- -rowSums(b)
  grad <- (v - b) %*% x
  for (s in 1:ndim) {
    is <- (s - 1) * nobj
    for (t in 1:ndim) {
      js <- (t - 1) * nobj
      for (k in 1:ndat) {
        i <- iind[k]
        j <- jind[k]
        aux <- (x[i, s] - x[j, s]) * (x[i, t] - x[j, t])
        gg[is + i, js + j] <- -(wght[k] * dhat[k] * aux) / d[i, j]^3
      }
    gg[is + 1:nobj, js + 1:nobj] <- gg[is + 1:nobj, js + 1:nobj] + t(gg[is + 1:nobj, js + 1:nobj])
    diag(gg[is + 1:nobj, js + 1:nobj]) <- -rowSums(gg[is + 1:nobj, js + 1:nobj])
    }
  }
  hh <- -gg
  dd <- gg
  for (s in 1:ndim) {
    is <- (s - 1) * nobj
    hh[is + 1:nobj, is + 1:nobj] <- hh[is + 1:nobj, is + 1:nobj] + b
    dd[is + 1:nobj, is + 1:nobj] <- dd[is + 1:nobj, is + 1:nobj] + v - b
  }
  return(list(v = v, b = b, grad = grad, gg = gg, hh = hh, hess = dd))
}

theFunc <- function(x) {
  xmat <- matrix(x, nobj, ndim)
  s <- 0.0
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    d <- sqrt(sum((xmat[i, ] - xmat[j, ])^2))
    s <- s + wght[k] * (delta[k] - d)^2
  }
  return(s / 2)
}

numDifDist <- function(x) {
  x <- as.vector(x)
  g <- grad(theFunc, x)
  h <- hessian(theFunc, x)
  return(list(grad = g, hess = h))
}

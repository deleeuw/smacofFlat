smacofSSM <- function(delta,
                      weights = NULL,
                      ndim = 2,
                      xold = torgersonM(delta, ndim),
                      ordinal = FALSE,
                      itmax = 1000,
                      eps = 1e-10,
                      verbose = TRUE) {
  delta <- as.matrix(delta)
  nobj <- nrow(delta)
  if (!is.null(weights)) {
    weights <- as.matrix(weights)
    vmat <- -weights
    diag(vmat) <- -rowSums(vmat)
    vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  }
  itel <- 1
  dold <- as.matrix(dist(xold))
  if (is.null(weights)) {
    dhat <- delta / sqrt(sum(delta^2))
    sdd <- sum(dold^2)
    sde <- sum(dhat * dold)
    lbd <- sde / sdd
    dold <- lbd * dold
    xold <- lbd * xold
    sold <- sum((dhat - dold)^2)
  } else {
    dhat <- delta / sqrt(sum(weights * delta^2))
    sdd <- sum(weights * dold^2)
    sde <- sum(weights * dhat * dold)
    lbd <- sde / sdd
    dold <- lbd * dold
    xold <- lbd * xold
    sold <- sum(weights * (dhat - dold)^2)
  }
  repeat {
    if (!is.null(weights)) {
      bmat <- -weights * dhat / (dold + diag(nobj))
      diag(bmat) <- -rowSums(bmat)
      xnew <- vinv %*% bmat %*% xold
    } else {
      bmat <- -dhat / (dold + diag(nobj))
      diag(bmat) <- -rowSums(bmat)
      xnew <- (bmat %*% xold) / nobj
    }
    dnew <- as.matrix(dist(xnew))
    if (is.null(weights)) {
      snew <- sum((dhat - dnew)^2)
    } else {
      snew <- sum(weights * (dhat - dnew)^2)
    }
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
    dold <- dnew
    itel <- itel + 1
  }
}

torgersonM <- function(delta, ndim) {
  dmat <- delta^2
  dr <- apply(dmat, 2, mean)
  dm <- mean(dmat)
  cmat <- -(dmat - outer(dr, dr, "+") + dm) / 2
  ev <- eigen(cmat)
  x <- ev$vectors[, 1:ndim] %*% diag(sqrt(ev$values[1:ndim]), nrow = ndim, ncol = ndim)
  return(x)
}
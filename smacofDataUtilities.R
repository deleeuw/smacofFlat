
makeMDSData <- function(delta, weights = NULL) {
  nobj <- attr(delta, "Size")
  if (is.null(weights)) {
    weights <- as.dist(1 - diag(nobj))
  }
  theData <- NULL
  k <- 1
  for (j in 1:(nobj - 1)) {
    for (i in (j + 1):nobj) {
      if ((weights[k] > 0) && (!is.na(weights[k])) && (!is.na(delta[k]))) {
        theData <- rbind(theData, c(i, j, delta[k], 0, weights[k]))
      }
      k <- k + 1
    }
  }
  colnames(theData) <- c("i", "j", "delta", "blocks", "weights")
  ndat <- nrow(theData)
  theData <- theData[order(theData[, 3]), ]
  dvec <- theData[, 3]
  k <- 1
  repeat {
    m <- length(which(dvec == dvec[k]))
    theData[k, 4] <- m
    k <- k + m
    if (k > ndat) {
      break
    }
  }
  return(theData)
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


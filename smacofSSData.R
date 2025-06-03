small <- as.dist(matrix(c(0, 1, 3, 2, 1, 0, 1, 3, 3, 1, 0, 1, 2, 3, 1, 0), 4, 4))
smallData <- makeMDSData(small)
smallLabels <- c("A", "B", "C", "D")

data(ekman, package = "smacof")
ekmanData <- makeMDSData((1 - ekman)^3)
ekmanLabels <- as.character(attr(ekman, "Labels"))

data(wish, package = "smacof")
wishData <- makeMDSData(7 - wish)
wishLabels <- as.character(attr(wish, "Labels"))

data(morse, package = "smacof")
morseData <- makeMDSData(morse, 1 / morse)
morseLabels <- as.character(attr(morse, "Labels"))


# data(iris)
# irisDist <- dist(iris[,1:4])
# irisData <- makeMDSData(irisDist)


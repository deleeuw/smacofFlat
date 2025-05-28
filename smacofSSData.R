small <- as.dist(matrix(c(0, 1, 3, 2, 1, 0, 1, 3, 3, 1, 0, 1, 2, 3, 1, 0), 4, 4))
smallData <- makeMDSData(small)

data(ekman, package = "smacof")
ekmanData <- makeMDSData((1 - ekman)^3)

data(wish, package = "smacof")
wishData <- makeMDSData(7 - wish, 1 / (7 - wish))

data(morse, package = "smacof")
morseData <- makeMDSData(morse, 1 / morse)


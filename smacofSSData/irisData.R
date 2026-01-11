source("smacofDataUtilities.R")

data(iris)
irisDist <- dist(iris[,1:4])
irisMatrix <- as.matrix(irisDist)
irisData <- makeMDSData(irisDist)


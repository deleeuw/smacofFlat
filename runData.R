source("smacofDataUtilities.R")

vec <- c(0, 1, 2, 3, 4, 5, 
             0, 0, 6, 7, 8, 9,
             0, 0, 0, 10, 11, 12,
             0, 0, 0, 0, 13, 14,
             0, 0, 0, 0, 0, 15,
             0, 0, 0, 0, 0, 0)
runDelta <- as.dist(matrix(vec, 6, 6))

runData <- makeMDSData(runDelta)
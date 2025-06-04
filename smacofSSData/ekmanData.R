data(ekman, package = "smacof")
ekmanData <- makeMDSData((1 - ekman)^3)
ekmanLabels <- as.character(attr(ekman, "Labels"))


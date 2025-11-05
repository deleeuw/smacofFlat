source("smacofDataUtilities.R")
source("smacofAuxiliaries.R")
source("smacofTorgerson.R")
source("smacofGuttman.R")
source("smacofElegant.R")
source("smacofPlots.R")
source("smacofSS.R")
source("smacofSSU.R")
source("smacofSSData/ekmanData.R")

h <- smacofSSU(
  ekmanData,
  ndim = 2,
  xinit = NULL,
  ties = 1,
  itmax = 100,
  eps = 1e-10,
  digits = 10,
  width = 15,
  verbose = TRUE,
  ordinal = FALSE
) 
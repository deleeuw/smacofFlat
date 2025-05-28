h <- .C(
  "smacofStepUnweighted",
  nobj = as.integer(nobj),
  ndim = as.integer(ndim),
  ndat = as.integer(ndat),
  iind = as.integer(iind),
  jind = as.integer(jind),
  edis = as.double(edis),
  dhat = as.double(dhat),
  xold = as.double(as.vector(xold)),
  xnew = as.double(as.vector(xnew))
)

h <- .C64(
  "smacofStepUnweighted",
  SIGNATURE = c(
    "integer",
    "integer",
    "integer",
    "integer",
    "integer",
    "double",
    "double",
    "double",
    "double"
  ),
  nobj = as.integer(nobj),
  ndim = as.integer(ndim),
  ndat = as.integer(ndat),
  iind = as.integer(iind),
  jind = as.integer(jind),
  edis = as.double(edis),
  dhat = as.double(dhat),
  xold = as.double(as.vector(xold)),
  xnew = as.double(as.vector(xnew)),
  INTENT = c("r", "r", "r", "r", "r", "rw", "r", "r", "rw"),
  NAOK = TRUE,
  verbose = 0
)
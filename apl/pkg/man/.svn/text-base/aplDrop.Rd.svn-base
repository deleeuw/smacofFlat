\name{aplDrop}
\alias{aplDrop}
\title{
APL Drop
}
\description{
Remove the first (or last) a elements of x according to sign of a
}
\usage{
aplDrop(a, x, drop = FALSE)
}
\arguments{
  \item{a}{an array}
  \item{x}{a vector of integers with length same as dimension of a.  Positive integer specifies number of elements from the head to drop, negative integer specifies number of elemets from the tail to drop.  0 will not drop any element.}
  \item{drop}{ Boolean flag to specify whether to convert the result into a vector when it only has one column.  Default if FALSE.}
}

\references{
G. Helzer. (1989):
\emph{An Encyclopedia of APL}, second edition, I-APL LTD, St. Albans, G.B.
}
\author{
  Jan de Leeuw \email{deleeuw@stat.ucla.edu}
}

\seealso{
    \code{\link[apl]{aplTake}}
}
\examples{
a<-array(1:24,c(2,3,4))
aplDrop(a, c(1,1,0))
}
\keyword{ manip }

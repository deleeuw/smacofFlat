\name{aplTake}
\alias{aplTake}
\title{
APL Take
}
\description{
Takes the elements out of an array
}
\usage{
aplTake(a, x, drop = FALSE)
}
\arguments{
  \item{a}{an array}
  \item{x}{a vector of integers with length same as dimension of a.  Positive integer specifies number of elements from the head to extract, negative integer specifies number of elemets to leave thus all the other elemets that is not specifed by the negative integer will be returned.  0 will return empty value.}
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
    \code{\link[apl]{aplDrop}}
}
\examples{
a<-array(1:24,c(2,3,4))
aplTake(a, c(1,1,2))

}
\keyword{ manip }

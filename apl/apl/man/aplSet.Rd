\name{aplSet}
\alias{aplSet}
\title{
APL Set
}
\description{
Sets a value to an array.
}
\usage{
aplSet(a, b, cell)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{an array}
  \item{b}{a value to set}
  \item{cell}{indicies of cell to set the value to.  Length of indices must match the dimension of the array.}
}

\references{
G. Helzer. (1989):
\emph{An Encyclopedia of APL}, second edition, I-APL LTD, St. Albans, G.B.
}
\author{
  Jan de Leeuw \email{deleeuw@stat.ucla.edu}
}
\seealso{
    \code{\link[apl]{aplGet}}
}
\examples{
a<-array(1:24,c(2,3,4))
aplSet(a,999,c(1,1,1))
}
\keyword{ manip }

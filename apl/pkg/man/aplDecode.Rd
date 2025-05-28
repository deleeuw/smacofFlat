\name{aplDecode}
\alias{aplDecode}
\title{
APL Decode
}
\description{
}
\usage{
    aplDecode(cell, dims)
}

\arguments{
  \item{cell}{}
  \item{dims}{}
}
\details{
If cell is scalar and dims is a vector, then LR is the polynomial 
evaluated at L. This means that if the ri are nonnegative integers less than L, then
LR gives the base-10 equivalent of the base-L number R.

If the arguments cell and dims are vectors of the same length. In the array context,
decode returns the index of element x[b] in an array x with dim(x)=a. Obviously
the R implementation, which uses colum-major ordering, will give results different
from the APL implementation. In APL the expression 3 3 31 2 3 evaluates to 18,
while aplDecode(1:3,rep(3,3)) gives 22.
}
\references{
G. Helzer. (1989):
\emph{An Encyclopedia of APL}, second edition, I-APL LTD, St. Albans, G.B.
}
\author{
  Jan de Leeuw \email{deleeuw@stat.ucla.edu},
  Masanao Yajima \email{yajima@stat.ucla.edu}
}

\seealso{
    \code{\link[apl]{aplEncode}}
}
\examples{

}
\keyword{ manip }

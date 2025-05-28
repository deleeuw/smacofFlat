\name{aplInnerProduct}
\alias{aplInnerProduct}
\title{
APL InnerProduct
}
\description{
Generalized inner product.
}
\usage{
aplInnerProduct(a, b, f = "*", g = "+")
}
\arguments{
  \item{a}{an array of dimension I * J * ... * K }
  \item{b}{an array of dimension K * L * ... * M}
  \item{f}{function to apply across.  Default is "*".}
  \item{g}{function to apply within.  Default is "+".}
}
\details{
This function is generalized version of array inner product \%*\% for multidimensional array.  Array a and b must have same size for the last and the first dimension.  By default it performs the usual inner product where each element of a and corresponding element of b is multiplied and the result is summed over.  But this operation can be generalized by specifying different function for f and g.
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
    \code{\link[apl]{aplOuterProduct}}
}
\examples{
a <- matrix(1:9,3,3)
b <- matrix(9:20,3,4)
aplInnerProduct(a,b)
a%*%b
}
\keyword{ manip }

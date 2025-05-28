\name{aplReduce}
\alias{aplReduce}
\alias{aplRDV}
\title{
APL Reduce
}
\description{
Application of function f across a in dimension specified by k.
}
\usage{
aplReduce(a, k, f = "+")
}
\arguments{
  \item{a}{a vector, an array, or a matrix}
  \item{k}{dimension to apply the function f.}
  \item{f}{function to apply.  Default is sum.}
}
\details{
Function f is applied accross dimension specified by k.  This function is similar to aplScan.  The output is same as the last element of aplScan.
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
    \code{\link[apl]{aplScan}}
}
\examples{
aplReduce(array(1:27,c(3,3,3)),1)
aplReduce(array(1:27,c(3,3,3)),2,"*")
aplReduce(array(1:27,c(3,3,3)),3,"-")
}
\keyword{ manip }

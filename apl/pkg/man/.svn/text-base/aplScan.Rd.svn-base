\name{aplScan}
\alias{aplScan}
\alias{aplSCV}
\title{
APL Scan
}
\description{
Running application of function f across a in dimension specified by k.
}
\usage{
aplScan(a, k, f = "+")
}
\arguments{
  \item{a}{a vector, an array, or a matrix}
  \item{k}{dimension to apply the function f.  If f is "+" and k is 1 then running sum for each column is calculated.}
  \item{f}{function to apply.  Default is sum.}
}
\details{
When f is sum aplScan is array version of cumsum.  It applies the function f sequentially in dimension specified by k.
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
    \code{\link[apl]{aplReduce}}
}

\examples{
 aplScan(array(1:27,c(3,3,3)),1)
 aplScan(array(1:27,c(3,3,3)),2,"*")
 aplScan(array(1:27,c(3,3,3)),3,"-")
}
\keyword{ manip }

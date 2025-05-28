\name{aplTranspose}
\alias{aplTranspose}
\title{
APL Transpose
}
\description{
Reverse the order of the dimensions.
}
\usage{
aplTranspose(a, x = rev(1:aplRank(a)))
}
\arguments{
  \item{a}{an array}
  \item{x}{index}
}
\details{
The monadic aplTranspose(a) and aperm(a) are always the same, they reverse the order of the dimensions.
If x is a permutation of 1:aplRank(a), then aperm(a,x) is actually equal to aplTranspose(a,order(x)). For permutations we could consequently define aplTranspose(a,x) simply as aperm(a,order(x)) (which would undoubtedly
be more efficient as well).  
If x is not a permutation, then aperm(a,x) is undefined, but aplTranspose(a,x) can still be defined in some cases.
If x has aplRank(a) elements equal to one of 1:m, with each of 1:m occurring a least once, then aplTranspose(a,x) has rank m. For obvious reasons dyadic transpose is not used a great deal.}
\value{
 \item{z }{transposed array}
}
\references{
G. Helzer. (1989):
\emph{An Encyclopedia of APL}, second edition, I-APL LTD, St. Albans, G.B.
}
\author{
  Jan de Leeuw \email{deleeuw@stat.ucla.edu},
  Masanao Yajima \email{yajima@stat.ucla.edu}
}

%\seealso{
%}
\examples{
aplTranspose(array(1:48,c(3,4,4)))
}
\keyword{ manip }

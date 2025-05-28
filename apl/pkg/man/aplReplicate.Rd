\name{aplReplicate}
\alias{aplReplicate}
\alias{aplCRV}
\title{
APL Replicate
}
\description{
Array Replication 
}
\usage{
aplReplicate(x, y, k)
}
\arguments{
  \item{x}{an array or a vector
}
  \item{y}{an integer vector of length the same length as the size of kth dimension of x
}
  \item{k}{a direction to replicate x y times.  When x is a vector this is not used.}
}
\details{
This function replicates elements of x sliced in the direction of k according to 
the number specified by vector y.  
}
\references{
G. Helzer. (1989):
\emph{An Encyclopedia of APL}, second edition, I-APL LTD, St. Albans, G.B.
}
\author{
  Jan de Leeuw \email{deleeuw@stat.ucla.edu}
}

%\seealso{}
\examples{
aplReplicate(1:5,c(1,2,3,4,5))
a<-array(1:24,c(2,3,4))
aplReplicate(a,c(2,1,1),2)
}
\keyword{ manip }

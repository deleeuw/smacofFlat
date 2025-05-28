\name{aplRotate}
\alias{aplRotate}
\alias{aplRTV}
\title{
APL Rotate
}
\description{
Rotates an array or a vector 
}
\usage{
aplRotate(a, x, k)
}
\arguments{
  \item{a}{an array or a vector}
  \item{x}{an integer array of dimension size same as a without the kth dimension.
A positive number specified represents the number of elemets in the head of the 
array, in the direction of kth dimension, to be placed at the end of the array in kth direction.  
Negative number represents
the number of elements from the tail, in the direction of kth dimension, of the array 
to be placed at the head of the array in the kth dimension.  
When a is a vector this will not be used and one .}
  \item{k}{When a is a vector, if a positive number is specified, it represents the 
number of elemets in the head of the vector to be placed at the end of the vector. 
On the other hand when a negative number is specified, it represents the number of 
elements from the tail of the vector to be placed at the head of the vector.}
}

%\details{}

\references{
G. Helzer. (1989):
\emph{An Encyclopedia of APL}, second edition, I-APL LTD, St. Albans, G.B.
}

\author{
  Jan de Leeuw \email{deleeuw@stat.ucla.edu}
}

%\seealso{}
\examples{
aplRotate(a=1:5,k=1)
aplRotate(a=1:5,k=-2)
a <- array(1:24,c(2,3,4))
b <- array(1, c(2,4))
aplRotate(a,b,2)
aplRotate(a,-b,2)
}
\keyword{ manip }

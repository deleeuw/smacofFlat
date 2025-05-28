\name{aplOuterProduct}
\alias{aplOuterProduct}
\title{
APL Outer Product
}
\description{
Perform general outer product.
}
\usage{
aplOuterProduct(x, y, f = "*")
}
\arguments{
  \item{x}{ an array}
  \item{y}{ an array}
  \item{f}{ a function that takes 2 parameters }
}

\details{ 
This function evaluates the function f for each combination of the elements in array x and y.
}

\references{
G. Helzer. (1989):
\emph{An Encyclopedia of APL}, second edition, I-APL LTD, St. Albans, G.B.
}
\author{
  Jan de Leeuw \email{deleeuw@stat.ucla.edu}
}
\note{
Currently its same as outer function.
}


\seealso{
    \code{\link[apl]{aplInnerProduct}}
}
\examples{
    x <- 1:9
    y <- 1:9
    # default outer product
    aplOuterProduct(x,y)
    # general outer addition
    aplOuterProduct(x,y,f="+")
    # self defined function
    af<- function(x,y){ x * x + y * y }
    aplOuterProduct(x,y,f="af")
}
\keyword{ manip }

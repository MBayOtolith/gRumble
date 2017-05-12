\name{vectorAngle}
\alias{vectorAngle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate angle between two vecotrs
}
\description{
Function to calculate the angle in radians between two different vectors.
}
\usage{
vectorAngle(vec1, vec2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec1}{
%%     ~~Describe \code{vec1} here~~
}
  \item{vec2}{
%%     ~~Describe \code{vec2} here~~
}
}
\details{
each vector is assumed represent one point in N dimension space (same as length of vector)
}
\value{
Returns a single number that is the angle between teh vectors in radians
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Connor F. White
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

angle in 2 dimension space
vectorAngle(c(0,1), c(1,1))

vec1<-c(0,1,1,4,6,-1)

vec2<-c(1,1,2,1,1,-3)

vectorAngle(vec1, vec2


}

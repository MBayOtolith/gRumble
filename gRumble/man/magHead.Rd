\name{magHead}
\alias{magHead}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heading from Magnetometer
}
\description{
based on the animals pitch and roll the heading directed to north, Pi might needed to be added to the heading for the correct orientation, Might just be the hard iron offset though, blah
}
\usage{
magHead(PR,magxyz,magoff)
}

\arguments{
  \item{PR}{
pitch and roll in radians
}
  \item{magxyz}{
magnetometer readings centered on zero
}
  \item{magOff}{
vector of length three representing the hard iron offset
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns a vector of headings in reference to north in radians
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

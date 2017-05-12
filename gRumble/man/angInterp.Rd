\name{angInterp}
\alias{angInterp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interpolate between angles in Radians using shoretest distance
}
\description{
need to write
}
\usage{
angInterp(estAngles, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estAngles}{
%%     ~~Describe \code{estAngles} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (estAngles, data) 
{
    rotX <- rep(NA, nrow(data))
    rotX[1] <- estAngles[1, 4]
    rotX[nrow(data)] <- estAngles[nrow(estAngles), 4]
    rotX[estAngles[, 1]] <- estAngles[, 4]
    rotY <- rep(NA, nrow(data))
    rotY[1] <- estAngles[1, 3]
    rotY[nrow(data)] <- estAngles[nrow(estAngles), 3]
    rotY[estAngles[, 1]] <- estAngles[, 3]
    rotZ <- rep(NA, nrow(data))
    rotZ[1] <- estAngles[1, 5]
    rotZ[nrow(data)] <- estAngles[nrow(estAngles), 5]
    rotZ[estAngles[, 1]] <- estAngles[, 5]
    rotR <- rep(NA, nrow(data))
    rotR[1] <- estAngles[1, 6]
    rotR[nrow(data)] <- estAngles[nrow(estAngles), 6]
    rotR[estAngles[, 1]] <- estAngles[, 6]
    rotX <- angleGapInterp(rotX, is.na(rotX))
    rotY <- angleGapInterp(rotY, is.na(rotY))
    rotZ <- eventInterp(rotZ, is.na(rotZ))
    rotR <- eventInterp(rotR, is.na(rotR))
    return(cbind(rotY, rotX, rotZ, rotR))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

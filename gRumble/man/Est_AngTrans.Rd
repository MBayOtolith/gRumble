\name{Est_AngTrans}
\alias{Est_AngTrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rotate  data
}
\description{
Deleting as it is currently redundant
}
\usage{
Est_AngTrans(dat, est_angles)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{est_angles}{
%%     ~~Describe \code{est_angles} here~~
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
function (dat, est_angles) 
{
    for (i in 1:nrow(est_angles)) {
        xyz <- xyzTrans(dat[est_angles[i, 1]:est_angles[i, 2], 
            (1:3)], Yr = est_angles[i, 3], Xr = est_angles[i, 
            4], Zr = est_angles[i, 5])
        if (est_angles[i, 6] < 0) {
            xyz[, 1] <- xyz[, 1] * -1
            xyz[, 2] <- xyz[, 2] * -1
        }
        dat[est_angles[i, 1]:est_angles[i, 2], (1:3)] <- xyz
    }
    return(dat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

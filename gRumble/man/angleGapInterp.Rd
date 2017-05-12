\name{angleGapInterp}
\alias{angleGapInterp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interpolate angels across gaps in angle data }

\description{
interpolate between the angles
}
\usage{
angleGapInterp(dat, events)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{events}{
%%     ~~Describe \code{events} here~~
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
function (dat, events) 
{
    loc <- event(events)
    leng <- loc[, 2] - loc[, 1] + 2
    for (i in 1:nrow(loc)) {
        d <- (dat[loc[i, 2] + 1] - dat[loc[i, 1] - 1])
        if (d > pi | d < -pi) {
            ang1 <- dat[loc[i, 1] - 1] - (pi/2)
            ang2 <- dat[loc[i, 2] + 1] - (pi/2)
            angmag <- pi - abs(ang1) + pi - abs(ang2)
            angs <- seq(0, angmag, length.out = leng[i] + 1)
            if (ang1 < ang2) {
                angs <- ang1 + angs * -1
                switc <- which(angs < -pi)
                angs[switc] <- (angs[switc] + 2 * pi)
            }
            else {
                angs <- ang1 + angs
                switc <- which(angs > pi)
                angs[switc] <- (angs[switc] - 2 * pi)
            }
            dat[loc[i, 1]:loc[i, 2]] <- angs[(2:length(angs) - 
                2)] + (pi/2)
        }
        else {
            dif <- d/leng[i]
            if (leng[i] == 2) {
                dat[loc[i, 1]] <- (dat[loc[i, 1] - 1] + dif)
            }
            else {
                dat[loc[i, 1]:loc[i, 2]] <- seq((dat[loc[i, 1] - 
                  1] + dif), (dat[loc[i, 2] + 1] - dif), length.out = length(loc[i, 
                  1]:loc[i, 2]))
            }
        }
    }
    return(dat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

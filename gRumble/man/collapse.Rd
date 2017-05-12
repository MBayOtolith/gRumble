\name{collapse}
\alias{collapse}

\title{
downsample data
}
\description{
Given a time series of data this summarizes the data into the mean value every x point
}
\usage{
collapse(dat, freq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
vector to be downsamples
}
  \item{freq}{
ever x point to be retained
}
}
\details{
Fast is you only want to downsample a small amount, This acutally calcualtes teh running average using filter and then only retains every 'freq' datapoint. For very large window sizes it is not efficient.

For instance if the data is collected at 5 hz, if 'freq=5' then the data would be downsampled to 1 second means, if 'freq=1500', then th it would be 5 minute means.

}
\value{
returns the summarized vector
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
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
val<-rnorm(3600*5,mean=0,sd=0.1)

#Mean value every second
collapse(val,5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

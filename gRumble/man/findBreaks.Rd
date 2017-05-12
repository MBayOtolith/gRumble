\name{findBreaks}
\alias{findBreaks}
\title{
Generate break points for orientation corrections
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
findBreaks(VV = NULL, corrFreq = 30, minsize = NULL, lims = c(-0.05, 0.05), dat_Freq = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{VV}{
vector of vertical velocities at same frequency as the XYZ data
}
  \item{corrFreq}{
frequency in seconds of how often to calulate break points and limits
}
  \item{minsize}{
minimum size of each window
}
  \item{lims}{
minimum and mazimum vertical velocities with in each window
}
  \item{dat_Freq}{
frequency that the data is collected at so convert from seconds to data points
}
}
\details{
This function is designed to make a series of windows over which to orient tags of unknown positions. The frequency of which to make windows and thus subsequency calculate orientation angles is the corrFreq argument. Windows that are calculated more frequently can be more accurate, however, this increases your introduced processing error. The minimal size of each window is set by the minsize argument. Smaller sizes can more rapidly respond to changes in the tag orientation, however, small window sizes can also attribute animal positional changes to stomach movement and the data is overcorrected, and animal positional information is lost. Small windows can be especially problematic for resolving orientation around the Z axis, which is derived from an FFT and relies on looking at signal strength. Windows that are too large are also problematic as stomach movents occur within a window and increase error.

However, by ensuring that the mean virtical velocity of each window is within the limits provided by the lims argument the animal should largely be prone. If the windows mean VV is outside the values provided by lims, the window size will increase until the vertical velocity is within the limits. This limits the possibility that while an animal is ascending it will correct its body pitch to horizontal. 

}
\value{
findBreaks returns a 4 column matrix.  The first column it the mid point of the window and where the correction will be applied, these are equally spaced, by the corrFreq argument. The second column will be the start of the window and the 3rd column will be the end of the window. These points are not eually spaced but represent a start and a end point that provide the midpoint with an VV within the limits. The 4th column is the average vertical velocity within the window.
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

depth<-filter(cumsum(rnorm(5*3600,mean=0,sd=.2)),filter=rep(1,10)/10,sides=2,circular=TRUE)
vv<-c(0,depth[2:length(depth)]- depth[1:(length(depth)-1)])
vv<-filter(vv,filter=rep(1,5)/5, sides=2,circular=TRUE)


breaks<-findBreaks(VV=vv,corrFreq=5,minsize=15,lims=c(-0.001,0.001),dat_Freq=5)
windSizes<-(breaks[,3]-breaks[,2]) / 5

plot(depth, type="l")
points(depth[breaks[,1]] ~ breaks[,1],cex=scal(windSizes, scal=c(1,5)))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

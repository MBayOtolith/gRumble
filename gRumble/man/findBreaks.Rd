\name{findBreaks}
\alias{findBreaks}
\title{
Generate break points for orientation corrections
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
An indexing function to identifying qualifying periods ("windows") of horizontal/prone behavior where the static acceleration of gravity can be isolated
}
\usage{
findBreaks(VV = NULL, corrFreq = 30, minsize = NULL, lims = c(-0.05, 0.05), dat_Freq = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{VV}{
vector of Vertical Velocities (VV) at same length as raw XYZ data
}
  \item{corrFreq}{
resolution in seconds of how often to calulate break points and limits
}
  \item{minsize}{
minimum size of each window for identifying horizontal/prone behavior
}
  \item{lims}{
minimum and maximum thresholds of VV limits for definining assumptions about presence of horizontal/prone posture within the window
}
  \item{dat_Freq}{
temporal resolution (in seconds) of raw data for converting from seconds to data points
}
}
\details{
This function indexes a series of behavioral windows over which average VV is near zero, the animal is in a relatively horizontal/prone posture, and the sensor data can be re-oriented.  The animal's posture is estimated by calculating the mean VV within windows of increasing size, until mean VV falls within the threshold limits set by the lims arguement.

The temporal resolution of this angle estimation & re-orientation procedure is controlled by the arguement, corrFreq, which determines how often (in seconds) to calculate correction angles.  While windows that are calculated more frequently can be more accurate, this however can also  increase introduced processing error.

The minimal possible size of a window, in seconds, is set by the minsize argument and should be set with significant consideration.  Smaller window sizes are desirable because they can respond more rapidly to changes in tag orientation.  However, overly small windows can lead to overcorrection if stomach movements are attributed to posture/positional changes and are additionally problematic for resolving orientation around the Z-axis.  Z-axis corrections are derived from an FFT which relies on signal strengths.  On the other hand, windows that are too large are problematic by introducing error from stomach movements that occur within a window.

The lims arguement defines the acceptable thresholds of assuming the animal's posture is prone or horizontal.  Window size is increased until the associated mean VV falls within these thresholds, thus limiting the possibility of mis-positioning a solely ascending or descedning animal to horizontal.

}
\value{
findBreaks returns a matrix with 4 columns.  Column 1 marks the midpoint of each window, where corrections are to be applied in an equally-spaced manner set by the corrFreq arguement.  Columns 2 and 3 mark the start/end indices of each window, respectively, but have variable spacing in order to ensure a mean VV over the window's duration within the specified VV threshold limits.  Column 4 is the average VV within this window.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Connor F. White (code, text)
Jerry H. Moxley (text)
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

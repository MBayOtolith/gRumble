\name{EstAngles}
\alias{EstAngles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generate angles to correct tag orientation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EstAngles(dat, vv, breaks, dat_Freq = 5, signal = c(0.3, 0.7), wind = 7, Zby = pi/180, lim = c(-pi/2, pi/2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
3 column xyz acceleration data
}
  \item{vv}{
vertical velocity the same length as dat
}
  \item{breaks}{
2 column matrix representing the start and end of each window
}
  \item{dat_Freq}{
the frequency the data is collected at in Hz
}
  \item{signal}{
Frequency range of the locomotion signal for the correction around Z
}
  \item{wind}{
window size to smooth the data by to isolate the dynamic acceleration from locomation (in seconds)
}
  \item{Zby}{
resolution of corrections around the z axis in radians
}
  \item{lim}{
between subsequent windows what are the limits the tag can shift by
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a data frame that is 6 columns wide and the same number of rows as the breaks argument. The columns 1 and 2 are teh start and end of the windows. column 3 is the rotation the Y axis, column 4 is the rotation around the x axis, and column 5 is the rotation around the Z axis. Column 6 is the coefficient of the regression between VV and body pitch, to indicate if teh tag was corrected in the wrong direction along the anterior posterior axis.
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

#Making Depth and VV data over a 1 hour period
depth<-filter(cumsum(rnorm(5*3600,mean=0,sd=.2)),filter=rep(1,10)/10,sides=2,circular=TRUE)
vv<-c(0,depth[2:length(depth)]- depth[1:(length(depth)-1)])
vv<-filter(vv,filter=rep(1,5)/5, sides=2,circular=TRUE)

#Making where to have break points be
breaks<-findBreaks(VV=vv,corrFreq=5,minsize=15,lims=c(-0.001,0.001),dat_Freq=5)

#Simulating some acceleration data
accelX<-rnorm(length(depth),mean=0,sd=0.05)
accelY<-rep(sin(seq(0,2*pi,length.out=15)),length.out=length(depth))/10 + rnorm(length(depth),mean=0,sd=0.01)
accelZ<-runif(n= length(depth), max = 1,min =.9)

xyz<-cbind(accelX, accelY, accelZ)

#Random walk of angles
angles<-cbind(cumsum(filter(rnorm(length(depth),0,sd=0.001),filter=rep(1,15)/15, sides=2,circular=TRUE)),
              cumsum(filter(rnorm(length(depth),0,sd=0.001),filter=rep(1,15)/15, sides=2,circular=TRUE)),
              cumsum(filter(rnorm(length(depth),0,sd=0.001),filter=rep(1,15)/15, sides=2,circular=TRUE)))

#rotating the acceleratoin data to simulate it in a stomach
xyz<-RotData(angCorr = cbind(angles,1),xyz=xyz)

#estimating the angle of rotation
est_ang <- EstAngles(dat=xyz,vv = vv,breaks=breaks[,c(2,3)],
                     dat_Freq = 5, signal = c(.3,.7), wind = 6, Zby = pi/90,lim=c(-pi/4,pi/4))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

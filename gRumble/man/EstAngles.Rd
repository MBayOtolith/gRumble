\name{EstAngles}
\alias{EstAngles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate angles to correct tag orientation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Estimation function to calculate the X/Y/Z angles within each window of horizontal/prone behavior that can correct data to the animal's frame of reference;
}
\usage{
EstAngles(dat, vv, breaks, dat_Freq = 5, signal = c(0.3, 0.7), wind = 7, Zby = pi/180, lim = c(-pi/2, pi/2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
3-column matrix of tri-axial acceleration data (X/Y/Z axes) at length of raw data
}
  \item{vv}{
vector of Vertical Velocities (VV) at length of raw XYZ data
}
  \item{breaks}{
2-column matrix of start/end indices of each window identified by findBreaks()
}
  \item{dat_Freq}{
the frequency (Hz) of raw data at time of collection
}
  \item{signal}{
range of frequencies (Hz??) for identifying locomotion signals in corrections around Z-axis
}
  \item{wind}{
window size to smooth the data by to isolate the dynamic acceleration from locomation (in seconds)
}
  \item{Zby}{
resolution of corrections around Z-axis (radians)
}
  \item{lim}{
max and min values setting limits of relative angular change between subsequent windows (radians)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns 6-column data frame that is the same length as the breaks arguement.  Columns 1 and 2 are the start/end indices of each window indentified in the breaks arguement.  Columns 3, 4, and 5 are the rotations around the X ("roll"), Y ("pitch"), and Z ("yaw") axes respectively.  Column 6 is the coefficient of regression between VV and body pitch that indicates alignment or malalignment of the tag with the animal's anterior/posterior axis.

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

#rotating the acceleration data to simulate it in a stomach
xyz<-RotData(xyz = xyz, Xr = angles[,1],Yr = angles[,2],Zr = angles[,3])

#estimating the angle of rotation
est_ang <- EstAngles(dat=xyz,vv = vv,breaks=breaks,
                     dat_Freq = 5, signal = c(.3,.7), wind = 6, Zby = pi/90,lim=c(-pi/4,pi/4))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

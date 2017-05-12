\name{ZrotAng}
\alias{ZrotAng}

\title{
Find the angle of rotation to maximize signal
}
\description{
perform a series of rotations around the z axis in order to find the angle that maximizes the desired signal in the y axis.
}
\usage{
ZrotAng(xyz, dat_Freq = 5, signal = c(0.3, 0.7), wind = 7, by = pi/360, head = 0, lim = c(0, pi))
}


\arguments{
  \item{xyz}{
3 column xyz data
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
  \item{by}{
resolution of corrections around the z axis in radians
}
  \item{lim}{
how far from the current heading shouldsignal strength be evaluated for, in radians
}
  \item{head}{
current heading in radians
}
  \item{weighting}{
TRUE/FALSE, should the estimate from the FFT be fused with the heading based on weighting factor
}
}
\details{
This fuction performs a series of rotations of the XYZ data around the Z axis. The number of rotations is determined by both the 'by' argument and the 'lim' argument. The XYZ data is rotated between headying - lims[1] through heading + lims[2], at a resolution provided by the 'by' argument. At each rotation, the Y axis acceleration is smoothed using a moving window (wind argument) and is then subtracted from the raw acceleration of the Y axis to produce the synamic acceleration which would indicate locomotion. An fast fourier transform (FFT) is run on the dynamic acceleration, and the magnitude of all frequencies within the 'signal' range are averaged. The rotation angle that maximizes the mean of the FFT magnitude within the signal range is returned. 

weighting is determined by
p<-(max(mags)-median(mags))/max(mags)
theta<- ((p*theta) + ((1-p)*head))
}
\value{
This function returns the angle (in radians) that maximizes the FFT magnitude for signals in the range provided by within the 'signal' argument.  
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

#Simulating some acceleration data
accelX<-rnorm(1000,mean=0,sd=0.05)
accelY<-rep(sin(seq(0,2*pi,length.out=15)),length.out=1000)/10 + rnorm(1000,mean=0,sd=0.01)
accelZ<-runif(n= 1000, max = 1,min =.9)

xyz<-cbind(accelX, accelY, accelZ)

#Rotate by 45 degrees
xyz<- xyz \%*\% Zb(pi/4) 

theta<-ZrotAng(xyz, dat_Freq = 5, signal = c(0.3, 0.7), wind = 7, by = pi/360, head = 0, lim = c(0, pi))

theta* (180/pi)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

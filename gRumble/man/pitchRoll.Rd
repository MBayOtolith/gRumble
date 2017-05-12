\name{pitchRoll}
\alias{pitchRoll}

\title{
Generate Pitch and Roll angles
}
\description{
From XYZ accleration data generate pitch and roll angles. The acceleration data is already expected to be only comprised of static acceleration due to gravity 
}
\usage{
pitchRoll(xyz, degrees= TRUE)
}

\arguments{
  \item{xyz}{
3 column XYZ acceleration data
}
  \item{degrees}{
TRUE/FALSE if you would like the output returned in degrees vs radians
}
}
\details{
X is assumed to be along the anterior posterior axis of the animal while Z is the up down axis assosciated with gravity. Uses atan so will only go to +- 90 degrees and cannot resolve in the between 2nd and 4th quadrant or the 1st and 3r quadrants.

If you want pitch and roll values greater than +- 90 use Pitch Roll 2 which uses atan2

}
\value{
returns a 2 column matrix with the same number of rows as the xyz input. The angles are returned in radians by default. 




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
#make some data with gravity pulling down
xyz<-c(0,0,-1)

#rotate the gravity vector by pi/6 in the x and y direction
xyzRot<-xyzTrans(matrix(c(0,0,1),nrow=1),Yr=pi/6,Xr=pi/6,Zr=0)

pitchRoll(xyzRot)



#generating some data and adding some noise which would be indicative of dynamic animal movement
x<-seq(0, 0.5,length.out=50) + rnorm(50,mean=0,sd=.1)
y<-seq(0.5, 0,length.out=50) + rnorm(50,mean=0,sd=.1)
z<-sqrt((1-(seq(0, 0.5,length.out=50)^2 + seq(0.5, 0,length.out=50)^2))) + rnorm(50, mean= 0, sd=.1)

xyz<-cbind(xyz)
#smooth the data over 5 data points to get rid of the activity added
xyzG<-Gsep(xyz,filt=rep(1,5)/5)

pitchRoll(xyzG[,1:3])

}

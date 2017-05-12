\name{Gsep}
\alias{Gsep}

\title{
Seperate Static and Dynamic Acceleration
}
\description{
given raw acceleration seperate static acceleration due to gravity and dynamic acceleration due to body movement based on a moving average. Additionally with generate overall dynamic body acceleration (ODBA) based on the dynamic acceleration
}
\usage{
Gsep(xyz, filt = rep(1, 5 * 5)/(5 * 5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xyz}{
  3 column xyz acceleration data
}
  \item{filt}{
  the coefficients for a moving window, should sum to 1.
}
}
\details{
See filter in the base stats package for information on construction, defaults to a 25 point moving average, which would be a 5 second window if the data is recorded at 5 hertz.
}
\value{
Returns a 7 columns matrix with the same number of rows as XYZ input. Columns 1:3 are the acceleration due to gravity along the x,y and z axis and are indicative of body orientation. Columns 4:6 are the dynamic acceleration and are indicative of acceleration due to the movement of the organism. Column 7 is ODBA, and is the sum of the absolute value of the dynamic acceleration 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Connor F. White
}


\seealso{
\code{\link{stats:filter}}
}
\examples{
#generating some data and adding some noise which would be indicative of dynamic animal movement
x<-seq(0, 0.5,length.out=50) + rnorm(50,mean=0,sd=.1)
y<-seq(0.5, 0,length.out=50) + rnorm(50,mean=0,sd=.1)
z<-sqrt((1-(seq(0, 0.5,length.out=50)^2 + seq(0.5, 0,length.out=50)^2))) + rnorm(50, mean= 0, sd=.1)

xyz<-cbind(xyz)
#smooth the data over 5 data points to get rid of the activity added
xyzG<-Gsep(xyz,filt=rep(1,5)/5)

#smooth the data over 10 data points to get rid of the activity added
xyzG<-Gsep(xyz,filt=rep(1,10)/10)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

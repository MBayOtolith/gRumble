\name{Xb}
\alias{Xb}

\title{
rotation matrix around X axis
}
\description{
generates a 3x3 rotation matrix to rotate XYZ data around the X axis. This is a Roll rotation
}
\usage{
Xb(b)
}
\arguments{
  \item{b}{
  angle of rotation in radians
}
}
\details{
This was created for rotating XYZ acceleration data to a global frame of reference, x value will  remain the same while Y and Z values will change
}
\value{
generates a 3x3 rotation matrix that can be multiplied by a 1x3 vector to rotate the vector around the first coordiate. 
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
\code{\link{Zb}}, \code{\link{Yb}}
}
\examples{
#The rotation Matrix
Xb(pi/2)

#point to rotate
xyz<-c(1,0,1)

#rotating xyz 90 degrees (pi/2) around the Z axis
xyz \%*\% Xb((pi/2))
}

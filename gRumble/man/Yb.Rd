\name{Yb}
\alias{Yb}

\title{
rotation matrix around Y axis
}
\description{
generates a 3x3 rotation matrix to rotate XYZ data around the Y axis. This is a pitch rotation
}
\usage{
Yb(b)
}
\arguments{
  \item{b}{
  angle of rotation in radians
}
}
\details{
This was created for rotating XYZ acceleration data to a global frame of reference, y value will  remain the same while X and Z values will change
}
\value{
generates a 3x3 rotation matrix that can be multiplied by a 1x3 vector to rotate the vector around the second coordiate. 
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
\code{\link{Zb}}, \code{\link{Xb}}
}
\examples{
#The rotation Matrix
Zb(pi/2)

#point to rotate
xyz<-c(1,0,1)

#rotating xyz 90 degrees (pi/2) around the Z axis
xyz \%*\% Yb((pi/2))
}


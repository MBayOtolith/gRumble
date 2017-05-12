\name{Zb}
\alias{Zb}

\title{
rotation matrix around Z axis
}
\description{
generates a 3x3 rotation matrix to rotate XYZ data around the Z axis. This is a yaw rotation
}
\usage{
Zb(b)
}
\arguments{
  \item{b}{
  angle of rotation in radians
}
}
\details{
This was created for rotating XYZ acceleration data to a global frame of reference, z value will  remain the same while X and Y values will change
}
\value{
generates a 3x3 rotation matrix that can be multiplied by a 1x3 vector to rotate the vector around the third coordiate.
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
\code{\link{Xb}}, \code{\link{Yb}}
}
\examples{
#The rotation Matrix
Zb(pi/2)

#point to rotate
xyz<-c(1,0,1)

#rotating xyz 90 degrees (pi/2) around the Z axis
xyz \%*\% Zb((pi/2))
}


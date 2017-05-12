\name{xyzTrans}
\alias{xyzTrans}

\title{
rotate xyz data
}
\description{
Rotate a vector of length three representing XYZ data around 0 in the X, Y and Z direction 
}
\usage{
xyzTrans(xyz,  Xr=0, Yr=0, Zr=0)
}

\arguments{
  \item{xyz}{
vector representing three coordinates to be rotated
}
  \item{Yr}{
rotation anlge in radians around the Y axis 
}
  \item{Xr}{
rotation anlge in radians around the X axis
}
  \item{Zr}{
rotation anlge in radians around the Z axis
}
}
\details{
this is simply a wrapper function for using the Yb, Xb, and Zb function. See RotData for applying through a large dataset. The order of rotations is roll, then pitch then yaw, which is rotatoins around X then around Y and then around X, which is the equivilant of Xb() \%*\% Yb() \%*\% Zb(). The order of rotations is important, so when determining where you want the final data point enure that you derive them in that order.

THIS USED TO BE IN A DIFFERENT ORDER SO WILL NEED TO BE DOUBLE CHECKED.
}
\value{
returns a vector of length three that represeents the transformed x y and z coordinates.
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
\code{\link{RotData}}
}
\examples{
xyz<-c(1,0,1)

#rotation of 0 produces original data
xyzTrans(xyz, Yb=0, Xb=0, Zb=0)

#Rotate aroun the by 45 degrees.
xyzTrans(xyz, Yb=0, Xb=0, Zb=pi/4)

xyzTrans(xyz, Yb=pi/4, Xb=0, Zb=pi/4)

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

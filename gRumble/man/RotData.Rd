\name{RotData}
\alias{RotData}

\title{
Rotate XYZ data
}
\description{
Rotate a series of XYZ data points by an equal number of XYZ rotation angles
}
\usage{
RotData(xyz,Xr,Yr,Zr)
}

\arguments{
  \item{xyz}{
3 column XYZ data to be rotated 
}
  \item{Xr}{
rotation around X, roll
}
  \item{Yr}{
rotation around Y, pitch
}
  \item{Zr}{
rotation around Z, Yaw
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns a three column xyz matrix of the rotated data.
}
\author{
Connor F. White
}


\examples{

#Generate XYZ location data
xyz<-matrix(rep(c(0,0,-1),10),ncol=3)

#Rotate around the Z axis by 45 degrees (pi/4)
angs<-cbind(rep(0,10),rep(0,10),rep(pi/4,10))

RotData(xyz=xyz,Xr=angs[,1],Yr=angs[,2],Zr=angs[,3])

}

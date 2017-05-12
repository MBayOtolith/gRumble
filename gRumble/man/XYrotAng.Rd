\name{XYrotAng}
\alias{XYrotAng}

\title{
Angles to rotat XYZ to level
}
\description{
function to calculate the angle of rotations around teh X and the Y directions to make the locations xyz level.

}
\usage{
XYrotAng(xyz)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xyz}{
%%     ~~Describe \code{xyz} here~~
}
}
\details{
if xyz is more than 1 row, then the mean of each column will be used to estimate the correction angles. This corrects so that X and Y are in the horizontal direction, and should have a value of 0, while the z value should be pointed down in the negative direction.
}
\value{
returns 2 angles Xr and Yr. Xr is the roll angle and Yr is the pitch angle
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

xyz<-c(0,0,-1)

Zr=0*pi/2
Xr=30*pi/180
Yr=-15*pi/180

xyz1<-xyzTrans(xyz,Xr,Yr,Zr)

angs<-XYrotAng(xyz1)

#should be 0,0,-1
xyzTrans(xyz1,Xr=angs[1],Yr=angs[2],Zr=Zr)

}

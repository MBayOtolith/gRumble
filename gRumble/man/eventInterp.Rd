\name{eventInterp}
\alias{eventInterp}

\title{
Interpret during events
}
\description{
Identify events such as strings of NA and then linearly interpolating across them
}
\usage{
eventInterp(dat, event, clip = "none", na.rm = FALSE, nas = FALSE, duration = NULL)
}

\arguments{
  \item{dat}{
data set that the interpolated dataset will go into
}
  \item{x}{
A logical statement in which the output is TRUE or FALSE
}
  \item{clip}{
Defaults to "none", If set to "buffer" a false if placed on each end of the string,if set to TRUE the data is clipped to the first and last time when the event did not occur
}
  \item{na.rm}{
TRUE/FALSE should NAs be handled
}
  \item{nas}{
if na.rm is set to true, what should NAs be replaced with TRUE/FALSE
}
  \item{duration}{
Numeric number selecting the minimum duration of an event
}
}
\details{
See description of event
}
\value{
Returns a vector the same length as the original dataset, yet values are linearly interpolated where the event is true
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
event, eventMean
}
\examples{
data<-c(.5,.6,.75,NA,NA,NA,.6,.9,.8,.7)

eventInterp(dat=data,x=is.na(data))

}
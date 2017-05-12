\name{event}
\alias{event}

\title{
Identify events
}
\description{
Finds the first and last location of an event
}
\usage{
event(x, ends = 0, nas = TRUE, duration = c(NA,NA))
}

\arguments{
  \item{x}{
A logical statement in which the output is TRUE or FALSE
}
  \item{ends}{
Defaults to 0, If set to 2, a false if placed on each end of the string,if set to 1 the data is clipped to the first and last time when the event did not occur
}
  \item{nas}{
What should NAs in x be replaced with TRUE/FALSE
}
  \item{duration}{
Numeric number selecting the minimum and max duration of selected events in datapoints, NA is infinity
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 Returns a matix of two columns
 
 \item{start }{the location in the original vector at which the event starts}
 \item{end }{the location in the original vector at which the event stops}
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
eventInterp
}
\examples{

x<-c(.5,.5,.75,.25,.2,.25,.6,.9,.1,.7)
x.e<-event(x<.3)

#Excluding the .1 where there is only 1 value below .3
x.e<-event(x<.3,duration=c(1,NA))


#If we were interested in values greater than .3 we would need
# to alter the data as the the first and last values are true.

#Starts searching for events from the first time a FALSE occurs
x.e<-event(x>.3,ends=1)

#Starts searching for events from the very beginning
x.e<-event(x>.3,ends=2)

}
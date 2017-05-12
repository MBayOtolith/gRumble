\name{eventFunc}
\alias{eventFunc}

\title{
executes function over events
}
\description{
Executes a funtion over every event and returns the output for each event
}
\usage{
eventFunc(dat, events, fun=function(x){mean(x,na.rm=TRUE)})
}


\arguments{
  \item{dat}{
the vector that you wish to take values from
}
  \item{event}{
an output from the event function
}
  \item{fun}{
a function to execute for each event
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns the output of the suppied funtion for each event
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
events
}
\examples{
data<-rnorm(100)
on<-rbinom(prob=.3,size=1,n=100)

instances<-event(on==0,clip="buffer")

eventFunc(dat=data,events=instances,fun=function(x){mean(x,na.rm=TRUE)})

eventFunc(dat=data,events=instances,fun=function(x){length(x)})
}



magHead<-function(PR,magxyz){
  yH<-magxyz[,3]*sin(PR[,2]) - magxyz[,2]*cos(PR[,2])
  xH<-magxyz[,1]*cos(PR[,1]) + magxyz[,2]*sin(PR[,1])*sin(PR[,2])+ magxyz[,3]*sin(PR[,1])*cos(PR[,2])
  Heading<-atan2(-yH,xH)
  return(Heading)
}
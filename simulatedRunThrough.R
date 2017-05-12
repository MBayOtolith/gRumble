library(devtools)
install_github("connorfwhite/gRumble/gRumble")
library(gRumble)



#Making Depth and VV data over a 1 hour period
depth<-filter(cumsum(rnorm(5*3600,mean=0,sd=.2)),filter=rep(1,10)/10,sides=2,circular=TRUE)
vv<-c(0,depth[2:length(depth)]- depth[1:(length(depth)-1)])
vv<-filter(vv,filter=rep(1,5)/5, sides=2,circular=TRUE)

#Making where to have break points be
breaks<-findBreaks(VV=vv,corrFreq=5,minsize=15,lims=c(-0.001,0.001),dat_Freq=5)

#Simulating some acceleration data
accelX<-vv/max(vv) *.5 + rnorm(length(depth),mean=0,sd=0.05)
accelY<-rep(sin(seq(0,2*pi,length.out=15)),length.out=length(depth))/10 + rnorm(length(depth),mean=0,sd=0.01)
accelZ<-runif(n= length(depth), max = -0.9 ,min =-1.1)

xyz0<-cbind(accelX, accelY, accelZ)

#Random walk of angles
angles<-cbind(cumsum(filter(rnorm(length(depth),0,sd=0.01),filter=rep(1,15)/15, sides=2,circular=TRUE)),
              cumsum(filter(rnorm(length(depth),0,sd=0.01),filter=rep(1,15)/15, sides=2,circular=TRUE)),
              cumsum(filter(rnorm(length(depth),0,sd=0.01),filter=rep(1,15)/15, sides=2,circular=TRUE)))

#rotating the acceleratoin data to simulate it in a stomach
xyz<-RotData(xyz = xyz0, Xr = angles[,1],Yr = angles[,2],Zr = angles[,3])

#estimating the angle of rotation
est_ang <- EstAngles(dat=xyz,vv = vv,breaks=breaks,
                     dat_Freq = 5, signal = c(.3,.7), wind = 6, Zby = pi/90,lim=c(-pi/4,pi/4))

#Not Strictly Nessesary, but smoothing the coefficients 
est_ang[,7]<-filter(est_ang[,7],filter=rep(1,5)/5,circular=TRUE,sides=2)

#If the coefficient between vertical velocity and pitch is negative then rotate Zr by 180 degrees
est_ang[which(est_ang[,7] < 0 ),6]<-est_ang[which(est_ang[,7] < 0 ),6] + pi


#backfill and interpolate the angles
angCorr<-angInterp(estAngles = est_ang,data=xyz)

#Rotate the messed up stomach data back to try and estimate the original data.
corrected<-RotData(xyz=xyz,Xr = angCorr[,1], Yr = angCorr[,2], Zr = angCorr[,3])


#Plotting the corrected data
plot(corrected[,1],type="l",col="black",ylim=c(-1.2,1.2),xlim=c(-1.2,1.2))
lines(corrected[,2],type="l",col="red")
lines(corrected[,3],type="l",col="blue")


plot(xyz[,1]~xyz0[,1],ylim=c(-1.2,1.2),xlim=c(-1.2,1.2))
points(corrected[,1]~xyz0[,1],col="red")

plot(xyz[,2]~xyz0[,2],ylim=c(-1.2,1.2),xlim=c(-.2,.2))
points(corrected[,2]~xyz0[,2],col="red")


plot(xyz[,3]~xyz0[,3],ylim=c(-1.2,1.2),xlim=c(-1.2,1.2))
points(corrected[,3]~xyz0[,3],col="red")





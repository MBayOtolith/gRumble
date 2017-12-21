
Xb<-function(b){
  mat.out<-rbind(c(1,0,0),
                 c(0,cos(b),sin(b)),
                 c(0,-sin(b),cos(b)))
  return(mat.out)
}


Yb<-function(b){
  mat.out<-rbind(c(cos(b),0,sin(b)),
                 c(0,1,0),
                 c(-sin(b),0,cos(b)))
  return(mat.out)
}


Zb<-function(b){
  mat.out<-rbind(c(cos(b),sin(b),0),
                 c(-sin(b),cos(b),0),
                 c(0,0,1))
  return(mat.out)
}


XYZb<-function(Xr,Yr,Zr){
  rotMat<-Xb(Xr) %*% Yb(Yr) %*% Zb(Zr)
  return(rotMat)
}

ZXYb<-function(Xr,Yr,Zr){
  rotMat<-Zb(Zr) %*% Yb(Yr) %*% Xb(Xr)
  return(rotMat)
}


xyzTrans<-function(xyz,Xr=0,Yr=0,Zr=0){
  xyz<-xyz %*% XYZb(Xr,Yr,Zr)
  return(xyz)
}

#rotate so that x and y are perpendicular to gravity
XYTrans<-function(xyz,Yr,Xr){  
  xyz<-Yb(Yr) %*% t(xyz)
  xyz<-Xb(Xr) %*% xyz
  xyz<-t(xyz)
  return(xyz)
}

XYrotAng<-function(xyz){
  if(length(xyz)>3){
    pos<-c(mean(xyz[,1],na.rm=TRUE),
         mean(xyz[,2],na.rm=TRUE),
         mean(xyz[,3],na.rm=TRUE))
  } else{pos=xyz}

  roll = atan2(pos[2],(pos[3]))
  pitch = atan2((-pos[1]),((pos[2] * sin(roll)) + (pos[3]*cos(roll))))
  Xr=roll
  Yr=pi-pitch

  return(c(Xr,Yr))
}




ZrotAng<-function(xyz,dat_Freq=5, signal=c(.3,.7),wind=7, by=pi/360, head=0,lim=c(-pi/2,pi/2),weighting=FALSE){
  if(any(is.na(xyz))){
    stop("NAs in data")
  }
  if(nrow(xyz)%%2 > 0){
    xyz<-xyz[-nrow(xyz),]
  }

  zrot<-head + seq(lim[1],lim[2],by=by)

  ###Output frequencies for the FFT
  fs<-(0:((nrow(xyz)/2) - 1)) * dat_Freq/nrow(xyz)
  locs<-which(fs>signal[1] & fs<signal[2])

  mags<-NULL
  for(x in 1:length(zrot)){
    xyz_cor<-Zb(zrot[x])%*%t(xyz)
    xyz_cor<-t(xyz_cor)

    filt<-filter(xyz_cor[,2],filter=rep(1,dat_Freq*wind)/(dat_Freq*wind), sides = 2,circular = TRUE)
    dyn<-xyz_cor[,2]-filt
    y.fft<-Mod(fft(dyn))
    mags<-cbind(mags,y.fft[locs])

  }
  theta<-zrot[which(mags== max(mags),arr.ind=TRUE)[2]]
  if(weighting==TRUE){
    p<-(max(mags)-median(mags))/max(mags)
    theta<- ((p*theta) + ((1-p)*head))
  }

  return(theta)
}



findBreaks<-function(VV = NULL,corrFreq= 30,minsize=NULL,lims=c(-.05,.05), dat_Freq=5){
  size<-30 #Freq

  if(is.null(minsize)){
    minsize <- CorrFreq
  }

  exp<-floor((minsize- corrFreq) /2)* dat_Freq

  loc<-seq(1,length(VV), by = corrFreq*dat_Freq)
  Breaks<-matrix(0,nrow=(length(loc)-1),ncol=4)

  for(i in 1:(length(loc)-1)){
    start<-loc[i] - exp
    end<-loc[i+1] + exp

    if(start<1){
      start<-1
    }
    if(end>length(VV)){
      end<-length(VV)
    }

    d<-mean(VV[start:end])
    while((d > lims[2]) | (d < lims[1])){
      start<-start-5
      end<-end+5
      if(start<1){
        start<-1
      }
      if(end>length(VV)){
        end<- length(VV)
      }
      d<-mean(VV[start:end])
    }
    Breaks[i,1]<-floor((loc[i] +loc[i+1] )/ 2)
    Breaks[i,2]<-start
    Breaks[i,3]<-end
    Breaks[i,4]<-d
  }
  colnames(Breaks)<-c("Location","Start","End","VV")
  return(Breaks)
}






EstAngles<-function(dat,vv, breaks,dat_Freq = 5,signal = c(.3,.7),wind = 7,Zby = pi/180, lim=c(-pi/2,pi/2),weighting=FALSE){
  est_angles<-matrix(0,nrow=nrow(breaks),ncol=7)

  #'Initialize' the model to calculate a heading
  dtemp<-dat[breaks[1,2]:breaks[1,3],]
  vv1<-vv[breaks[1,2]:breaks[1,3]]
  dtemp<-dtemp[complete.cases(dtemp[,1:3]),]
  vv1<-vv1[complete.cases(dtemp[,1:3])]

  XYr<-XYrotAng(dtemp[,1:3])
  Zr<-ZrotAng(xyzTrans(as.matrix(dtemp[,1:3]),Xr=XYr[1],Yr=XYr[2],Zr=0),
              dat_Freq = dat_Freq,signal = signal,wind = wind,by = Zby,head=0,lim=c(0,pi),weighting=weighting)

  for(i in 1:nrow(breaks)){

    dtemp<-dat[breaks[i,2]:breaks[i,3],]
    vv1<-vv[breaks[i,2]:breaks[i,3]]
    dtemp<-dtemp[complete.cases(dtemp[,1:3]),]
    vv1<-vv1[complete.cases(dtemp[,1:3])]

    XYr<-XYrotAng(dtemp[,1:3])
    Zr<-ZrotAng(xyzTrans(as.matrix(dtemp[,1:3]),Xr=XYr[1],Yr=XYr[2],Zr=0),
                dat_Freq = dat_Freq,signal = signal,wind = wind,by = Zby,head=Zr,lim=lim,weighting=weighting)

    xyz<-xyzTrans(as.matrix(dtemp[,1:3]),Xr=XYr[1],Yr=XYr[2], Zr=Zr)

    p<-filter(xyz[,1],filter = rep(1,50)/(50), sides = 2,circular = TRUE)
    m<-lm(p ~ vv1)

    est_angles[i,1]<-breaks[i,1]
    est_angles[i,2]<-breaks[i,2]
    est_angles[i,3]<-breaks[i,3]
    est_angles[i,4]<-XYr[1]
    est_angles[i,5]<-XYr[2]
    est_angles[i,6]<-Zr
    est_angles[i,7]<-coef(m)[2]
  }
  colnames(est_angles)<-c("Location","Start","End","Xr","Yr","Zr","P_VV")
  return(est_angles)
}




angleGapInterp<-function(dat,events){
  loc<-event(events)
  leng<-loc[,2] - loc[,1]+2
  for(i in 1 : nrow(loc)){
    d<-(dat[loc[i,2]+1]-dat[loc[i,1]-1])
    if(d > pi| d < -pi){
      ang1<-dat[loc[i,1]-1]- (pi/2)
      ang2<-dat[loc[i,2]+1]- (pi/2)

      angmag<-pi-abs(ang1) + pi-abs(ang2)
      angs<-seq(0,angmag,length.out=leng[i]+1)

      if(ang1 < ang2){
        angs<-ang1+angs*-1
        switc<-which(angs < -pi)
        angs[switc]<-(angs[switc] + 2*pi)
      }else{
        angs<-ang1+angs
        switc<-which(angs > pi)
        angs[switc]<-(angs[switc] - 2*pi)
      }
      dat[loc[i,1]:loc[i,2]]<-angs[(2:length(angs)-2)] +(pi/2)
    }else{
      dif<-d/leng[i]
      if(leng[i]==2){
        dat[loc[i,1]]<-(dat[loc[i,1]-1]+dif)
      }else{
        dat[loc[i,1]:loc[i,2]]<-seq((dat[loc[i,1]-1]+dif),(dat[loc[i,2]+1]-dif), length.out=length(loc[i,1]:loc[i,2]))
      }
    }
  }
  return(dat)
}


angInterp<-function(estAngles, data){
  Xr<-rep(NA,nrow(data))
  Xr[1]<-estAngles[1,4]
  Xr[nrow(data)]<-estAngles[nrow(estAngles),4]
  Xr[estAngles[,1]] <- estAngles[,4]

  Yr<-rep(NA,nrow(data))
  Yr[1]<-estAngles[1,5]
  Yr[nrow(data)]<-estAngles[nrow(estAngles),5]
  Yr[estAngles[,1]] <- estAngles[,5]

  Zr<-rep(NA,nrow(data))
  Zr[1]<-estAngles[1,6]
  Zr[nrow(data)]<-estAngles[nrow(estAngles),6]
  Zr[estAngles[,1]] <- estAngles[,6]

  Xr<-angleGapInterp(Xr,is.na(Xr))
  Yr<-angleGapInterp(Yr,is.na(Yr))
  Zr <-eventInterp(Zr,is.na(Zr))

  return(cbind(Xr,Yr,Zr))
}


pitchRoll<-function(xyz,degrees=FALSE){
  roll = atan(xyz[,2]/(xyz[,3]))
  pitch = atan((-xyz[,1])/((xyz[,2] * sin(roll)) + (xyz[,3]*cos(roll))))
  if(degrees==TRUE){
      pitch<- pitch*(180/pi)
      roll <- roll*(180/pi)
  }
  return(cbind(pitch,roll))
}



pitchRoll2<-function(xyz,degrees=FALSE){
  roll = atan2(xyz[,2],(xyz[,3]))
  pitch = atan2((-xyz[,1]),((xyz[,2] * sin(roll)) + (xyz[,3]*cos(roll))))
  if(degrees==TRUE){
    pitch<- pitch*(180/pi)
    roll <- roll*(180/pi)
  }
  return(cbind(pitch,roll))
}



RotData<-function(xyz,Xr,Yr,Zr){
  out<-apply(cbind(xyz,Xr,Yr,Zr),MARGIN=1,FUN = function(x){
    xyzTrans(xyz=c(x[1],x[2],x[3]),Xr=x[4],Yr=x[5],Zr=x[6])
  })
  return(t(out))
}


vectorAngle<-function(vec1,vec2){
  theta <- acos( sum(vec1*vec2) / (sqrt(sum(vec1 * vec1)) * sqrt(sum(vec2 * vec2)) ) )
  return(theta)
}


Gsep<-function(xyz,filt=rep(1,5*5)/(5*5)){
  X_Static<-filter(xyz[,1],filter = filt, sides = 2,circular = TRUE)
  Y_Static<-filter(xyz[,2],filter = filt, sides = 2,circular = TRUE)
  Z_Static<-filter(xyz[,3],filter = filt, sides = 2,circular = TRUE)

  X_Dynamic<-xyz[,1] - X_Static
  Y_Dynamic<-xyz[,2] - Y_Static
  Z_Dynamic<-xyz[,3] - Z_Static

  ODBA<- (abs(X_Dynamic) + abs(Y_Dynamic) + abs(Z_Dynamic))

  return(cbind(X_Static,Y_Static,Z_Static,X_Dynamic,Y_Dynamic,Z_Dynamic,ODBA))
}



tag<-function(x, y ,ang=0,exp=1,col="orange"){

  xlocs<-c(-0.25,-.25,.25,.25)

  ylocs<- c(-.25,.25,.25,-.25)

  shapePlot(x=x, y=y, shapeX=xlocs, shapeY=ylocs, ang=ang, exp=exp, col=col)
}


sharkback<-function(x, y, ang=0,exp=1,col="lightblue4"){

  ylocs<-c(1.002,  0.988,  0.920,  0.857,  0.752,  0.643,  0.594,  0.494,
           0.462, 0.430, 0.371, 0.321, 0.258, 0.172, 0.117, 0.131, 0.172, 0.185,
            0.199,  0.117,  0.018, -0.159, -0.282, -0.363,
           -0.445, -0.540, -0.604, -0.672, -0.708, -0.762, -0.817, -0.894,
           -0.966, -1.05, -0.966, -0.894, -0.817, -0.762, -0.708, -0.672, -0.604,
           -0.540, -0.445, -0.363, -0.282, -0.159,  0.018,  0.117,  0.199,
           0.185, 0.172, 0.131, 0.117, 0.172, 0.258, 0.321, 0.371, 0.430, 0.462,
           0.494, 0.594,  0.643,  0.752,  0.857,  0.920, 0.988,  1.002)

  xlocs<-c(-0.002, -0.025, -0.051, -0.086, -0.103, -0.128, -0.135, -0.135,
           -0.133, -0.171, -0.258, -0.326, -0.378, -0.397, -0.371, -0.303, -0.203, -0.138,
           -0.132, -0.132, -0.128, -0.122, -0.116, -0.112,
           -0.112, -0.116, -0.109, -0.070, -0.051, -0.015, -0.012, -0.012,
           -0.006, 0, 0.006,  0.012,  0.012,  0.015,  0.051,  0.070,  0.109,
           0.116,  0.112,  0.112,  0.116,  0.122,  0.128,  0.132,  0.132,
           0.138,  0.203,  0.303,  0.371,  0.397,  0.378, 0.326,  0.258, 0.171,  0.133,
           0.135,  0.135,  0.128,  0.103,  0.086,  0.051, 0.025,  0.002)
  shapePlot(x=x, y=y, shapeX=xlocs, shapeY=ylocs, ang=ang, exp=exp, col=col)


  ylocs<-c(0.24,  0.122,  0.058, -0.06, -0.06,  0.058,  0.122,  0.24)
  xlocs<-c(-0.005, -0.031, -0.027, -0.004,  0.004,  0.027,  0.031,  0.005)

  shapePlot(x=x, y=y, shapeX=xlocs, shapeY=ylocs, ang=ang, exp=exp, col=col)
}


shark<-function(x, y, ang=0,exp=1,col="lightblue4"){
  xlocs<-c(-0.540, -0.695, -0.705, -0.645, -0.750, -0.745, -0.535,
           -0.075, -0.105, -0.100, 0.075, 0.295, 0.470, 0.575,
           0.605, 0.615, 0.535, 0.345, 0.015, -0.175)

  ylocs<- c(-0.030, -0.270, -0.225, -0.005, 0.315, 0.385, 0.030,
            0.125, 0.245, 0.275, 0.130, 0.115, 0.090, 0.040, 0.005,
            -0.030, -0.070, -0.110, -0.125, -0.105)

  shapePlot(x=x, y=y, shapeX=xlocs, shapeY=ylocs, ang=ang, exp=exp, col=col)
}


sharkhead<-function(x,y,ang=0,exp=1,col="lightblue4"){
  ylocs<-c(-0.27705, -0.26344, -0.23170, -0.17276, -0.18636,
           -0.22263, -0.26344, -0.33, -0.23624, -0.16822,
           -0.11834, -0.07753,  0.03582,  0.15371,  0.23986,
           0.30788,  0.33509,  0.45751,  0.62982,  0.78852,
           0.87014,  0.87014,  0.78852,  0.62982,  0.45751,
           0.33509,  0.30788,  0.23986,  0.15371,  0.03582,
           -0.07753, -0.11834, -0.16822, -0.23624, -0.33,
           -0.26344, -0.22263, -0.18636, -0.17276, -0.23170,
           -0.26344, -0.27705)

  xlocs<-c(0.001213, -0.12418, -0.269817, -0.346677, -0.423537,
           -0.528713, -0.625799, -0.75, -0.641980, -0.544894,
           -0.431627, -0.358813, -0.338587, -0.285998, -0.217229,
           -0.124188, -0.087781, -0.071600, -0.047329, -0.027102,
           0.001213, -0.001213,  0.027102,  0.047329,  0.071600,
           0.087781,  0.124188,  0.217229,  0.285998,  0.338587,
           0.358813,  0.431627,  0.544894,  0.641980,  0.75,
           0.625799,  0.528713,  0.423537,  0.346677,  0.269817,
           0.124188, -0.001213)
  shapePlot(x=x, y=y, shapeX=xlocs, shapeY=ylocs, ang=ang, exp=exp, col=col)
}


shapePlot<-function(x, y, shapeX=NULL,shapeY=NULL,ang=0,exp=1,col="lightblue4"){
  op<-par()
  xdif<-op$usr[2]-par()$usr[1]
  ydif<-op$usr[4]-par()$usr[3]

  windasp<-op$pin[1]/op$pin[2]


  ang<-(ang/180)*pi
  xy<-Zb(ang) %*% rbind(shapeX,shapeY,rep(0,length(x)))

  xlocs<-(xy[1,]*exp*(xdif/2)) + x
  ylocs<-(xy[2,]*exp*windasp*(ydif/2)) + y

  polygon(x=xlocs,y=ylocs, col=col)

}


collapse<-function(dat,freq){
  d<-stats::filter(dat,filter = rep(1,freq)/freq,sides=2,circular=TRUE)
  d<-d[seq(1,length(dat), by=freq)]
  return(d)
}

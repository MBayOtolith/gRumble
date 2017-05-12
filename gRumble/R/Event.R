event<-function(x,ends=0, nas=TRUE,duration = c(NA,NA)){
  NAs<-is.na(x)
  if(sum(NAs)==length(x)){
    stop("All input values are NA")
  }
  if(any(NAs==TRUE)){
    NAloc<-which(NAs==TRUE)
    x[NAloc]<-nas
  }
  last<-length(x)

  if(ends==0){
    firstv<-x[1]
    lastv<-x[last]
    first<-1
    last<-length(x)
    if(firstv==TRUE){stop("First Value cannot be true, try ends = 1 or 2 ")
    }else if(lastv==TRUE){stop("Last Value cannot be true,try ends = 1 or 2")
    }
  }else if(ends==2){
    x<-c(FALSE,x,FALSE)
    first<-1
    last<-length(x)
  }else{
    falses<-which(x==FALSE)
    first<-falses[1]
    last<-falses[length(falses)]
  }
  
  edge<-which(((x[first:(last-1)]==TRUE)+(x[(first+1):last]==TRUE))==1)
  edge<-edge+first-1
  if(length(edge)==0){stop("There are no events where this occurs")}
  start<-edge[seq(1,length(edge),by=2)]
  end<-edge[seq(2,length(edge),by=2)]
  start<-start+1
  
  if(ends==2){
    start<-start-1
    end<-end-1
  }
  
  if(any(!is.na(duration))){
    dur<-end-start
    minD<-min(dur)
    maxD<-max(dur)
    if(!is.na(duration[1])){
      minD<-duration[1]
    }
    if(!is.na(duration[2])){
      maxD<-duration[2]
    }
    loc<-which(dur>=min & dur<=max)
    start<-start[loc]
    end<-end[loc]
  }
  return(cbind(start,end))
}




eventInterp<-function(dat,x=NULL,events=NULL,ends="none",na.rm=FALSE,nas=FALSE,duration=c(NA,NA)){
  if(is.numeric(dat)==FALSE){stop("data must be numeric to be interpolated")}
  
  if(!is.null(x)){
    loc<-event(x,ends=ends,nas=nas,duration=duration)
  }else{
    loc<-events
  }
  
  leng<-loc[,2] - loc[,1]+1
  for(i in 1 : nrow(loc)){
    if(loc[i,1]==1){
      dat[loc[i,1]:loc[i,2]]<-(dat[loc[i,2]+1])
    }else if( loc[i,2]==length(dat)){
      dat[loc[i,1]:loc[i,2]]<-(dat[loc[i,1]-1])
    }else{
      vals<-seq((dat[loc[i,1]-1]),(dat[loc[i,2]+1]), length.out=(leng[i]+2))
      dat[loc[i,1]:loc[i,2]]<-vals[2:(length(vals)-1)]
    }
  }
  return(dat)
}



eventFunc<-function(dat,events,fun=function(x){mean(x,na.rm=TRUE)}){
  out<-apply(events,1,FUN = function(rows){
    sub<-dat[rows[1]:rows[2]]
    return(fun(sub))
  })
  return(out)
}



length(smooth1)

PSM.plot <- function(Smooth, Data, indiv = NULL, type = c('Xs','YsObs')) {
  stateplot <- function(type,obs,len) {
    for(i in 1:len) {
      if(type!='Yp') {
        plot(Smooth[[j]]$Time,Smooth[[j]][[type]][i,],type="l",xlab="",ylab="")
      } else {
        plot(Data[[j]]$Time,na.omit(Smooth[[j]][[type]][i,]),type="l",xlab="",ylab="")
      }
      if(obs) {
        points(Data[[j]]$Time,Data[[j]]$Y[i,])
        rug(Data[[j]]$Time)
      }
      if(j==indiv[1])
        title(ylab=paste(type,i,sep=""))
      if(i==1 && k==1)
        title(main=paste('Individual',j))
    }
  }
  resfun <- function(type) {
    subs <- ceiling(length(Smooth[[j]]$Time)/length(Data[[j]]$Time)) - 1
    idx <- (1:dimT)*(subs+1)-subs
    Data[[j]]$Y-Smooth[[j]][[type]][,idx]
  }
  resplot <- function(type) {
    res <- resfun(type)
    for(i in 1:dimY) {
      plot(Data[[j]]$Time,res[i,],xlab="",ylab="")
      abline(h=0)
      if(j==indiv[1])
        title(ylab=paste('Y',i," - ",type,i,sep=""))      
      if(i==1 && k==1)
        title(main=paste('Individual',j))
    }
  }
  acfplot <- function(type) {
    res <- resfun(type)
    for(i in 1:dimY) {
      tmpylab <- ifelse(j==indiv[1],paste('ACF (Y',i," - ",type,i,')',sep=""),'')
      acf(res[i,],main="",ylab=tmpylab)
      if(i==1 && k==1)
        title(main=paste('Individual',j))
    }
  }
  dimS <- length(Smooth)
  dimX <- dim(Smooth[[1]]$Xs)[1]
  dimY <- ifelse(is.null(Data),0,dim(Data[[1]]$Y)[1])
  if(is.null(indiv))
    indiv = 1:dimS
  numrows <- 0
  for(k in 1:length(type)) {
    numrows <- numrows + switch(type[k], Xp = dimX, Xf = dimX, Xs = dimX, dimY)
  }
  plot.new()
  par(mfcol=c(numrows,length(indiv)),mar=c(2,4,2,0)+.1)
  for(j in indiv) {
    dimT <- dim(Data[[j]]$Y)[2]
    for (k in 1:length(type)) {
      switch(type[k],
             Xp = stateplot(type[k],obs=FALSE,len=dimX),
             Xf = stateplot(type[k],obs=FALSE,len=dimX),
             Xs = stateplot(type[k],obs=FALSE,len=dimX),
             Yp = stateplot(type[k],obs=FALSE,len=dimY),
             Ys = stateplot(type[k],obs=FALSE,len=dimY),
             YpObs = stateplot('Yp',obs=TRUE,len=dimY),
             YsObs = stateplot('Ys',obs=TRUE,len=dimY),
             res.p = resplot('Yp'),
             acf.p = acfplot('Yp'),
             res.s = resplot('Ys'),
             acf.s = acfplot('Ys')
             )       
    }
  }
}

PSM.plot(smooth1,Cpep,indiv=4:7, type = c('Xs','Ys','YpObs','res.p','acf.p'))
PSM.plot(smooth1,Cpep,indiv=4:5, type = c('YpObs','res.p'))
PSM.plot(smooth1,Cpep)



par(mfcol=c(3,2))
for(j in 2:3)
  for(i in 1:3) {
    plot(smooth1[[j]]$Time,smooth1[[j]]$Xs[i,],type="l",
         ylab=paste('state',i), xlab=paste('individual',j))
    if(i==1) points(Cpep[[i]]$Time,Cpep[[j]]$Y)
    rug(Cpep[[i]]$Time)
  }

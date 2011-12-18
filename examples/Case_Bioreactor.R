# Test Case
# Originating from CTSM
# Non-Linear Case

rm(list=ls())

library(PSM)
#library(nlme)

# Load the Data and Variables
tmpData <- read.table("sde0_1.csv",sep=";", col.names=c("Time","F","y1","y2","y3"))
tmpData <- read.table("~/PSM/PSM/examples/sde0_1.csv",sep=";", col.names=c("Time","F","y1","y2","y3"))
Time=tmpData$Time
Y=t(tmpData[,c("y1","y2","y3")])
U=t(as.matrix(tmpData[,c("F")]))
Data <- list(Time=Time,Y=Y,U=U)

MyModel <- list()
MyModel$ModelPar = function(THETA) {
  list(theta=list(
         mumax=THETA['mumax'],
         k1   =THETA['k1'],
         k2   =0.5,
         Y    =0.5,
         sf   =10,
         sig11=9.6704E-28,
         sig22=1.7471E-06,
         sig33=1.0903E-08,
         s11  =THETA['s11'],
         s22  =THETA['s22'],
         s33  =THETA['s33'],
         x0   =1.0095E+00,
         s0   =2.3835E-01,
         v0   =1.0040E+00))
}
MyModel$X0 = function(Time,phi,U) {
  matrix(c(phi$x0,phi$s0,phi$v0),ncol=1)
}
MyModel$h = function(eta,theta,covar) {
  phi = theta
}
MyModel$Functions$f = function(x,u,time,phi) {
  X <- x[1]; S <- x[2]; V <- x[3]; F <- u[1]
  mu <- phi$mumax*S/(phi$k2*S^2+S+phi$k1)
  matrix(c(
           mu*X-F*X/V,
           -mu*X/phi$Y+F*(phi$sf-S)/V,
           F
           ),ncol=1)
}
MyModel$Functions$df = function(x,u,time,phi) {
  X <- x[1]; S <- x[2]; V <- x[3]; F <- u[1]
  kssk = (phi$k2*S^2+S+phi$k1);
  matrix(c(
           phi$mumax*S/kssk-F/V,
           phi$mumax/kssk*X-phi$mumax*S/(kssk)^2*X*(2*phi$k2*S+1),
           F*X/V^2,
           -phi$mumax*S/(kssk)/phi$Y,
           -phi$mumax/(kssk)*X/phi$Y+phi$mumax*S/(kssk)^2*X/phi$Y*(2*phi$k2*S+1)-F/V,
           -F*(phi$sf-S)/V^2,
           0,
           0,
           0
           ),nrow=3,ncol=3,byrow=TRUE)
}
MyModel$Functions$g = function(x,u,time,phi) {
  x
}
MyModel$Functions$dg = function(x,u,time,phi) {
  diag(3)
}
MyModel$S = function(u,time,phi) {
  diag(c(phi$s11,phi$s22,phi$s33))
}
MyModel$SIG = function(u,time,phi) {
  diag(c(phi$sig11,phi$sig22,phi$sig33))
}
names(MyModel)

ctsmTHETA <- c(mumax=1.0022E+00,k1=3.1629E-02,
             s11=7.5248E-03,s22=1.0636E-03,s33=1.1388E-02)

phi <- MyModel$ModelPar(ctsmTHETA)$theta
ExtKalmanFilter(phi,MyModel,Data) #CTSM: -388.4857 #PSM: -388.4689


if(FALSE) { #Numerical gradients
  MyModel$Functions$df = function(x,u,time,phi) {
    jacobian(MyModel$Functions$f,x=x,u=u,time=time,phi=phi)
  }
  MyModel$Functions$dg = function(x,u,time,phi) {
    jacobian(MyModel$Functions$g,x=x,u=u,time=time,phi=phi)
  }
}

newData <- list()
for (i in 1:9) {
  newData[[i]] <- list(Time=seq(0,3.6,.9), U=matrix(rnorm(5)^2,nrow=1))
}

sim <- PSM.simulate(MyModel, newData, ctsmTHETA,deltaTime=0.01)
sim[[2]] <- Data
PSM.plot(sim,type=c('Y','U','longX'))

MyPar <- list(LB = 0.5*ctsmTHETA,
              Init = ctsmTHETA*1,
              UB = 1.5*ctsmTHETA)

system.time(
fit <- PSM.estimate(MyModel,sim,MyPar,CI=FALSE,trace=1)
            )
fit[1:5]

system.time(
fit <- PSM.estimate(MyModel,list(Data),MyPar,CI=FALSE,trace=1)
            )
fit[1:5]

###########
# Parameter uncertainty comparison
ctsm <- matrix(0,nrow=2,ncol=length(ctsmTHETA))
rownames(ctsm) <- c('MLE','SD')
colnames(ctsm) <- Pname <- names(ctsmTHETA)
ctsm['MLE',] <- ctsmTHETA
ctsm['SD',] <- c(3.3419E-03,1.5977E-03,9.7618E-04,1.3763E-04,1.4960E-03)
ctsm


par(mfrow=c(5,1),mar=c(2,4,2,2)+.1)
for(i in 1:length(Pname)) {
  psmCI <- fit$CI[c(1,3),Pname[i]]
  ctsmCI <- ctsm['MLE',Pname[i]]+ctsm['SD',Pname[i]]*1.96*c(-1,1)
  plot.new()
  plot.window(ylim=c(0.5,3),xlim=range(c(psmCI,ctsmCI)))
  axis(2,1:2,c('PSM','CTSM'),tick=FALSE,las=1,line=-1); axis(1);
  #box(bty='l')
  title(main=Pname[i],line=-.5)
  # CI limits
  lines(rep(psmCI,each=3),1+c(-.2,.2,0,0,.2,-.2),col=2,lwd=3)
  lines(rep(ctsmCI,each=3),2+c(-.2,.2,0,0,.2,-.2),col=4,lwd=3)
  # MLE
  points(fit$CI[2,Pname[i]],1,pch=23,bg=2,col=2,cex=1.5)
  points(ctsm['MLE',Pname[i]],2,pch=23,bg=4,col=4,cex=1.5)
}





nameY <- paste('y',1:3,sep='')
nameXp <- paste('Xp',1:3,sep='')
nameXf <- paste('Xf',1:3,sep='')
nameXs <- paste('Xs',1:3,sep='')
CTSMPred = read.table("BioReactorCTSMPred.csv",sep=";",
  col.names=c("Time",nameXp,paste('SD',nameXp),nameY, paste('SD',nameY)))
CTSMFilter = read.table("BioReactorCTSMFilter.csv",sep=";",
  col.names=c("Time",nameXf,paste('SD',nameXf)))
CTSMSmooth = read.table("BioReactorCTSMSmooth.csv",sep=";",
  col.names=c("Time",nameXs,paste('SD',nameXs)))
CTSMPred[1:2,];CTSMFilter[1:2,];CTSMSmooth[1:2,]


sm <- PSM.smooth(MyModel,list(Data),ctsmTHETA)
sm1 <- sm[[1]]


#save(fit,sm1,file="Case_Bioreactor.RData")

par(mfrow=c(3,1))
for(i in 1:3) {
  plot(sm1$Time,sm1$Xp[i,],type="l",col=2,xlab='Time',ylab=nameXp[i])
  title(main=paste('State prediction est.',nameXp[i],' - PSM vs. CTSM'))
  lines(CTSMPred[,'Time'],CTSMPred[,i+1],lty=2,lwd=2)
}


par(mfrow=c(3,1))
for(i in 1:3) {
  plot(sm1$Time,sm1$Xf[i,],type="l",col=2,xlab='Time',ylab=nameXf[i])
  title(main=paste('State filter est.',nameXp[i],' - PSM vs. CTSM'))
  lines(CTSMFilter[,'Time'],CTSMFilter[,i+1],lty=2,lwd=2)
}

par(mfrow=c(3,1))
for(i in 1:3) {
  plot(sm1$Time,sm1$Xs[i,],type="l",col=2,xlab='Time',ylab=nameXs[i])
  title(main=paste('State smoother est.',nameXp[i],' - PSM vs. CTSM'))
  lines(CTSMSmooth[,'Time'],CTSMSmooth[,i+1],lty=2,lwd=2)
}


#CALCULATE DEVIATION IN PCT
(pct <- max(abs(t(CTSMSmooth[,2:4])-sm1$Xs)/sm1$Xs*100))


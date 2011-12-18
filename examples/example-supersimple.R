library(PSM)


mod1 <- vector(mode="list")
mod1$Matrices=function(phi) {
  list(
       matA=matrix(c(-1*phi$k),ncol=1),
       matC=matrix(c(1/phi$V),nrow=1)
       )
}
mod1$X0 = function(Time,phi,U) {
  A0 <- phi$A0
  matrix(c(A0),ncol=1)
}
mod1$S = function(phi) {
  matrix(phi[["S"]])
}
mod1$ModelPar = function(THETA) {
  S <- 20
  A0 <- 1500
  list(theta=list(k=THETA['k'],V=THETA['V'], A0 = A0, S = S),
       OMEGA=matrix(.2))
}
mod1$SIG = function(phi) {
  matrix(50)
}
mod1$h = function(eta,theta,covar) {
  phi = theta
  phi$V = theta$V*exp(eta[1])
  phi
}



names(mod1)
TimeVec <- c(0,1,2,3,5,7,10,13,16,20,25)
DL = list(Time=10,State=1,Amount=500 )
PrepData = list( list(Time = TimeVec, Dose=DL),list(Time = TimeVec, Dose=DL) )

MyTHETA = c(k = 0.08, V = 15)

SimData <- PSM.simulate(mod1,PrepData,MyTHETA,deltaTime = .1)
SimData[[1]]$Y[c(1,6,7,8)] <- NA #MISSING OBS!
SimData[[2]]$Y[c(2,3,4)] <- NA #MISSING OBS!

names(SimData[[1]])
plot(SimData[[1]]$Time,SimData[[1]]$Y,ylim=c(10,150))
lines(SimData[[1]]$longTime,SimData[[1]]$longX/MyTHETA['V']/exp(SimData[[1]]$eta))
points(SimData[[2]]$Time,SimData[[2]]$Y,col=2)
lines(SimData[[2]]$longTime,SimData[[2]]$longX/MyTHETA['V']/exp(SimData[[2]]$eta),col=2)


par <- list(LB = MyTHETA/20,
            Init = c(k=.005,V=100),
            UB = MyTHETA*20)

con <- list(optimizer='ucminf',trace=T)
fit <- PSM.estimate(mod1,SimData,par,CI=TRUE,control=con)
fit[1:3]

sm <- PSM.smooth(mod1,SimData,fit$THETA,subs=50)
names(sm[[1]])
lines(sm[[1]]$Time,sm[[1]]$Ys,col=3)



##########################################################
### Non-linear model


mod2 <- mod1[c("X0","h","ModelPar")]
mod2$Functions <- 
  list(
       f = function(x,u,time,phi) {
         -1*phi$k*x
       },
       df = function(x,u,time,phi) {
         -1*phi$k
       },
       g = function(x,u,time,phi) {
         x/phi$V
       },
       dg = function(x,u,time,phi) {
         1/phi$V
       }
       )
mod2$S = function(u,time,phi) {
  as.matrix.default(phi[["S"]])
  #phi[["S"]]
}
mod2$SIG = function(u,time,phi) {
  as.matrix.default(50)
  #10
}
MyPar <- mod2$ModelPar(MyTHETA)
myphi <- mod2$h(0,MyPar$theta,NA)


source(file="~/PSM/PSM/R/ExtKalmanFilter.R")
source(file="~/PSM/PSM/R/PSM.estimate.R")
source(file="~/PSM/PSM/R/ModelCheck.R")
source(file="~/PSM/PSM/R/APL.KF.R")
source(file="~/PSM/PSM/R/APL.KF.gr.R")
source(file="~/PSM/PSM/R/IndividualLL.KF.R")
source(file="~/PSM/PSM/R/APL.KF.individualloop.R")
source(file="~/PSM/PSM/R/IndividualLL.KF.gr.R")
source(file="~/PSM/PSM/R/PSM.smooth.R")
source(file="~/PSM/PSM/R/PSM.simulate.R")
source(file="~/PSM/PSM/R/LinKalmanSmoother.R")
source(file="~/PSM/PSM/R/ExtKalmanSmoother.R")

library(compiler)
help(cmpfun)
setwd("~/PSM/PSM/R/")

cmpfile("ExtKalmanFilter.R")
loadcmp("ExtKalmanFilter.Rc")
source("ExtKalmanFilter.Rc")

#Rprof(tmp <- tempfile())
system.time(
            for(i in 1:1000) ExtKalmanFilter( myphi, mod2, SimData[[1]] )
            )
#Rprof()
#summaryRprof(tmp)
ekf2 <- cmpfun(ExtKalmanFilter)

#Rprof(tmp <- tempfile())
system.time(
            for(i in 1:1000) ekf2( myphi, mod2, SimData[[1]] )
            )

mod2c <- mod2
mod2c$Functions$f <- cmpfun(mod2$Functions$f)
mod2c$Functions$df <- cmpfun(mod2$Functions$df)
mod2c$Functions$g <- cmpfun(mod2$Functions$g)
mod2c$Functions$dg <- cmpfun(mod2$Functions$dg)
mod2c$SIG <- cmpfun(mod2$SIG)
mod2c$S <- cmpfun(mod2$S)
mod2c$h <- cmpfun(mod2$h)
mod2c$ModelPar <- cmpfun(mod2$ModelPar)

system.time(
            for(i in 1:1000) ekf2( myphi, mod2c, SimData[[1]] )
            )

#Rprof()
#summaryRprof(tmp)

ExtKalmanFilter( myphi, mod2, SimData[[1]] )


#Rprof()
#summaryRprof(tmp)
source(file="~/PSM/PSM/R/LinKalmanFilter.R")

LinKalmanFilter( myphi, mod1, SimData[[1]])


sm1 <- PSM.smooth(mod1,SimData,fit$THETA,subs=20)
sm2 <- PSM.smooth(mod2,SimData,fit$THETA,subs=20)
lines(sm2[[1]]$Time,sm2[[1]]$Ys,lty=2,lwd=2)


#as.vector(sm1[[1]]$Xs-sm2[[1]]$Xs)/as.vector(sm2[[1]]$Xs)*100

mod1.no.omega <- mod1
mod2.no.omega <- mod2
mod1.no.omega$ModelPar = function(THETA) {
  list(theta=list(k=THETA['k'],V=THETA['V'], A0 = 1500, S = 20), OMEGA=NULL)
}
mod2.no.omega$ModelPar <- mod1.no.omega$ModelPar

# Linear, no-OMEGA
fit1.no.omega <- PSM.estimate(mod1.no.omega,SimData,par)
fit1.no.omega[1:2]

# Non-Linear, no-OMEGA
fit2.no.omega <- PSM.estimate(mod2.no.omega,SimData,par)
fit2.no.omega[1:2]


sm <- PSM.smooth(mod1.no.omega,SimData,fit$THETA,subs=50)
sm2 <- PSM.smooth(mod2.no.omega,SimData,fit$THETA,subs=50)
lines(sm[[1]]$Time,sm[[1]]$Ys,col=5)
lines(sm2[[1]]$Time,sm2[[1]]$Ys,lty=2)


# Linear, with-OMEGA
fit1 <- PSM.estimate(mod1,SimData,par)
fit1[1:2]

# Non-Linear, with-OMEGA
#Rprof(tmp <- tempfile())
fit2 <- PSM.estimate(mod2,SimData,par,trace=1)
fit2[1:2]
fit2$sec
#Rprof()
#summaryRprof(tmp)


fit2c <- PSM.estimate(mod2c,SimData,par,trace=1)
fit2c[1:2]
fit2c$sec

#APL.KF.individualloop <- cmpfun(APL.KF.individualloop)
APL.KF <- cmpfun(APL.KF)
APL.KF.gr <- cmpfun(APL.KF.gr)
CutThirdDim <- cmpfun(CutThirdDim)
ExtKalmanFilter <- cmpfun(ExtKalmanFilter)
IndividualLL.KF <- cmpfun(IndividualLL.KF)
IndividualLL.KF.gr <- cmpfun(IndividualLL.KF.gr)
ModelCheck <- cmpfun(ModelCheck)
PSM.estimate <- cmpfun(PSM.estimate)
invlogit <- cmpfun(invlogit)
logit <- cmpfun(logit)
rm(list=c('APL.KF','APL.KF.gr','APL.KF.individualloop','CutThirdDim',
     'ExtKalmanFilter','IndividualLL.KF','IndividualLL.KF.gr',
     'ModelCheck','PSM.estimate','invlogit','logit'))
fit2cc <- PSM.estimate(mod2c,SimData,par,trace=1)
fit2cc[1:2]
fit2cc$sec
(fit2$sec-fit2cc$sec)/fit2$sec*100

source(file="~/PSM/PSM/R/ExtKalmanSmoother.R")
source(file="~/PSM/PSM/R/ExtKalmanFilter.R")

sm1 <- PSM.smooth(mod1,SimData,fit$THETA,subs=20)
sm2 <- PSM.smooth(mod2,SimData,fit$THETA,subs=50)
lines(sm2[[1]]$Time,sm2[[1]]$Ys,lty=2,lwd=2,col=2)


V1 <- fit$THETA['V']*exp(sm2[[1]]$eta)
plot(sm2[[1]]$Time,sm2[[1]]$Xs,col=3,type="l",lwd=2)
st <- sqrt(sm2[[1]]$Ps[1,1,])
lines(sm2[[1]]$Time,sm2[[1]]$Xs+1.96*st,col=3,type="l")
lines(sm2[[1]]$Time,sm2[[1]]$Xs-1.96*st,col=3,type="l")


##### Penalty function.

if(FALSE) {
  mini = .1
  maxi = 10
  xvec = seq(mini,maxi,length=1000)
  lambda = 1e-4
  dx = 1e-30
  dx=0
  p = lambda*(mini/(xvec-mini+dx) + maxi/(maxi-xvec+dx))
  range(p)
  plot(xvec,p,type="l",ylim=c(0,.001))
}


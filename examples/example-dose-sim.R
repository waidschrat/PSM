rm(list=ls())

detach(package:PSM)
library(PSM)

REDO = FALSE

# This example demonstrates the simulation of a bolus dose
# in a 2-compartment system.
V1 <- 5      #[L]
V2 <- 10     #[L]
CLd <- 0.005 #[L/min]
CL <- 0.002  #[L/min]

Model.SimDose <- list(
                      Matrices = function(phi) {
                        V1i <- phi[["V1i"]]
                        matA <- matrix(c(-(CL+CLd)*V1i ,  CLd*V2 ,   
                                           CLd*V1i , -CLd*V2 ) ,nrow=2,byrow=T)
                        matC <- matrix(c(1/V1i,0),nrow=1)
                        list(matA=matA,matC=matC)
                      },
                      X0 = function(Time=Na,phi,U=Na) {
                        matrix(0,nrow=2)
                      },
                      SIG = function(phi) { 
                        sig11 <- phi[["sig11"]]
                        matrix(c( sig11,0,
                                 -sig11,0), nrow=2, byrow=T)
                      },
                      S = function(phi) {
                        matrix(phi[["S"]])
                      },
                      h = function(eta,theta,covar) {
                        phi <- theta
                        phi[["V1i"]] <- V1*exp(eta[1])
                        phi                      },
                      ModelPar = function(THETA){
                        list(theta=list(sig11=THETA[1], S=THETA[2]),
                             OMEGA=matrix(THETA[3]) )
                      }
                      )

# Create Simulation Timeline 
SimDose.Subj <- 2
if(REDO) {
  SimDose.Data <- vector(mode="list",length=SimDose.Subj)
  for (i in 1:SimDose.Subj) {
    SimDose.Data[[i]]$Time <- seq(from=10,by=10,to=400)
    SimDose.Data[[i]]$Dose <- list(
                        Time = c(30,180),
                        State = c(1, 1),
                        Amount = c(1500,1500)
                        )
  }
}

#############
# Simulation
#############

#                   sig11  S  OMEGA
SimDose.THETA <-  c(10 , 40 , .5)


if(REDO) {
  SimDose.Data <- PSM.simulate(Model.SimDose, SimDose.Data, SimDose.THETA, deltaTime=.1 , individuals=SimDose.Subj)
} else {
  load("simdose.RData")
}
  # View the Simulated datastructure
names(SimDose.Data[[1]])


PSM.plot(SimDose.Data,type=c('Y','longX','eta'))
# Plot of the simulations
par(mfcol=c(3,SimDose.Subj),mar = c(2, 4, 2, 2)+.1)
for(id in 1:SimDose.Subj) {
  plot(SimDose.Data[[id]]$Time , SimDose.Data[[id]]$Y,
         ylab="Observations", main=paste('individual',id,', eta:',round(SimDose.Data[[id]]$eta,3)))
  for(i in 1:2) {
    plot(SimDose.Data[[id]]$longTime , SimDose.Data[[id]]$longX[i,],type="l",
         ylab=paste('state',i))
    rug(SimDose.Data[[id]]$Time)
  }
  print(SimDose.Data[[id]]$eta)
}
var(c(SimDose.Data[[1]]$eta,SimDose.Data[[2]]$eta))


## TEST KALMAN FILTER
source("../R/LinKalmanFilter.R",echo=F)
theta <- Model.SimDose$ModelPar(SimDose.THETA)$theta
OMEGA <- Model.SimDose$ModelPar(SimDose.THETA)$OMEGA
indv = 2
DataI <- SimDose.Data[[indv]]
etaI =.5
phi <- Model.SimDose$h(eta=etaI,theta=theta)
of <- LinKalmanFilter(phi=phi, Model=Model.SimDose , Data=DataI , output=T,fast=TRUE)
os <- LinKalmanFilter(phi=phi, Model=Model.SimDose , Data=DataI , output=T,fast=FALSE)
names(of)
names(os)
plot(as.vector(DataI$Y*V1*exp(etaI)))
plot(as.vector(DataI$Y))
lines(as.vector(os$Yp))
lines(as.vector(of$Yp),lty=2,col=2)
plot(as.vector(of$Xf[1,])/as.vector(os$Xf[1,]))
of$negLogLike
os$negLogLike
system.time(
            for(i in 1:100)
            o <- LinKalmanFilter(phi=phi, Model=Model.SimDose , Data=DataI , output=T,fast=TRUE)
            )
o$negLogLike
par(mfrow=c(3,1))
for(i in 1:3)
  plot(DataI$T,o$Xf[i,],type="l")
rug(DataI$T)





###########
# Estimate
###########

parA <- list(LB=SimDose.THETA/50, Init=SimDose.THETA , UB=SimDose.THETA*50 )
fitA <- PSM.estimate(Model=Model.SimDose,Data=SimDose.Data,Par=parA,CI=T,trace=1)
fitA[1:3]

Model.SimDoseB <- Model.SimDose
Model.SimDoseB$ModelPar = function(THETA){
                        list(theta=list(sig11=0, S=THETA[1]), OMEGA=matrix(THETA[2]) ) }
parB <- list(LB=SimDose.THETA[2:3]/50, Init=SimDose.THETA[2:3] , UB=SimDose.THETA[2:3]*50 )

fitB <- PSM.estimate(Model=Model.SimDoseB,Data=SimDose.Data,Par=parB,CI=T,trace=1)
fitB[1:3]
#TEST for SDE
pchisq(2*(fitB$NegLogL-fitA$NegLogL), 3-2, lower.tail = FALSE) #p=0.005 - significant


#Plot individual LL with init param.
par(mfcol=c(3,2))
mp <- Model.SimDose$ModelPar(SimDose.THETA)
etavec <- seq(from=-2,to =2,length=50)
outvec <- vector(length=length(etavec))
gradvec <- vector(length=length(etavec))
for (id in 1:SimDose.Subj) {
  for (i in 1:length(etavec)) {
    outvec[i] <- IndividualLL.KF(etavec[i],mp$theta,mp$OMEGA,Model.SimDose,SimDose.Data[[id]])
    a <- IndividualLL.KF(etavec[i]+.001,mp$theta,mp$OMEGA,Model.SimDose,SimDose.Data[[id]])
    gradvec[i] <- (a-outvec[i])/.001
  }
  plot(etavec,outvec)
  abline(h=0)
}


###########
# Smoothing
###########
#APL.KF(THETA=fitA$THETA,Model=Model.SimDose,Pop.Data=SimDose.Data,GUIFlag=1,longOutput=TRUE)
#source("../R/PSM.smooth.R",echo=F)
out <- PSM.smooth(Model = Model.SimDose, Data = SimDose.Data, THETA = fitA$THETA, subsample = 20)
#save(list=c('SimDose.Data','fitA','out'),file='simdose.RData')
# View the data structure
names(out[[1]])
PSM.plot(SimDose.Data,out)
PSM.plot(Data=SimDose.Data,out,type=c('X','longX'))
1

#Plot of the smoothed estimates
par(mfcol=c(3,SimDose.Subj),mar = c(2, 4, 2, 2)+.1)
for(id in 1:SimDose.Subj) {
  plot(SimDose.Data[[id]]$Time , SimDose.Data[[id]]$Y,
         ylab="Observations", main=paste('Individual',id))
  legend(400, y = max(SimDose.Data[[id]]$Y), legend=c('Smooth est.'),lty=c('71A1'),col=c(4),
           xjust=1,yjust=1,bty="n")
  lines(out[[id]]$Time , out[[id]]$Xs[1,]/(V1*exp(out[[id]]$eta)),lty='71A1',col=4)
  for(i in 1:2) {
    plot(SimDose.Data[[id]]$longTime , SimDose.Data[[id]]$longX[i,],type="l",
         ylab=paste('Smooth state',i))
    lines(out[[id]]$Time , out[[id]]$Xs[i,],lty='71A1',col=4)
    rug(SimDose.Data[[id]]$Time)
    legend(400, y = max(out[[id]]$Xs[i,]), legend=c('Simulation','Smooth est.'),lty=c('solid','71A1'),col=c(1,4),
           xjust=1,yjust=1,bty="n",y.intersp=.8)
  }
}

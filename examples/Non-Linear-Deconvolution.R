# Exercise 3-2
# Deconvolution of Infusion rates
# 1 compartment system

rm ( list=ls() )
graphics.off()

library(PSM)


### Non-Linear Version
NL.mod <- vector(mode="list")
NL.mod$Functions <- list(
       f = function(x,u,time,phi) { 
       ka <- (phi["Emax"]*x[1]) / (phi["ED50"]+x[1])
       c(-ka*x[1] , ka*x[1]-phi["KE"]*x[2])
       },
       df = function(x,u,time,phi) {
         tmp <- -(phi["Emax"]*x[1]*2)/(phi["ED50"] + x[1])-(phi["Emax"]*x[1]^2)/(phi["ED50"]+x[1])^2
         matrix(c(tmp,-tmp,0,-phi["KE"]),nrow=2)  },
       g = function(x,u,time,phi) {  matrix( log(x[2]),nrow=1)},
       dg = function(x,u,time,phi) { matrix(c(0,1/x[2],nrow=1))       }
       )
NL.mod$S = function(u,time,phi) {  matrix(phi["S"],nrow=1)}
NL.mod$SIG = function(u,time,phi) {  matrix(c(phi["SIG"],-phi["SIG"],0,0),nrow=2)}
NL.mod$X0 = function(Time,phi,U) { matrix(c(1e-4,1e-4),ncol=1)}
NL.mod$ModelPar = function(THETA) {
  list(theta=c(THETA["Emax"], THETA["ED50"], THETA["KE"], THETA["S"], THETA["SIG"]),
       OMEGA=NULL)
}
NL.mod$h = function(eta,theta,covar) { phi = theta }
NL.mod$Dose = list(Time=0, State=1, Amount=500 )

simTHETA <- c( Emax=.01 , ED50=100, KE=.05, S=.01, SIG=1)

# Create Simulation Timeline 
NID   <- 1
Sim.Data <- vector("list",NID)
for(i in 1:NID) Sim.Data[[i]] <- list(Time = seq(from=0,by=2,to=200))
Sim[[1]]$Y[1:5]

# Simulated Data
NL.mod$ModelPar(simTHETA)$theta
Sim <- PSM.simulate(NL.mod, Sim.Data, simTHETA, deltaTime=.1 )

# Plot simulated profiles
YR <- c(0,max(unlist(lapply(Sim , function(x) max(exp(x$Y)) ))))
par(mfrow=c(3,1))
plot(Sim[[1]]$Time, exp(Sim[[1]]$Y), type="b", col=1,ylim=YR)
plot(Sim[[1]]$Time, (simTHETA['Emax']*Sim[[1]]$X[1,])/(100+Sim[[1]]$X[1,]) ,main="KA")
plot(Sim[[1]]$Time, (Sim[[1]]$X[1,]) ,main="Abs comp")

if(NID!=1)
  for(i in 2:NID) lines( Sim[[i]]$Time, exp(Sim[[i]]$Y), type="b", col=i)



# Subsample Data
Est <- vector("list",NID)
NoObs <- 7
for(i in 1:NID) {
  IDX <- unique(sort(c(1,sample(seq_along(Sim[[i]]$Time) , 7))))
  Est[[i]]$Time <- Sim[[i]]$Time[IDX]
  Est[[i]]$Y    <- matrix(Sim[[i]]$Y[IDX],nrow=1) }

# Plot subsampled data
YR <- c(0,max(unlist(lapply(Est , function(x) max(exp(x$Y)) ))))
plot(Est[[1]]$Time, exp(Est[[1]]$Y), type="b", col=1,xlim=c(0,200),ylim=YR)
if(1 != NID)
  for(i in 2:NID) lines( Est[[i]]$Time, exp(Est[[i]]$Y), type="b", col=i)
  
# Define Non-Linear Deconvolution Model
NL.deconv <- NL.mod
NL.deconv$Functions <- list(
       f = function(x,u,time,phi) { 
       matrix(c(-x[3]*x[1] , 
                x[3]*x[1]-phi["KE"]*x[2],
                    0 ),nrow=3)
       },
       df = function(x,u,time,phi) {
         matrix(c(-x[3],    0      , -x[1],
                   x[3], -phi["KE"],  x[1],
                   0   ,    0      ,   0 ),nrow=3,byrow=TRUE)  },
       g  = function(x,u,time,phi) { matrix(log(abs(x[2])),nrow=1)},
       dg = function(x,u,time,phi) { matrix(c(0,1/x[2],0),nrow=1)       }
)
NL.deconv$S     = function(u,time,phi) {  matrix(phi["S"],nrow=1)}
NL.deconv$SIG   = function(u,time,phi) {  matrix(c(rep(0,8),phi["SIG"]),nrow=3)}
NL.deconv$X0    = function(Time,phi,U) { matrix(c(1e-4,1e-4,phi["KA0"]),ncol=1)}
NL.deconv$ModelPar = function(THETA) {
  THETA["KE"] <- 0.05
  list(theta=c(THETA["KE"], THETA["S"], THETA["SIG"], THETA["KA0"]
         ),
       OMEGA=NULL)
}
simTHETA['Emax']*500/(100+500) #KA0
myTHETA <- c(S=.01,SIG=0.0001, KA0=.00833)
myPar   <- list(LB=.1*myTHETA, Init=myTHETA, UB=10*myTHETA)

NL.deconv$ModelPar(myTHETA)$theta

out <- PSM.estimate(Model=NL.deconv, Data=Sim, Par=myPar,trace=1)
#out <- PSM.estimate(Model=NL.deconv, Data=Est, Par=myPar,trace=1)
out[1:2]
#source("C:/SKLI/Code/R/R-Forge/PSM/R/ExtKalmanSmoother.R")

smo <- PSM.smooth(Model=NL.deconv, Data=Sim , THETA=myTHETA, subsample=0,trace=1)
#Est[[1]]$Y
#Sim[[1]]$Y

smo1 <- smo[[1]]
Sim1 <- Sim[[1]]
par(mfrow=c(3,1))
plot(smo1$Time,smo1$Xp[1,],type="l",ylim=c(0,500))
lines(Sim1$longTime,Sim1$longX[1,],col=2)

plot(smo1$Time,smo1$Xp[2,],type="l",ylim=c(0,80))
lines(Sim1$longTime,Sim1$longX[2,],col=2)
points(Sim1$Time,exp(Sim1$Y),col=2)

plot(smo1$Time,smo1$Xp[3,],type="l",ylim=range(smo1$Xp[3,]))
KAvec <- (simTHETA['Emax']*Sim[[1]]$X[1,])/(100+Sim[[1]]$X[1,])
lines(Sim[[1]]$Time, KAvec,col=3)

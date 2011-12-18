# Non-Linear Simulation file
rm(list=ls())

library(PSM)

###############################################################################
# Model 1:
#     Simple 1 compartment A is NON-Singular
#
# Model 2:
#     Simple 1 compartment with infusion and addition from Stochastic compartment
#     A singular
#
# Model 3:
#     Non-Linear Model with input.
###############################################################################

# Model 1
mod1 <- vector(mode="list")
mod1$Matrices=function(phi) {  list(
       matA=matrix(c(-1*phi$k),ncol=1),
       matC=matrix(c(1/phi$V),nrow=1)        ) }
mod1$X0 = function(Time,phi,U) {  matrix(c(phi$A0),ncol=1)}
mod1$S = function(phi) {  matrix(phi[["S"]])}
mod1$ModelPar = function(THETA) {  list(theta=list(k=THETA['k'],V=THETA['V'], A0 = 1500, S = 10),
OMEGA=NULL)}
mod1$SIG = function(phi) {  matrix(20)}
mod1$h = function(eta,theta,covar) { phi = theta ; phi$V = theta$V*exp(eta[1]) ; phi }
mod1$Dose = list(Time=10,State=1,Amount=500 )

TimeVec <- c(0,1,2,3,5,7,10,10.1,13,16,20,25)
PrepData = list( list(Time = TimeVec),list(Time = TimeVec) )
MyTHETA = c(k = 0.08, V = 15)

SimData1 <- PSM.simulate(mod1,PrepData,MyTHETA,deltaTime = .1)
plot(SimData1[[1]]$Time , SimData1[[1]]$Y , type="b",lwd=2)


# Model 2
mod2 <- vector(mode="list")
mod2$Matrices=function(phi) {  list(
       matA=matrix(c(-1*phi$k,0,1,0),nrow=2),
       matB=matrix(c(1,0),nrow=2),
       matC=matrix(c(1/phi$V,0),nrow=1),
       matD=matrix(0,nrow=1)      ) }
mod2$X0 = function(Time,phi,U) {  matrix(c(phi$A0,0),nrow=2)}
mod2$S = function(phi) {  matrix(phi[["S"]])}
mod2$SIG = function(phi) {  matrix(c(0,0,0,0),nrow=2)}
mod2$ModelPar = function(THETA) {  list(theta=list(k=THETA['k'],V=THETA['V'], A0 = 0, S = 10),
OMEGA=NULL)}
mod2$h = function(eta,theta,covar) { phi = theta ; phi$V = theta$V*exp(eta[1]) ; phi }
mod2$Dose = list(Time=10,State=1,Amount=500 )

TimeVec <- c(0,1,2,3,5,7,10,10.1,13,16,20,25)
U <- matrix( rep(0, length(TimeVec)) ,nrow=1)
U[TimeVec<=5] <- 50
PrepData = list( list(Time = TimeVec,U=U),list(Time = TimeVec,U=U) )
MyTHETA = c(k = 0.08, V = 15)

SimData2 <- PSM.simulate(mod2,PrepData,MyTHETA,deltaTime = .1)
par(mfrow=c(3,1))
plot(SimData2[[1]]$longTime , SimData2[[1]]$longX[1,] , type="b",lwd=2)
plot(SimData2[[1]]$longTime , SimData2[[1]]$longX[2,] , type="b",lwd=2)
plot(SimData2[[1]]$Time , SimData2[[1]]$Y , type="b",lwd=2)

# With scaling diffusion term
mod2$SIG = function(phi) {  matrix(c(0,0,0,10),nrow=2)}
SimData2 <- PSM.simulate(mod2,PrepData,MyTHETA,deltaTime = .1)
par(mfrow=c(3,1))
plot(SimData2[[1]]$longTime , SimData2[[1]]$longX[1,] , type="b",lwd=2)
plot(SimData2[[1]]$longTime , SimData2[[1]]$longX[2,] , type="b",lwd=2)
plot(SimData2[[1]]$Time , SimData2[[1]]$Y , type="b",lwd=2)





#Model 3:  Non-linear model
mod3 <- mod1[c("X0","h","ModelPar","Dose")]
mod3$Functions <- list(
       f = function(x,u,time,phi) { -1*phi$k*x},
       df = function(x,u,time,phi) {-1*phi$k},
       g = function(x,u,time,phi) {x/phi$V},
       dg = function(x,u,time,phi) {1/phi$V})
mod3$S = function(u,time,phi) {  as.matrix.default(phi[["S"]]) }
mod3$SIG = function(u,time,phi) {  as.matrix.default(20) }
MyPar <- mod2$ModelPar(MyTHETA)
myphi <- mod2$h(0,MyPar$theta,NA)
unlist(myphi)

SimData3 <- PSM.simulate(mod3,PrepData,MyTHETA,deltaTime = .1)

# Visual Inspection
par(mfrow=c(2,1),mar=c(4,4,1,1))
plot(   SimData1[[1]]$longTime ,SimData1[[1]]$longX,type="l")
lines(  SimData3[[1]]$longTime ,SimData3[[1]]$longX,col="red")
# MAthematical solution to 1-compartment model
Time <- SimData1[[1]]$longTime
DTime <- 10 ; Dose <- 500
matConc <- rep( 0 , length(Time))
matConc[Time<=DTime] <- 1500*exp(-MyTHETA["k"]*Time[Time<=DTime])
matConc[Time>DTime]  <- (matConc[Time==DTime]+Dose)*exp(-MyTHETA["k"]*(Time[Time>DTime]-DTime))
lines(Time,matConc,col="blue",lwd=1)

plot(   SimData1[[1]]$Time ,SimData1[[1]]$Y,type="b",ylim=range(SimData1[[1]]$Y,SimData3[[1]]$Y))
points( SimData3[[1]]$Time ,SimData3[[1]]$Y,type="b",col="red")





########################################################################
#### Deleted

# Testing multivariate rnorm
if(FALSE) {
library(MASS)
help(mvrnorm)
tmp  <- mvrnorm(n=1000, mu=c(2,5), Sigma=matrix(c(1,1.5,1.5,4),2))
plot(tmp,pch=".")
colMeans(tmp)
cov(tmp)
cor(tmp)
}



source("../R/PSM.simulate.R")
Rprof(tmp <- tempfile())
for(i in 1:100) PSM.simulate(mod2,PrepData,MyTHETA,deltaTime = .1)
Rprof()
summaryRprof(tmp)

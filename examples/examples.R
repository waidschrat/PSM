# This example file demonstrates the process of 
# Simulation, estimation and smoothing
# The general model used is the C-peptide model from
# Cauter et. al.


rm(list=ls())

detach(package:PSM)
library(PSM,lib.loc="~/PSM/Rpackages/gridterm")


#######################
# Model Intervention Model
#######################

k1 = 0.053; k2 = 0.051; ke = 0.062;

Model.Sim <- list(
               Matrices = function(phi) {
                 a1  <- phi[["a1"]]
                 a2  <- phi[["a2"]]
                 B  <- phi[["B"]]
                 K  <- phi[["K"]]
                 matA <- matrix( c(-(k1+ke) ,  k2 ,   1 ,   0,
                                         k1 , -k2 ,   0 ,   0,
                                          0 ,   0 , -a1 ,  a1,
                                          0 ,   0 ,   0 , -a2),nrow=4,byrow=T)                
                 matB <- matrix( c(0           , 0      ,
                                   0           , 0      ,
                                   B           , 0      ,
                                   0           , a2*K     ),byrow=T,nrow=4)
                 matC <- matrix(c(1,0,0,0),nrow=1)
                 matD <- matrix(c(0,0),nrow=1)
                 return(list(matA=matA,matB=matB,matC=matC,matD=matD))},
               X0 = function(Time=NA,phi,U=NA) {
                 C0 <- phi[["C0"]]
                 tmp    <- C0
                 tmp[2] <- C0*k1/k2
                 tmp[3] <- C0*ke
                 tmp[4] <- 0
                 return(matrix(tmp,ncol=1) )} ,
               SIG = function(phi) {
                 return( diag( c(0,0,phi[["SIG33"]],0)) ) } ,
               S = function(phi) {
                 return( matrix(phi[["S"]])) } ,
               h = function(eta,theta,covar) {
                 phi <- theta
                 phi[["B"]] <- theta[["B"]]*exp(eta[1])
                 phi[["K"]] <- theta[["K"]]*exp(eta[2])
                 phi[["C0"]] <- theta[["C0"]]*exp(eta[3])
                 return(phi) } ,
               ModelPar = function(THETA){
                 return(list(theta=list(C0=900,S=8500,
                               a1=THETA[1],a2=THETA[2],SIG33=THETA[3],
                               K = THETA[4], B = THETA[5]),
                             OMEGA=diag(c(.2,2,.2))))}
               )



# Create Simulation Timeline and Simulation Input
NoOfSubjects <- 2
Sim.Data <- vector(mode="list",length=NoOfSubjects)
for (i in 1:NoOfSubjects) {
  Sim.Data[[i]]$Time <- c( 0,15,30,45,60,75,90,120,150,180,210,240,270,300,330,360,420,480,600,615,630,645,660,675,690,720,750,780,810,840,960,1140,1320,1410,1440)
  Sim.Data[[i]]$U <- matrix(c( rep(1,35) , 
                    as.numeric( Sim.Data[[i]]$Time %in% c(0,15,240,600,615)) ),byrow=T,nrow=2)
}

Sim.Data



                                        # Create Simulation THETA parameter
# a1=THETA[1],a2=THETA[2],SIG33=THETA[3], K = THETA[4], B = THETA[5]),
(Sim.THETA <-  c(0.02798 , 0.01048 , 6.9861 , 427.63 , 1.7434))
# (Sim.THETA <-  c(0.02798 , 0.01048 , 0 , 427.63 , 1.7434))

Sim.Data <- PSM.simulate(Model.Sim, Sim.Data, Sim.THETA, dt=.1 ,individuals=NoOfSubjects)

Sim.Data

# Plot the Result
par(mfrow=c(2,2))
for(id in 1:NoOfSubjects) {
  for(i in 1:4) {
    plot(Sim.Data[[id]]$Time , Sim.Data[[id]]$X[i,],type="l",
         ylab=paste('state',i), xlab=paste('individual',id))
    rug(Sim.Data[[id]]$Time)
  }
}


#######################################################
### Estimation on Simulation Data
#######################################################

# Transform Simulation Data -> Estimation Data
Pop.Data <- vector( mode="list" , length=NoOfSubjects)
for(i in 1:NoOfSubjects) {
  Pop.Data[[i]] <- list(Time=Sim.Data[[i]]$Time,
                Y=matrix( Sim.Data[[i]]$Y , 1 , length(Sim.Data[[i]]$Y)),
                U=NULL)
}

Pop.Data

# Insulin Secretion Rates - Deconvolution.
k1 = 0.053; k2 = 0.051; ke = 0.062;

Model.Est <- list(
            Matrices=function(phi) { list(
              matA=matrix(c(-(k1+ke),  k2, 1,
                             k1     , -k2, 0,
                             0      ,   0, 0  ),ncol=3,byrow=T),
              matB=NA,
              matC=matrix(c(1,0,0),nrow=1),
              matD=NA )  },
            X0 = function(Time=NA,phi=NA,U=NA) {
              C0 <- phi[["C0"]]
              tmp    <- C0
              tmp[2] <- C0*k1/k2
              tmp[3] <- C0*ke
                    return(matrix(tmp,ncol=1) )} ,
            SIG = function(phi) {
                      return( diag( c(1e-3,1e-3,phi[["SIG33"]])) ) } ,
            S = function(phi) {
                      return( matrix(phi[["S"]])) } ,
            h = function(eta,theta,covar) {
              phi <- theta
              phi[["C0"]] <- theta[["C0"]]*exp(eta[1])
              return(phi) } ,
            ModelPar = function(THETA){
              return(list(theta=list(C0=THETA[1],S=THETA[2],SIG33=THETA[3]),
                          OMEGA=matrix(THETA[4])))}
            )


# -------------------------------------------------------------
# PSM - Minimization
# -------------------------------------------------------------

#                       C0     S    SIG   OMEGA
par1 <- list(LB   = c(  200,  50^2,   0,  .0 ),
             Init = c( 1000, 100^2,  10,  .25),
             UB   = c( 3000, 150^2,  15,  .50))



obj1 <- PSM.estimate(Model=Model.Est, Data=Pop.Data, Par=par1,CI=T,trace=2)

obj2 <- PSM.estimate(Model=Model.Est, Data=Pop.Data, Par=par1,CI=T,trace=2,optimizer="nlm")

obj3 <- PSM.estimate(Model=Model.Est, Data=Pop.Data, Par=par1,CI=T,trace=2,optimizer="blaa")

obj1

(THETA <- obj1$THETA)
# THETA <- c(1.885498e+03, 2.584864e+03, 9.862890e+00, 8.981350e-03)

# -------------------------------------------------------------
# Smoother 
# -------------------------------------------------------------

Data.Sm <- PSM.smooth( Model=Model.Est , Data=Pop.Data, THETA=THETA, subsample=5,trace=1)


ID <- NoOfSubjects
CMT <- 3
Data <- Data.Sm[[ID]]
plot( Data$Time, Data$Xs[CMT,] , type="n", ylim=c(0,250),
     xlab="Min",ylab="pmol/min",main="Insulin Secretion Rate")
polygon( c(Data$Time,rev(Data$Time)) ,
        c(Data$Xs[CMT,]+sqrt(abs(Data$Ps[CMT,CMT,])) ,
          rev(Data$Xs[CMT,]-sqrt(abs(Data$Ps[CMT,CMT,])))),col=4,density=50)
lines( Data$Time, Data$Xs[CMT,], type="l",lwd=2)
points( Sim.Data[[ID]]$Time, Sim.Data[[ID]]$X[3,] , col="red")
rug(Data$Time)


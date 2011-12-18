rm(list=ls())

PC <- TRUE

if (PC) {
  library(PSM)
  DataPath  <- "C:/Data/"
} else {
  detach(package:PSM); detach(package:mexp)
  library(PSM,lib.loc="~/PSM/Rpackages/gridterm")
  DataPath  <- "~/PSM/isr/"
}



###############################################
# Load C-peptide data for 12 individuals
###############################################

tmpData <- read.table(paste(DataPath,"CPEP",sep=""),sep="\t", col.names=paste('ID',1:13))
Cpep <- vector(mode="list")
for(i in 1:12)
  Cpep[[i]] <-
  list(Time= read.table(paste(DataPath,"T",sep=""))$V1, Y=matrix(tmpData[[i]],nrow=1))



###############################################
# Model 1 - ISR as random walk.
###############################################

k1 = 0.053; k2 = 0.051; ke = 0.062;
Model1 <- list(
               Matrices=function(phi) {
                 list(
                      matA=matrix(c(-(k1+ke),k2,1,k1,-k2,0,0,0,0.00001),ncol=3,byrow=T),
                      matB=NA,
                      matC=matrix(c(1,0,0),nrow=1),
                      matD=NA )
               },
               X0 = function(Time,phi,U) {
                 
                 C0 <- phi[["C0"]]
                 matrix(c(C0, C0*k1/k2, C0*ke),ncol=1)
               },
               SIG = function(phi) {
                 diag( c(0,0,phi[["SIG33"]]))
               },
               S = function(phi) {
                 matrix(phi[["S"]])
               },
               h = function(eta,theta,covar) {
                 phi <- theta
                 phi[["C0"]] <- theta[["C0"]]*exp(eta[1])
                 phi
               },
               ModelPar = function(THETA) {
                 list(theta=list(C0=THETA[1],S=THETA[2],SIG33=THETA[3]),
                      OMEGA=matrix(THETA[4]))
               }
               )

# From Matlab version of PSM
# THETA = (912.50196725972,   8574.60240381849,   6.06399638066, 0.16988998197)
par1 <- list(LB   = c(  200,  50^2,   0,  .0 ),
             Init = c( 1000, 100^2,  10,  .25),
             UB   = c( 1800, 150^2,  15,  .50))

Rprof(filename = "Rprof1.out")
fit1 <- PSM.estimate(Model=Model1,Data=Cpep,Par=par1,CI=T,trace=1,control=list(optimizer='ucminf',trace=TRUE))
Rprof(NULL)

Rprof(filename = "Rprof2.out")
fit2 <- PSM.estimate(Model=Model1,Data=Cpep,Par=par1,CI=T,trace=1,fast=F)
Rprof(NULL)

summaryRprof(filename = "Rprof1.out")
summaryRprof(filename = "Rprof2.out")

fit1[1:5] 
                                        #$NegLogL [1] 3052.865,
                                        #Runtime:  23:23.02 (linux04, CI=F)

                                        #med L-BFGS-B p? indvendig opt:
                                        #$NegLogL [1] 3052.866
                                        #Runtime:  18:15.79 (linux06, CI=T)

smooth1 <- PSM.smooth(Model=Model1,Data=Cpep,THETA=fit1$THETA,sub=10)



PSM.plot(Cpep,smooth1,indiv=3:5,type=c('Xs','Ys.Y','res','acf','eta'))




par(mfcol=c(3,2))
for(j in 2:3)
  for(i in 1:3) {
    plot(smooth1[[j]]$Time,smooth1[[j]]$Xs[i,],type="l",
         ylab=paste('state',i), xlab=paste('individual',j))
    if(i==1) points(Cpep[[i]]$Time,Cpep[[j]]$Y)
    rug(Cpep[[i]]$Time)
  }




###############################################
# Model 1b - w/o random effect on C0
###############################################

Model1b <- Model1
Model1b$ModelPar = function(THETA){
              return(list(theta=list(C0=THETA[1],S=THETA[2],SIG33=THETA[3]),
                          OMEGA=NULL))}
Model1b$h <- function(eta,theta,covar) {theta}
par1b <- list(LB  = c(  200,  50^2,   0 ),
             Init = c( 1000, 100^2,  10 ),
             UB   = c( 1800, 150^2,  20 ))

fit1b <- PSM.estimate(Model=Model1b,Data=Cpep,Par=par1b,CI=TRUE,trace=1)
fit1b[1:3]
                                        #final  value 3071.008968 

smooth1b <- PSM.smooth(Model=Model1b,Data=Cpep,THETA=fit1b$THETA,sub=10,trace=0)

PSM.plot(Cpep,smooth1b,indiv=3:8,type=c('Xs','Yp.Y','res','acf','eta'))

par(mfcol=c(3,2))
for(j in 2:3)        # Note worse fit for initial obs.
for(i in 1:3) { 
  plot(smooth1b[[j]]$Time,smooth1b[[j]]$Xs[i,],type="l",
       ylab=paste('state',i), xlab=paste('individual',j))
  if(i==1) points(Cpep[[i]]$Time,Cpep[[j]]$Y)
  rug(Cpep[[1]]$Time)
} 




##################################################################
# Model 2 - ISR as intervention model with
#           random effects for peaks, baseline and C0.
##################################################################

Cpep2 <- Cpep
for(i in 1:12) {
  Cpep2[[i]]$U <- matrix(c(1,1,rep(0,33),
                               rep(0,11), 1, rep(0,23),
                               rep(0,18), 1, 1, rep(0,15),
                               rep(1,35)),byrow=T,nrow=4)
}
k1 = 0.053; k2 = 0.051; ke = 0.062;
Model2 <- list(
               Matrices = function(phi) {
                 a1  <- phi[["a1"]]
                 a2  <- phi[["a2"]]
                 B  <- phi[["B"]]
                 K  <- phi[["K"]]
                 a2K1 <- phi[["a2K1"]]
                 a2K2 <- phi[["a2K2"]]
                 a2K3 <- phi[["a2K3"]]
                 matA <- matrix( c(-(k1+ke) ,  k2 ,   1 ,   0,
                                   k1 , -k2 ,   0 ,   0,
                                   0 ,   0 , -a1 ,  a1,
                                   0 ,   0 ,   0 , -a2),nrow=4,byrow=T)
                 matB <- matrix( c(0           , 0          , 0               , 0,
                                   0           , 0          , 0               , 0,
                                   0           , 0          , 0               , B,
                                   a2K1        , a2K2       , a2K3            , 0),
                                byrow=T,nrow=4)
                 matC <- matrix(c(1,0,0,0),nrow=1)
                 matD <- matrix(c(0,0,0,0),nrow=1)
                 list(matA=matA,matB=matB,matC=matC,matD=matD)
               },
               X0 = function(Time,phi,U) {
                 C0 <- phi[["C0"]]
                 matrix(c(C0,C0*k1/k2,C0*ke,0),ncol=1)
               },
               SIG = function(phi) {
                 diag( c(0,0,phi[["SIG33"]],0)) 
               },
               S = function(phi) {
                 matrix(phi[["S"]])
               },
               h = function(eta,theta,covar) {
                 phi <- theta
                 phi[["B"]] <- theta[["B"]]*exp(eta[1])
                 phi[["a2K1"]] <- theta[["a2"]]*theta[["K"]]*exp(eta[2])
                 phi[["a2K2"]] <- theta[["a2"]]*theta[["K"]]*exp(eta[3])
                 phi[["a2K3"]] <- theta[["a2"]]*theta[["K"]]*exp(eta[4])
                 phi[["C0"]] <- theta[["C0"]]*exp(eta[5])
                 phi
               },
               ModelPar = function(THETA){
                 list(theta=
                      list(C0=900, S=8500, a1=THETA[1], a2=THETA[2],
                           SIG33=THETA[3], K=THETA[4], B=THETA[5]),
                      OMEGA=diag(c(.2,2,2,2,.2)))
               }
               )


par2 <- list(LB   = c(   0.01, 0.005,  2     ,  200     ,    1     ),
             Init = c( 0.0263, 0.0097, 4.7853,  422.2910,    1.7505), 
             UB   = c(   0.05, 0.020,  8     ,  600     ,    2.5   ))
#ovenst?ende startv?rdier genskaber AAA matrice, der giver problemer for Matrix:expm.
par2$Init <- c( 0.03, 0.01, 5,  400,  2)
par2$LB <- par2$Init/7
par2$UB <- par2$Init*7

fit2 <- PSM.estimate(Model=Model2,Data=Cpep2,Par=par2,CI=T,trace=1)
fit2[1:3]
                                        # Runtime:  566:1.2 > $NegLogL [1] 2950.34
                                        # Runtime:  235:4.6 > $NegLogL [1] 2950.357 (2. s?t par) $THETA [1] 2.605852e-02 9.804498e-03 4.815092e+00 3.973301e+02 1.655686e+00
                                        # med L-BFGS-B p? intern opt:
                                        # Runtime:  261:44.0 > $NegLogL 2950.348 (factr 1e10) (diag cov neg..)

smooth2 <- PSM.smooth(Model=Model2,Data=Cpep2,THETA=fit2$THETA,sub=10)

PSM.plot(Cpep,smooth2,indiv=3:5,type=c('Xs','Ys.Y','res','acf','eta'))

par(mfcol=c(4,2))
for(j in 2:3)
for(i in 1:4) {
  plot(smooth2[[j]]$Time,smooth2[[j]]$Xs[i,],type="l",
       ylab=paste('state',i), xlab=paste('individual',j))
  if(i==1) points(Cpep2[[j]]$Time,Cpep2[[j]]$Y)
  rug(Cpep[[1]]$Time)
}



##################################################################
# Model 2b - ISR as intervention model with
#            random effects for peaks, baseline and C0.
#            *NAs inserted in data.*
##################################################################

Cpep2b <- Cpep2
for (j in 1:12) {
  Cpep2b[[j]]$Y <- matrix(NA,ncol=length(Cpep2b[[j]]$Time)*2,nrow=2)
  Cpep2b[[j]]$U <- matrix(NA,ncol=length(Cpep2b[[j]]$Time)*2,nrow=4)
  for(i in 1:length(Cpep2b[[j]]$Time)) {
    Cpep2b[[j]]$Y[(i%%2)+1,i*2-1] <- Cpep2[[j]]$Y[1,i]
    Cpep2b[[j]]$U[,i*2-1] <- Cpep2b[[j]]$U[,i*2] <- Cpep2[[j]]$U[,i]
  }
  Cpep2b[[j]]$Y <- Cpep2b[[j]]$Y[,-2] #remove to preserve tau=t2-t1 for P0
  Cpep2b[[j]]$U <- Cpep2b[[j]]$U[,-2]
  Cpep2b[[j]]$Time <- c(0,15,22,30,37,45,52,60,67,75,82,90,105,120,135,150,165,180,195,
                        210,225,240,255,270,285,300,315,330,345,360,390,420,450,480,540,600,607,
                        615,622,630,637,645,652,660,667,675,682,690,705,720,735,750,765,780,795,
                        810,825,840,850,959,1050,1140,1230,1320,1365,1410,1425,1440,1455)
}
# Show the new data with NAs
rbind(Cpep2b[[1]]$Y[,1:8],
      Cpep2b[[1]]$U[,1:8],
      Cpep2b[[1]]$Time[1:8])
rbind(Cpep2[[1]]$Y[,1:4],
      Cpep2[[1]]$U[,1:4],
      Cpep2[[1]]$Time[1:4])

# Modify model to fit 2-dim Y
Model2b <- Model2
Model2b$Matrices <- function(phi) {
  a1  <- phi[["a1"]]
  a2  <- phi[["a2"]]
  B  <- phi[["B"]]
  K  <- phi[["K"]]
  a2K1 <- phi[["a2K1"]]
  a2K2 <- phi[["a2K2"]]
  a2K3 <- phi[["a2K3"]]
  matA <- matrix( c(-(k1+ke) ,  k2 ,   1 ,   0,
                    k1 , -k2 ,   0 ,   0,
                    0 ,   0 , -a1 ,  a1,

                    0 ,   0 ,   0 , -a2),nrow=4,byrow=T)
  matB <- matrix( c(0           , 0          , 0               , 0,
                    0           , 0          , 0               , 0,
                    0           , 0          , 0               , B,
                    a2K1        , a2K2       , a2K3            , 0),
                 byrow=T,nrow=4)
  matC <- matrix(rep(c(1,0,0,0),each=2),nrow=2)
  matD <- matrix(rep(0,8),nrow=2)
  list(matA=matA,matB=matB,matC=matC,matD=matD)
}
Model2b$S = function(phi) {
                 diag(phi[["S"]],2)
               }
ModelCheck( Model=Model2b , Data=Cpep2b[[1]], Par=par2)


# Test Kalman-Filter
theta2b <- Model2b$ModelPar(par2$Init)$theta
phi2b <- Model2b$h(eta=rep(0,5),theta=theta2b,covar=NULL)
Ob1 <- LinKalmanFilter( phi=phi2b , Model=Model2b , Data=Cpep2b[[1]] , output = TRUE, echo=FALSE,fast=TRUE)
Ob1$negLogLike
Ob2 <- LinKalmanFilter( phi=phi2b , Model=Model2  , Data=Cpep2[[1]] , output = TRUE, echo=FALSE,fast=TRUE)
Ob2$negLogLike

Ob2$R[,,1:3]
Ob1$R[,,1:4]
Ob2$KfGain[,,1:3]
Ob1$KfGain[,,1:4]
Ob2$Yp[,1:4]
Ob1$Yp[,1:4]
Cpep2b[[1]]$Y[,1:4]

# TEST APL-evaluering 
# Uden NAs
APL.KF(par2$Init,Model2,Cpep2,GUIFlag=2)
# Med NAs
APL.KF(par2$Init,Model2b,Cpep2b,GUIFlag=2)



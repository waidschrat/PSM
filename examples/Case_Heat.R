# Test Case
# Originating from CTSM
# Linear Time Invarient - Heat Case

rm(list=ls())

PC <- TRUE

if(PC){
  if( is.loaded('PSM')) detach(package:PSM) 
  library(PSM)
} else {
  detach(package:PSM)

  library(PSM,lib.loc="~/PSM/Rpackages/gridterm")
  source("../R/PSM.estimate.R",echo=F)
  source("../R/invlogit.R",echo=F)
  source("../R/LinKalmanFilter.R",echo=F)
}


# Load the Data and Variables
tmpData <- read.table("Heat_Data.csv",sep=";", col.names=c("TIME","Te","Ti","Q"))

Time=tmpData$TIME
Y=t(matrix(tmpData[,c("Q")]))
U=t(as.matrix(tmpData[,c("Te","Ti")]))
Pop.Data <- list( list(Time=tmpData$TIME,
                Y=t(matrix(tmpData[,c("Q")])),
                U=t(as.matrix(tmpData[,c("Te","Ti")]))) )

HeatModel <- list(
                  Matrices = function(phi=NA) {
                    G1  <- phi[["G1"]] ; G2  <- phi[["G2"]]
                    H1  <- phi[["H1"]] ; H2  <- phi[["H2"]] ; H3  <- phi[["H3"]]
                    tmp <- list(
                                matA = matrix( c(-1*(1/H1+1/H2)/G1,1/(G1*H2),1/(G2*H2) , -1*(1/H2+1/H3)/G2 ) , ncol=2, byrow=T),
                                matB = diag( c(1/(G1*H1) , 1/(G2*H3) ) ),
                                matC = matrix( c(0,-1/H3) ,nrow=1),
                                matD = matrix( c(0,1/H3)  ,nrow=1))
                    return(tmp)
                  },
                  X0 = function(Time=NA,phi=NA,U=NA) {
                    tmp    <- phi[["X01"]]
                    tmp[2] <- phi[["X02"]]
                    return(matrix(tmp,ncol=1) )} ,
                  SIG = function(phi=NA) {
                    return( diag( c(phi[["SIG11"]],phi[["SIG22"]])))} ,
                  S = function(phi=NA) {
                    return( matrix(phi[["S"]])) } ,
                  h = function(eta,theta,covar=NULL) {
                    phi <- theta
                    return(phi) } ,
                  ModelPar = function(THETA){
                    return(list(theta=list( G1=THETA[1],G2=THETA[2],
                                  H1=THETA[3],H2=THETA[4],H3=THETA[5],
                                  SIG11=THETA[6], SIG22=THETA[7], S=THETA[8],
                                  X01=THETA[9], X02=THETA[10]),
                                OMEGA=NULL))},
                  Dose=NULL
                  )

names(HeatModel)

# -------------------------------------------------------------
# Test of Fortran Code
# -------------------------------------------------------------

Testphi <- c( 13 ,25 , 100 , 1 , 2 , 49 , .5 , .2 , .2 , 0.01)
names(Testphi) <- c("X01","X02","G1","H1","H2","G2","H3","SIG11","SIG22","S")

Ob1 <- LinKalmanFilter( phi=Testphi , Model=HeatModel , Data=Pop.Data[[1]] , echo=FALSE, outputInternals=TRUE,fast=TRUE)
Ob2 <- LinKalmanFilter( phi=Testphi , Model=HeatModel , Data=Pop.Data[[1]] , echo=FALSE, outputInternals=TRUE,fast=FALSE)
names(Ob1)
Ob1$negLogLike
Ob2$negLogLike

IDX = 1:7

Ob1$Yp[,1]
Ob2$Yp[,1]

Ob1$Pp[,,1]
Ob2$Pp[,,1]

(TestmatC <- matrix( c(0,-2) ,nrow=1))
(TestS <- Testphi["S"])
TestmatC %*% Ob1$Pp[,,1] %*% t(TestmatC) + TestS

Ob1$R[,,1]
Ob2$R[,,1]



# -------------------------------------------------------------
# Test of Fortran Code
# -------------------------------------------------------------

# Parameter estimation
# Initial guess from CTSM
# THETA OBJ             G1,   G2,  H1,  H2,  H3, SIG11,SIG22,    S, X01,  X02
# CTSM starting guess fails in this implementation
par1 <- list(LB   = c(  10,   10,1e-1,1e-1,1e-2,  1e-8, 1e-8, 1e-4,   10,   20),
             Init = c( 100,   50,   1,   2,  .5,   .01,  .01,  .01,   15,   25),
             UB   = c( 200,  100,   2,   5,   1,     1,    1,    1,   20,   30)
             )
APL.KF(par1$Init,HeatModel,Pop.Data)
#par1$Init <-        c( 100 ,  50,   1,   2,  .5,  .001, .001, .001,   13,   25)
#par1$UB <- par1$LB <- NULL


# Check the Model
ModelCheck( Model=HeatModel , Data=Pop.Data[[1]], Par=par1)


# -------------------------------------------------------------
# Test Linear Kalman Filter with CTSM estimated parameters
# -------------------------------------------------------------

# CTSM returns -LL= -623 
CTSMphi <- c( 1.3134E+01,2.5330E+01,1.0394E+02,9.6509E-01,2.0215E+00,4.9320E+01,5.0929E-01,7.3779E-08,2.6951E-09,1.0330E-02)
names(CTSMphi) <- c("X01","X02","G1","H1","H2","G2","H3","SIG11","SIG22","S")
CTSMTHETA=c(CTSMphi[["G1"]],   CTSMphi[["G2"]],  CTSMphi[["H1"]],  CTSMphi[["H2"]],  CTSMphi[["H3"]], CTSMphi[["SIG11"]],CTSMphi[["SIG22"]],    CTSMphi[["S"]], CTSMphi[["X01"]],  CTSMphi[["X02"]])
Ob1 <- LinKalmanFilter( phi=CTSMphi , Model=HeatModel , Data=Pop.Data[[1]] , echo=FALSE, outputInternals=TRUE,fast=FALSE)
Ob1$negLogL #[1,] -623.3564 in R
APL.KF(CTSMTHETA,HeatModel,Pop.Data)


# Validation plot versus Data
D <- Pop.Data[[1]]
plot(D$Time , D$Y,col="red")
points( D$Time , Ob1$Yp, pch="+",col="blue")

# -------------------------------------------------------------
# Minimizers
# -------------------------------------------------------------
                                        # Test Run of initial parameters

# Perform minimization with 2 different optimizers
#Rprof()
Min1 <- PSM.estimate(Model=HeatModel,Data=Pop.Data,Par=par1,CI=TRUE,trace=2,optimizer="nlm")
#Rprof(NULL)
#summaryRprof()

Min2 <- PSM.estimate(Model=HeatModel,Data=Pop.Data,Par=par1,CI=TRUE,trace=2,optimizer="optim")

cat( "nlm: "   , Min1$sec, "\t ", Min1$opt$minimum, "\n")
cat( "optim: " , Min2$sec, "\t ", Min2$opt$value, "\n")


# -------------------------------------------------------------
# Smoother 
# -------------------------------------------------------------

SmoothObj <- PSM.smooth(Model=HeatModel, Data=Pop.Data, THETA=CTSMTHETA, subsample=0,trace=1)


D <- SmoothObj[[1]]
names(D)
D$negLogL

Idx <- 200:500
plot( D$Time[Idx], D$Xs[1,Idx] , type="n" )
for(i in 1:2) polygon( c(D$Time[Idx],rev(D$Time[Idx])) , c(D$Xs[i,Idx],rev(D$Xs[i,Idx]))+sqrt( abs(c(D$Ps[i,i,Idx], - rev(D$Ps[i,i,Idx])))),col=4)
for(i in 1:2) lines( D$Time[Idx], D$Xs[i,Idx], type="l",lwd=2)



# Load the predicted Data from CTSM
CTSM.Val.Data <- read.table("CTSM_Heat_Xp.csv",sep=";", col.names=c("Time","Xp1","Xp2","SDX1","SDX2","Y1","SDY1"))

names(D)


# Measurements Validation
plot(Pop.Data[[1]]$Time, Pop.Data[[1]]$Y  )
lines( CTSM.Val.Data$Time, CTSM.Val.Data$Y1, col="red")
Idx <- !is.na(D$Yp)
lines( D$Time[Idx] , D$Yp[Idx] , col="blue")
plot(CTSM.Val.Data$Y1[Idx]/D$Yp[Idx])

# State validation
plot(D$Time[-(1:10)], D$Xp[1,-(1:10)],type="l",col="red")
lines(D$Time[-(1:10)], D$Xp[2,-(1:10)],type="l",col="red")

plot(CTSM.Val.Data$Time , CTSM.Val.Data$Xp2,type="l",col="blue",ylim=range(CTSM.Val.Data[,c("Xp1","Xp2")]))
lines(CTSM.Val.Data$Time , CTSM.Val.Data$Xp1,type="l",col="blue")




# -------------------------------------------------------------
# Include NAs in data 
# -------------------------------------------------------------


Pop.DataNA <- Pop.Data
Pop.DataNA[[1]]$Y <- matrix(NA,ncol=length(Pop.DataNA[[1]]$Time)*2,nrow=2)
Pop.DataNA[[1]]$U <- matrix(NA,ncol=length(Pop.DataNA[[1]]$Time)*2,nrow=2)
for(i in 1:length(Pop.DataNA[[1]]$Time)) {
  Pop.DataNA[[1]]$Y[(i%%2)+1,i*2-1] <- Pop.Data[[1]]$Y[1,i]
  Pop.DataNA[[1]]$U[,i*2-1] <- Pop.DataNA[[1]]$U[,i*2] <- Pop.Data[[1]]$U[,i]
}
Pop.DataNA[[1]]$Time <- seq(0,718.5,.5)

# Show new Y, U and Time
rbind(Pop.DataNA[[1]]$Y[,1:8],
      Pop.DataNA[[1]]$U[,1:8],
      Pop.DataNA[[1]]$Time[1:8])

# Update model til 2-dim Y
HeatModelNA <- HeatModel
HeatModelNA$Matrices <- function(phi=NA) {
  G1  <- phi[["G1"]] ; G2  <- phi[["G2"]]
  H1  <- phi[["H1"]] ; H2  <- phi[["H2"]] ; H3  <- phi[["H3"]]
  list(
       matA = matrix( c(-1*(1/H1+1/H2)/G1,1/(G1*H2),1/(G2*H2) , -1*(1/H2+1/H3)/G2 ) , ncol=2, byrow=T),
       matB = diag( c(1/(G1*H1) , 1/(G2*H3) ) ),
       matC = matrix( rep(c(0,-1/H3),each=2) ,nrow=2),
       matD = matrix( rep(c(0,1/H3),each=2)  ,nrow=2))
}
HeatModelNA$S <- function(phi=NA) {
   diag(rep(phi[["S"]],2)) }

ModelCheck( Model=HeatModelNA , Data=Pop.DataNA[[1]], Par=par1)
APL.KF(CTSMTHETA,HeatModelNA,Pop.DataNA) #[1] -623.3564

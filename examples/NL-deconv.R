library(PSM)

k1 = 0.053; k2 = 0.051; ke = 0.062;
MyModel <- vector(mode="list")
MyModel$Functions <- 
  list(
       f = function(x,u,time,phi) {
         matrix(c(-(k1+ke),k2,1,k1,-k2,0,0,0,0), nrow=3, ncol=3,byrow=TRUE)%*%x
       },
       df = function(x,u,time,phi) {
         matrix(c(-(k1+ke),k2,1,k1,-k2,0,0,0,0), nrow=3, ncol=3,byrow=TRUE)
       },
       g = function(x,u,time,phi) {
         x[1,1,drop=FALSE]
       },
       dg = function(x,u,time,phi) {
         matrix(c(1,0,0), nrow=1, ncol=3)
       }
       )
MyModel$h = function(eta,theta,covar) {
  phi <- theta
  phi
}
MyModel$S = function(u,time,phi) {
  matrix(phi$S, nrow=1, ncol=1)
}
MyModel$SIG = function(u,time,phi) {
  matrix(c(0,0,0, 0,0,0, 0,0,phi$SIG33  ), nrow=3, ncol=3)
}
MyModel$X0 = function(Time,phi,U) {
  matrix(c(0,0,50), nrow=3, ncol=1)
}
MyModel$ModelPar = function(THETA) {
  list(theta=list(S=THETA['S'],#SIG33=0.1
         SIG33=THETA['SIG33']
         ))
}

#TEST numerical df + dg
MyModel$Functions$df = function(x,u,time,phi) {
  jacobian(MyModel$Functions$f,x=x,u=u,time=time,phi=phi)
}
MyModel$Functions$dg = function(x,u,time,phi) {
  jacobian(MyModel$Functions$g,x=x,u=u,time=time,phi=phi)
}


TimeData <- list(list( Time = seq(0,100,10,
                       Dose = list(Time=c(10,20,40,50,50), State=c(3,3,3,2,1),
                         Amount=c(100,100,-100,300,4000)))),
                 list( Time = seq(0,100,10,
                       Dose = list(Time=c(10,20,40,50,50), State=c(3,3,3,2,1),
                         Amount=c(100,100,-100,300,4000))))
                 )

MyTHETA <- c(S = 5000,SIG33=5)
Data <- PSM.simulate(MyModel,TimeData,MyTHETA,deltaTime=.1)

PSM.plot(Data,type=c('Y','longX'))

par(mfrow=c(3,1))
for(i in 1:3) {
  plot(Data[[1]]$longTime,Data[[1]]$longX[i,],type="l")
  if(i==1) points(Data[[1]]$Time,Data[[1]]$Y)
}

par <- list(LB = MyTHETA*.1,
            Init = MyTHETA,
            UB = MyTHETA*2)

system.time(
            fit <- PSM.estimate(MyModel,Data,par,trace=1,CI=TRUE)
            )
fit[1:5]

smooth <- PSM.smooth(MyModel,Data,fit$THETA,subs=40)
sm <- smooth[[1]]
PSM.plot(Data,smooth,type=c('Ys.Y','Xs'))

par(mfrow=c(3,1))
for(i in 1:3) {
  plot(Data[[1]]$longTime,Data[[1]]$longX[i,],type="l")
  if(i==1) points(Data[[1]]$Time,Data[[1]]$Y)
  lines(sm$Time,sm$Xs[i,],col=2)
}

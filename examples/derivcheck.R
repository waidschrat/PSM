# Test problems for the checkgrad implementation in R

rm(list=ls())

library(numDeriv)

# Extracts from the BioReactor 
ModelPar = function(THETA) {
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

f = function(x,u,time,phi) {
  X <- x[1]; S <- x[2]; V <- x[3]; F <- u[1]
  mu <- phi$mumax*S/(phi$k2*S^2+S+phi$k1)
  matrix(c(
           mu*X-F*X/V,
           -mu*X/phi$Y+F*(phi$sf-S)/V,
           F
           ),ncol=1)
}
df = function(x,u,time,phi) {
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
X0 = function(Time,phi,U) {
  matrix(c(phi$x0,phi$s0,phi$v0),ncol=1)
}
           
ctsmTHETA <- c(mumax=1.0022E+00,k1=3.1629E-02,
             s11=7.5248E-03,s22=1.0636E-03,s33=1.1388E-02)
u0 = 0.2129483
phi <- ModelPar(ctsmTHETA)$theta
x <- X0(phi=phi)


# Test of the Jacobian
f(x=x,u=u0,phi=phi)

df1 <- df(x=x,u=u0,phi=phi)
jc1 <- jacobian(func=f,x=x,phi=phi,u=u0)

all.equal2(df1,1.0001*jc1)
all.equal2(df1,jc1)



all.equal2 <- function(target,current,...) {
  out <- all.equal(target,current,...)
  b <- isTRUE(out)
  attributes(b) <- list(msg=as.character(out))
  b
}

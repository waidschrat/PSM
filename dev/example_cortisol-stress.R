library(MASS)
files <- list.files("C:/Users/Robert Miller/Dropbox/R Packages/PSM/PSM/R")
for(i in 1:length(files)) source(file = paste0("C:/Users/Robert Miller/Dropbox/R Packages/PSM/PSM/R/",files[i]), local = T)
source("PSM_aux.R", local = T)


setwd("C:/Users/Robert Miller/Google Drive/Studien/2023/EUSTRESS/Simdata")
PKdata <- nlme::groupedData(DV ~ TAD | ID, read.csv("CortStressResponse.csv")[,-1])


#format data for PSM modeling functions
N <- length(unique(PKdata$ID))
PKdat <- vector(mode="list",length=N)
for(i in 1:N){
  PKdat[[i]]$Time <- with(PKdata[as.numeric(PKdata$ID) == i,], TAD)
  PKdat[[i]]$Y <- with(PKdata[as.numeric(PKdata$ID) == i,], matrix(DV, nrow=1))
  if(PKdata$FLAG[as.numeric(PKdata$ID) == i][1] == 2){
    PKdat[[i]]$Time <- c(PKdat[[i]]$Time,0,1)
    PKdat[[i]]$Y <- matrix(c(PKdat[[i]]$Y,NA,NA)[order(PKdat[[i]]$Time)], nrow=1)
    PKdat[[i]]$Time <- sort(PKdat[[i]]$Time)
  }
  if(PKdata$FLAG[as.numeric(PKdata$ID) == i][1] == 1){
    PKdat[[i]]$Time <- c(PKdat[[i]]$Time,1)
    PKdat[[i]]$Y <- matrix(c(PKdat[[i]]$Y,NA)[order(PKdat[[i]]$Time)], nrow=1)
    PKdat[[i]]$Time <- sort(PKdat[[i]]$Time)
  }
  PKdat[[i]]$U <- matrix(c(ifelse(PKdat[[i]]$Time != 0, 0, 1),
                           rep(1, length(PKdat[[i]]$Time))
  ), nrow=2, byrow=T)
  
  PKdat[[i]]$covar <- c(STUDY = PKdata$FLAG[as.numeric(PKdata$ID) == i][1] - 1,
                        FEMALE = PKdata$SEX[as.numeric(PKdata$ID) == i][1])
}

#compile parameter estimates and fit of population ODE/SDE models
modt <- readRDS("PSM_transits_ODE.RDS")
mods <- readRDS("PSM_population_ODE.RDS")
i <- length(mods)
pars <- c(round(modt[[4]]$THETA,1), "OMEGA_stress"=1.5, "OMEGA_ke"=.05, "OMEGA_kt"=.15, "OMEGA_init"=.15)
pars[c("init","sigma")] <- c(.001,.001)
parA <- list(LB=pars*.2, Init=pars, UB=pars*2.5) #bounds + inits
npars <- Vectorize(function(x) sum(parA$Init != round(mods[[x]]$THETA,3)))(1:i)

res <- Vectorize(function(x) round(mods[[x]]$THETA,3))(1:i) #parameter estimates for each fitted model
res <- rbind(res, "Mtt"= Vectorize(function(x) round(4/mods[[x]]$THETA["kt"],3))(1:i)) #add mean transit time for each model
res <- rbind(res, "LL"= -Vectorize(function(x) mods[[x]]$NegLogL)(1:i))
res <- rbind(res, "AIC"= 2*npars +2*Vectorize(function(x) mods[[x]]$NegLogL)(1:i))
res <- rbind(res, "R2"= Vectorize(function(x) 1 - mods[[x]]$THETA["S"]/var(PKdata$DV))(1:i))
colnames(res) <- letters[1:i]
print(res)


#visual predictive checks (Figure 6; model e vs model f)
N <- 2000
simdat <- list()
for(i in 1:(N/2))  simdat[[i]] <- list("Time" = -10:75, "U" = matrix(c(ifelse(-10:75 == 1, 1, 0), ifelse(-10:75 == 100, 0, 1)), nrow=2, ncol=length(-10:75), byrow = T),
                                       "covar" = c("STUDY" = 1, "FEMALE" = 0))
for(i in (N/2+1):N)  simdat[[i]] <- list("Time" = -10:75, "U" = matrix(c(ifelse(-10:75 == 1, 1, 0), ifelse(-10:75 == 100, 0, 1)), nrow=2, ncol=length(-10:75), byrow = T),
                                         "covar" = c("STUDY" = 1, "FEMALE" = 1))
PKsimODE <- PSM.simulate(Model = ODE_lme(3,4), THETA = res[1:14,"e"], Data = simdat, deltaTime = .1)
PKsimSDE <- PSM.simulate(Model = SDE_lme(3,4), THETA = res[1:14,"f"], Data = simdat, deltaTime = .1)

predquant <- Vectorize(function(time, female, perc, simdata){
  if(female == 0) return(quantile(unlist(Vectorize(function(x) simdata[[x]]$Y[time+11])(1:(N/2))), perc))
  if(female == 1) return(quantile(unlist(Vectorize(function(x) simdata[[x]]$Y[time+11])((N/2+1):N)), perc))
}, vectorize.args = c("time"))

par(mfrow=c(1,2), mar=c(3.5,3.5,1,1.5), mgp=c(2.2,0.8,0))

plot(c(-10,75), c(0,70), type="n", xlab="time relative to TSST onset (min)", ylab="salivary cortisol (nM)", axes=F) #females
grid()
axis(1, at = seq(-10,70,10))
axis(2)
poly_y <- c(mA(predquant(time = -10:75, female = 1, perc = .025, PKsimODE)), mA(predquant(time = 75:-10, female = 1, perc = .975, PKsimODE)))
polygon(x = c(-10:75,75:-10), y = poly_y, col = rgb(.5,.5,.5,.5), border = NA) #ODE
lines(-10:75, mA(predquant(time = -10:75, female = 1, perc = .5, PKsimODE)), lwd=2, lty=2, col="darkgrey")
poly_y <- c(mA(predquant(time = -10:75, female = 1, perc = .025, PKsimSDE)), mA(predquant(time = 75:-10, female = 1, perc = .975, PKsimSDE)))
polygon(x = c(-10:75,75:-10), y = poly_y, col = rgb(135/255,206/255,235/255,.5), border = NA) #SDE
#with(PKdata[PKdata$SEX == 1 & PKdata$FLAG == 2,], points(TAD, DV, col=rgb(0,0,0,.1), pch=19))
with(PKdata[PKdata$SEX == 1 & PKdata$FLAG == 2,], vioplotm(DV[TAD == -5], DV[TAD == 11], DV[TAD == 20], DV[TAD == 30], DV[TAD == 40], DV[TAD == 55], DV[TAD == 70], at=c(-5,11,20,30,40,55,70), range = 0, wex=4, add=T, axes=F))
text(0, 60, "A", cex=2.5)
legend("topright", bty = "n", col=c("darkgrey","skyblue"), lwd=10, lty=1, legend=c("ODE","SDE"))

plot(c(-10,75), c(0,70), type="n", xlab="time relative to TSST onset (min)", ylab="salivary cortisol (nM)", axes=F) #females
grid()
axis(1, at = seq(-10,70,10))
axis(2)
poly_y <-c(mA(predquant(time = -10:75, female = 0, perc = .05, PKsimODE)), mA(predquant(time = 75:-10, female = 0, perc = .95, PKsimODE)))
polygon(x = c(-10:75,75:-10), y = poly_y, col = rgb(.5,.5,.5,.5), border = NA) #ODE
lines(-10:75, mA(predquant(time = -10:75, female = 0, perc = .5, PKsimODE)), lwd=2, lty=2, col="darkgrey")
poly_y <- c(mA(predquant(time = -10:75, female = 0, perc = .05, PKsimSDE)), mA(predquant(time = 75:-10, female = 0, perc = .95, PKsimSDE)))
polygon(x = c(-10:75,75:-10), y = poly_y, col = rgb(135/255,206/255,235/255,.5), border = NA) #SDE
#with(PKdata[PKdata$SEX == 1 & PKdata$FLAG == 2,], points(TAD, DV, col=rgb(0,0,0,.1), pch=19))
with(PKdata[PKdata$SEX == 0 & PKdata$FLAG == 2,], vioplotm(DV[TAD == -5], DV[TAD == 11], DV[TAD == 20], DV[TAD == 30], DV[TAD == 40], DV[TAD == 55], DV[TAD == 70], at=c(-5,11,20,30,40,55,70), range = 0, wex=4, add=T, axes=F))
text(0, 60, "B", cex=2.5)



#simulate artificial cortisol data using the population SDE model (model f)
N <- 10
samps <- -20:80
simdat <- list()
for(i in 1:(N/2))  simdat[[i]] <- list("Time" = samps, "U" = matrix(c(ifelse(samps == 1, 1, 0), ifelse(samps == 100, 0, 1)), nrow=2, ncol=length(samps), byrow = T),
                                       "covar" = c("STUDY" = 0, "FEMALE" = 0))
for(i in (N/2+1):N)  simdat[[i]] <- list("Time" = samps, "U" = matrix(c(ifelse(samps == 1, 1, 0), ifelse(samps == 100, 0, 1)), nrow=2, ncol=length(samps), byrow = T),
                                         "covar" = c("STUDY" = 0, "FEMALE" = 1))

set.seed(1234)
PKsim <- PSM.simulate(Model = SDE_lme(3,4), THETA = res[1:14,"f"], Data = simdat, deltaTime = 1)
pars <- matrix(Vectorize(function(x) c(res["stress","f"], res["ke","f"], res["kt","f"], res["ks","f"]/res["ke","f"]) * exp(PKsim[[x]]$eta[1:4]) * exp(c(-res["sex","f"]*simdat[[x]]$covar["FEMALE"],0,0,0)))(1:N), ncol=4, byrow=T)


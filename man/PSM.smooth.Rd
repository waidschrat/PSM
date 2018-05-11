\name{PSM.smooth}
\alias{PSM.smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smoothing of model states based on estimated population parameters.}
\description{
  Gives estimates of model states and random effects \eqn{\eta}. The
  function is intended to be used based on population parameters found
  using \code{\link{PSM.estimate}} or to check initial values before
  parameter estimation.
}

\usage{
PSM.smooth(Model, Data, THETA, subsample = 0, trace = 0, etaList = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Model}{
    Model list.*
  }
  \item{Data}{
    Data list.*
  }
  \item{THETA}{
    Vector of population parameters used for the state estimation.
  }
  \item{subsample}{
    Number of points to estimate states in between
    measurements. The extra points are linearly spaced.
  }
  \item{trace}{
    Non-negative integer. If positive, tracing
    information on the progress of the optimization is produced. Higher
    values produces more tracing information.
  }
  \item{etaList}{
    Matrix where each column contains an etimate of
    \eqn{\eta_i}. \code{etaList} has the same format as the output of 
    \code{\link{PSM.estimate}}. If ommitted, the function will evalutate
    the population likehood function to find estimates of \eqn{\eta_i}{eta}
    for all individuals.
  }
  * See description in \link{PSM.estimate}.
}
\details{
  The function produces three types of estimates.

  \describe{
    \item{Predicted}{Only past measurements are used for the state
      estimate at time t. }
    
    \item{Filtered}{ Only past and the current measurements are used for
      the state estimate at time t.}
    
    \item{Smoothed}{ All measurements (both past and future) are used to
      form the state estimate at time t. This is usually the prefered type
      of state estimate.}
  }
  If \code{subsample}>0 then the data is automatically subsampled to
  provide estimated of the model states between observation time points.
  
}

\value{
    An unnamed list with one element for each individual. Each element
    contains the following elements:
      \item{Time}{Possibly subsampled time-vector corresponding
	to the estimated states }
      \item{Xs, Ps}{Smoothed state and state co-variance
	estimate}
      \item{Ys}{Response based on smoothed state: Ys = g(Xs). }
      \item{Xf, Pf}{Filtered state and state co-variance
	estimate} 
      \item{Xp, Pp}{Predicted state and state co-variance
	estimate} 
      \item{Yp, R}{Predicted observations and observation
	variances}
      \item{eta}{Estimated eta}
      \item{etaSE}{Standard errors of eta}
      \item{negLogL}{Value of the negative log-likelihood function at
	\code{THETA} (thus same value for all individuals).} 


}

\note{
  For further details please also read the package vignette pdf-document
  by writing \code{vignette("PSM")} in R.
}

\references{Please visit \url{http://www.imm.dtu.dk/psm} or refer to the
  help page for \code{\link{PSM}}. } 

\author{ Stig B. Mortensen, Soeren Klim, and Robert Miller}

\seealso{   \code{\link{PSM}}, \code{\link{PSM.estimate}},
  \code{\link{PSM.simulate}}, \code{\link{PSM.plot}}, \code{\link{PSM.template}}
}
\examples{
#detailed examples are provided in the package vignette

#Theophylline data from Boeckmann et al (1994)
#objective: recover the administered doses

library(datasets)
data(Theoph) 


#reshape data to PSM format

TheophPSM = list()
for(i in 1:length(unique(Theoph$Subject))){
  TheophPSM[[i]] = with(
    Theoph[Theoph$Subject == i,],
    list(Y = matrix(conc, nrow=1), Time = Time)
  )
}


#specify a simple pharmacokinetic model comprised of
#2 state equations and 1 observation equation
#initial value of 1 state eq. varies randomly across individuals

mod = vector(mode="list")
mod$Matrices = function(phi) {
  list(
    matA=matrix(c(-phi$ka, 0, phi$ka, -phi$ke), nrow=2, ncol=2, byrow=TRUE),
    matC=matrix(c(0, 1), nrow=1, ncol=2)
  )
}
mod$h = function(eta, theta, covar) {
  phi = theta
  phi$dose = theta$dose * exp(eta[1])
  phi
}
mod$S = function(phi) {
  matrix(c(phi$sigma), nrow=1, ncol=1)
}
mod$SIG = function(phi) {
  matrix(c(0, 0, 0, phi$omega), nrow=2, ncol=2, byrow=TRUE)
}
mod$X0 = function(Time, phi, U) {
  matrix(c(phi$dose, 0), nrow=2, ncol=1)
}
mod$ModelPar = function(THETA) {
  list(theta=list(dose = THETA["dose"], ka = THETA["ka"], ke = THETA["ke"],
                  omega = THETA["omega"], sigma = THETA["sigma"]),
       OMEGA=matrix(c(THETA["BSV_dose"]), nrow=1, ncol=1)
  )
}


#specify the search space of the fitting algorithm

parM = c(ka = 1.5, ke = 0.1, dose = 10, omega = .3, sigma = 1,
         BSV_dose = 0.015)
pars = list(LB=parM*.25, Init=parM, UB=parM*1.75)

#fit model and predict data

fit = PSM.estimate(mod, TheophPSM, pars, trace = 1, fast = TRUE,
  control=list(optimizer="optim", maxit=1))
pred = PSM.smooth(mod, TheophPSM, fit$THETA)

#visualize recovery performance

true_dose = tapply(Theoph$conc, Theoph$Subject, mean)
true_dose = true_dose[order(as.numeric(names(true_dose)))]
est_dose = fit$THETA["dose"] * exp(unlist(lapply(pred, function(x) x$eta)))
plot(true_dose, est_dose,
  xlab="actually administered dose", ylab= "recovered dose")
abline(lm(est_dose ~ true_dose), lty=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}        %Statistical Inference
\keyword{models}       %Statistical Models
\keyword{multivariate} %Multivariate Techniques
\keyword{ts}           %Time Series

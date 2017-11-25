\name{PSM.simulate}
\alias{PSM.simulate}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Create simulation data for multiple individuals}
\description{
  Simulates data for multiple individuals in a mixed effects model based
  on stochastic differential equations using an euler scheme.
}
\usage{
PSM.simulate(Model, Data, THETA, deltaTime, longX=TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Model}{
    A list containing the model components either Linear or Non-Linear
    Model list.*
  }
  \item{Data}{
    List with elements described below. No \code{Data$Y} is needed as it
    is generated through the simulation. The number of individuals
    simulated is equal to length(Data).
    \describe{
      \item{\code{Time}}{
	Time vector
      }
      \item{\code{U}}{
	Input list for the Model
      }
      \item{\code{covar}}{
	Covariates list 
      }
    }
  }
  \item{THETA}{
    Vector of population parameters
  }
  \item{deltaTime}{
    Time Step in the Euler scheme
  }
  \item{longX}{
    Boolean. Toggles output of the entire simulated outcome of the states
  }
  * See description in \link{PSM.estimate}.
}
\details{
  
  The \eqn{\eta_i}{eta} is drawn from the multivariate normal
  distribution \eqn{N(0,\Omega)}{N(0,OMEGA)}. The simulation is an euler
  based method but for every time interval dt the model is predicted and
  the states affected by system noise (\eqn{\sigma}{SIG}).
  
  The measurements are added an normal error term belonging to \eqn{N(0,S)}.
  
  The function \code{mvrnorm} from the MASS pacakge is used to to
  generate random numbers fra multivariate normal distributions.  

}
\value{

  The simulated outcome of the model is returned in a list, where each
  element is the data for an individual.
  
  \item{X}{
    Simulated states sampled at time points for measurements
  }
  \item{Y}{
    Simulated measurements
  }
  \item{Time}{
    Time points for measurements
  }
  \item{U}{
    Input vector used in the simulation
  }
  \item{eta}{
    The random effects used in the simulation
  }
  \item{Dose}{
    The dose list used in the simulation
  }
  \item{longX}{
    Entire outcome of simulated states
  }
  \item{longTime}{
    Time points for \code{longX}.
  }
}

\note{
  For further details please also read the package vignette pdf-document
  by writing \code{vignette("PSM")} in R.
}

\references{ Please visit \url{http://www.imm.dtu.dk/psm} or refer to
  the help page for \code{\link{PSM}}. } 

\author{ Stig B. Mortensen and Soeren Klim }

\seealso{  \code{\link{PSM}}, \code{\link{PSM.estimate}},
  \code{\link{PSM.smooth}}, \code{\link{PSM.plot}}, \code{\link{PSM.template}}
}

\examples{
#specify pharmacokinetic model
#2 state equations, 1 observation equation, 1 random effect

mod = vector(mode="list")
mod$Matrices = function(phi) {
  list(
    matA=matrix(c(-phi$ka, 0, phi$ka, -phi$ke), nrow=2, ncol=2, byrow=T),
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
  matrix(c(0, 0, 0, phi$omega), nrow=2, ncol=2, byrow=T)
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


#specify sampling scheme and RNG

TheophPSM <- list()
TheophPSM[[1]] <- list(Time = seq(0,25,5))
set.seed(12345)

\dontrun{#simulate and visualize ODE model (no volatility)

parM <- c(ka = 1.58, ke = 0.08, dose = 9.54, omega = 0, sigma = 1.05, BSV_dose = 0)
TheophSim <- PSM.simulate(mod, TheophPSM, THETA = parM, deltaTime = 0.01)
plot(TheophSim[[1]]$longTime, TheophSim[[1]]$longX[2,],
     type="l", ylab="concentration", xlab="time")


#contrast it to SDE model

parM <- c(ka = 1.58, ke = 0.08, dose = 9.54, omega = 0.34, sigma = 1.05, BSV_dose = 0)
TheophSim <- PSM.simulate(mod, TheophPSM, THETA = parM, deltaTime = 0.01)
lines(TheophSim[[1]]$longTime, TheophSim[[1]]$longX[2,],
     ylab="concentration", xlab="time")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}        %Statistical Inference
\keyword{models}       %Statistical Models
\keyword{multivariate} %Multivariate Techniques
\keyword{ts}           %Time Series

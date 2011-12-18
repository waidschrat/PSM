\name{PSM-package}
\alias{PSM}
\docType{package}
\title{
  Population Stochastic Modelling
}
\description{
  Mixed-effects models using Stochastic Differential Equations
  
  This package provides functions for estimation and simulation of
  multivariate linear and non-linear mixed-effects models using stochastic
  differential equations. The package allows for multidimensional input,
  specification of dosage regimen and is able to return smoothed
  estimates of model states.
}

\details{
  Function overview:
  
  \code{PSM.estimate}\cr
  Estimate population parameters for any linear or non-linear model.
  
  \code{PSM.smooth}\cr
  Optimal estimates of model states based on estimated parameters.
  
  \code{PSM.simulate}\cr
  Simulate data for multiple individuals.

  \code{PSM.plot}\cr
  Plot data, state estimates ect. for multiple individuals.

  \code{PSM.template}\cr
  Creates a template with R-syntax to help setup a model in PSM.

}

\note{
  For further details please also read the package vignette pdf-document
  by writing \code{vignette("PSM")} in R.
}

\author{
  Stig B. Mortensen and Søren Klim
  
  Maintainer: Søren Klim <skl@imm.dtu.dk>
}

\references{
  Stig B. Mortensen, Søren Klim, Bernd Dammann, Niels R. Kristensen,
  Henrik Madsen, Rune V. Overgaard. A matlab framework for estimation of
  NLME models using stochastic differential equations: Application for
  estimation of insulin secretion rates. J Pharmacokinet Pharmacodyn
  (2007) 34:623-642.

  Web:  \url{http://www.imm.dtu.dk/psm}
  

}
\keyword{ package }
\seealso{
  \code{\link{PSM.estimate}}, \code{\link{PSM.smooth}},
  \code{\link{PSM.simulate}}, \code{\link{PSM.plot}}, \code{\link{PSM.template}}
}
%\examples{
%print("There are no examples.")
%}

\name{PSM-package}
\alias{PSM}
\docType{package}
\title{
  Population Stochastic Modelling
}
\description{
  Non-Linear Mixed-Effects Modelling using Stochastic Differential Equations
  
  Functions for fitting linear and non-linear mixed-effects models using
  stochastic differential equations (SDEs). The package allows for any multivariate
  non-linear time-variant model to be specified, and it also handles multidimensional
  input, covariates, missing observations, and specification of dosage regimen.
  The provided pipeline relies on the coupling of the FOCE algorithm and Kalman filtering
  as outlined by Klim et al. (2009) and has been validated against the proprietary
  software 'NONMEM' (Tornoe et al., 2005). Further functions are provided for finding smoothed
  estimates of model states and for simulation.
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
  Stig B. Mortensen <stigbm@gmail.com>, Soeren Klim <soren@klimens.dk>, and Robert Miller <robert.miller@tu-dresden.de>
}

\references{
  Klim, S., Mortensen, S. B., Kristensen, N. R., Overgaard, R. V., & Madsen, H. (2009). Population stochastic
  modelling (PSM)—an R package for mixed-effects models based on stochastic differential equations. Computer
  Methods and Programs in Biomedicine, 94:279-289.
  
  Moler, C., & Van Loan, C. (2003). Nineteen dubious ways to compute the exponential of a matrix, twenty-five
  years later. SIAM Review, 45:3-49.
  
  Mortensen, S. B., Klim, S., Dammann, B., Kristensen, N. R., Madsen, H., Overgaard, R. V. (2007).
  A matlab framework for estimation of NLME models using stochastic differential equations: Application for
  estimation of insulin secretion rates. Journal of Pharmacokinetics and Pharmacodynamics, 34:623-642.
  
  Tornoe, C. W., Overgaard, R. V., Agersoe, H., Nielsen, H. A., Madsen, H., & Jonsson, E. N. (2005). Stochastic
  differential equations in NONMEM: implementation, application, and comparison with ordinary differential
  equations. Pharmaceutical Research, 22:1247-1258.
  
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

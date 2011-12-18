\name{PSM.estimate}
\alias{PSM.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate population parameters}
\description{
  Estimates population parameters in a linear or non-linear mixed
  effects model based on stochastic differential equations by use of
  maximum likelihood and the Kalman filter.
}
\usage{
PSM.estimate(Model, Data, Par, CI = FALSE, trace = 0, control=NULL, fast=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Model}{
    A list containing the following elements:

    \describe{
      \item{\code{Matrices = function(phi)}}{
	\emph{Only in linear models.}

	Defines the matrices \eqn{A}, \eqn{B},
	\eqn{C} and \eqn{D} in
	the model equation. Must return a list of matrices named
	\code{matA}, \code{matB}, \code{matC} and \code{matD}.
	If there is no input, \code{matB} and \code{matD} may be omitted
	by setting them to \code{NULL}. Note, if the matrix \eqn{A} is
	singular the option \code{fast} is set to \code{FALSE}, as
	this is not supported in the compiled Fortran code.
      }
      \item{\code{Functions}}{
	\emph{Only in non-linear models.}

	A list containing the functions
	\code{f(x,u,time,phi)}, \code{g(x,u,time,phi)},
	\code{df(x,u,time,phi)} and \code{dg(x,u,time,phi)}.

	The functions \code{f} and \code{g} defines the
	system and \code{df}
	and \code{dg} are the Jacobian matrices with first-order partial
	derivatives for \code{f(x)} and \code{g(x)} which is needed to
	evaluate the model. A warning is issued if \code{df}
	or \code{dg} appear to be incorrect based on a numerical
	evaluation of the Jacobians of \code{f(x)} and \code{g(x)}.

	It is possible to avoid specifying the Jacobian functions in the model
	and use numerical approximations instead, but this will increase
	estimation time at least ten-fold. See the section
	\sQuote{Numerical Jacobians of f and g} below for more information.
      }
      \item{\code{X0 = function(Time, phi, U)}}{
	Defines the model state at \code{Time[1]} before
	update. \code{Time[1]} and \code{U[,1]} can be used in the
	evaluation of X0. Must return a column matrix.
      }
      \item{\code{SIG = function(phi)}}{ in linear models and \code{SIG
	  = function(u,time,phi)} in non-linear models. It defines the
	matrix \eqn{\sigma}{SIG} for the diffusion 
	term. Returns a square matrix. 
      }
      \item{\code{S = function(phi)}}{ in linear models and \code{S =
	  function(u,time,phi)} in non-linear models. It defines a
	covariance matrix for the observation noise. Returns a 
	square matrix. 
      }
      \item{\code{h = function(eta,theta,covar)}}{
	Second stage model. Defines how random effects (\code{eta}) and
	covariates (\code{covar}) affects the fixed effects parameters
	(\code{theta}). In models where \code{OMEGA=NULL} (no random-effects)
	\code{h} must still be defined with the same argument list to
	allow for covariates to affect theta, but the function \code{h}
	is evaluated with \code{eta=NULL}.
	Must return a list (or vector) \code{phi} of individual
	parameters which is used as input argument in the other user-defined functions. 
      }
      \item{\code{ModelPar = function(THETA)}}{
	Defines the population parameters to be optimized. Returns a
	list containing 2 elements, named: 
	\describe{
	  \item{\code{theta}}{A list of fixed effects parameters
	    \eqn{\theta} which are used as input to the function
	    \code{h} listed above.} 
	  \item{\code{OMEGA}}{A square covariance matrix
	    \eqn{\Omega}{OMEGA} for the random effects. If
	    \code{OMEGA} is missing or \code{NULL} then no 2nd stage
	    model is used. However, the function \code{h} must still be
	    defined, see above. }
	}
      }
    }
  }
  \item{Data}{An unnamed list where each element contains
    data for one individual. Each element in \code{Data} is a list
    containing: 

    \describe{
      \item{\code{Time}}{A vector of timepoints for measurements}

      \item{\code{Y}}{A matrix of multivariate observations for each timepoint, where
	each column is a multivariate measurement. Y may contain \code{NA} for missing
	observations and a column may consist of both some or only
	\code{NA}s. The latter is useful if a dose is given when no measurement
	is taken. }

      \item{\code{U}}{A matrix of multivariate input to the model for each
	timepoint. \code{U} is assumed constant between measurements and may not
	contain any \code{NA}. If \code{U} is ommitted, the model is 
	assumed to have no input and \code{matB} and \code{matD} need no
	to be specified.}
      
      \item{\code{Dose}}{
	A list containing the 3 elements listed below. If the element
	\code{Dose} is missing or \code{NULL}, no dose is assumed. 
	\describe{
	  \item{\code{Time}}{ A vector of timepoints for the dosing. Each must
	    coinside with a measurement time. Remember to insert a
	    missing measurement in \code{Y} if a corresponding
	    timepoint is not present. Dose is considered added to the system
	    just after the measurement.}
	  \item{\code{State}}{ A vector with indexes of the state for dosing.}
	  \item{\code{Amount}}{ A vector of amounts to be added.}
	}
      }
    }
  }
  \item{Par}{A list containing the following elements:

    \describe{
      \item{\code{Init}}{A vector with initial estimates for \code{THETA}, vector of
	population parameters to be optimized.}

      \item{\code{LB}, \code{UB}}{: Two vectors with lower and upper bounds for
	parameters. If ommitted, the program performs unconstrained
	optimization. It is highly recommended to specify bounds to ensure
	robust optimization. }
    }
  }
  \item{CI}{ Boolean. If true, the program estimates 95\% confidence
    intervals, standard deviation and correlation matrix for the
    parameter estimates based on the Hessian of the likelihood function. The
    Hessian is estimated by \code{hessian} in the \code{numDeriv} package.
  }
  \item{trace}{ Non-negative integer. If positive, tracing
    information on the progress of the optimization is produced. Higher
    values produces more tracing information.
  }
  \item{control}{ A list of control parameters for the optimization of
    the likelihood function. The list has one required component, namely:
    \describe{
      \item{\code{optimizer}}{ A string value equal to either
	\code{'optim'} or \code{'ucminf'}. This gives the choise of optimizer. Default
	is \code{optimizer = 'optim'}.
      }
    }
    The remaining components in the list are given as the control argument
    for the chosen optimizer. See corresponding help file for further detail.
  }
  \item{fast}{  Boolean. Use compiled Fortran code for faster
    estimation.
  }
}

\details{
  The first stage model describing intra-individual variations is
  for linear models defined as

  \deqn{ dx_t = (A(\phi_i)x_t + B(\phi_i)u_t)dt + \sigma(\phi_i) d\omega_t}{ dx
    = (A(phi)*x + B(phi)*u)dt + SIG(phi)*dw} 
  \deqn{ y_{ij} = C(\phi_i)x_{ij} + D(\phi_i)u_{ij} + e_{ij}}{y = C(phi)*x + D(phi)*u + e}

  and for  non-linear models as

  \deqn{ dx_t = f(x_t,u_t,t,\phi_i)dt + \sigma(u_t,t,\phi_i) d\omega_t}{
    dx = f(x,u,t,phi)dt + SIG(u,t,phi)dw}
  \deqn{y_{ij} = g(x_{ij},u_{ij},t_{ij},\phi_i) + e_{ij}}{y = g(x, u, t,
    phi) + e}

  where \eqn{e_{ij} \sim N(0,S(u_{ij},t_{ij},\phi_i))}{e ~ N(0,S(x, u, t))} and \eqn{\omega_t}{w} is a
  standard Brownian motion. 

  The second stage model describing inter-individual variations is
  defined as:

  \deqn{\phi_i = h(\eta_i,\theta,Z_i)}{phi = h(eta,theta,Z)}

  where \eqn{\eta_i \sim N(0,\Omega)}{eta ~ N(0,OMEGA)}, \eqn{\theta}
  are the fixed effect parameters and \eqn{Z_i}{Z} are covariates for
  individual i. In a model without random-effects the function \eqn{h}
  is only used to include possible covariates in the model.
  
}
\value{
  A list containing the following elements:
  
  \item{NegLogL}{Value of the negative log-likelihood function at optimum.}
  \item{THETA}{Population parameters at optimum}
  \item{CI}{95\% confidence interval for the estimated parameters}
  \item{SD}{Standard deviation for the estimated parameters}
  \item{COR}{Correlation matrix for the estimated parameters}
  \item{sec}{Time for the estimation in seconds}
  \item{opt}{Raw output from \code{optim}}
}

\section{Numerical Jacobians of f and g}{
  Automatic numerical approximations of the Jacobians of \code{f} and
  \code{g} can be used in PSM. In the folliwing, the name of the model
  object is assumed to be \code{MyModel}.

  First define the
  functions \code{MyModel$Functions$f} and
  \code{MyModel$Functions$g}. When these are defined in MyModel the
  functions \code{df} and \code{dg} can be added to the model object by
  writing as below: 

  \preformatted{
  MyModel$Functions$df = function(x,u,time,phi) {
    jacobian(MyModel$Functions$f,x=x,u=u,time=time,phi=phi)
  }
  MyModel$Functions$dg = function(x,u,time,phi) {
    jacobian(MyModel$Functions$g,x=x,u=u,time=time,phi=phi)
  }

  }

  This way of defining \code{df} and \code{dg} forces a numerical
  evaluation of the Jacobians using the \pkg{numDeriv} package. It may
  be usefull in some cases, but it should be stressed that it will
  probably give at least a ten-fold increase in estimation times.

}

\note{
  For further details please also read the package vignette pdf-document
  by writing \code{vignette("PSM")} in R.
}

\references{ Please visit \url{http://www.imm.dtu.dk/psm} or refer to
  the main help page for \code{\link{PSM}}. }

\author{ Stig B. Mortensen and Søren Klim}

\seealso{
  \code{\link{PSM}}, \code{\link{PSM.smooth}},
  \code{\link{PSM.simulate}}, \code{\link{PSM.plot}}, \code{\link{PSM.template}}
} 

\examples{
cat("\nExamples are included in the package vignette.\n")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}        %Statistical Inference
\keyword{models}       %Statistical Models
\keyword{multivariate} %Multivariate Techniques
\keyword{ts}           %Time Series

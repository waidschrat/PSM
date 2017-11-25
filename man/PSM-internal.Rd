\name{Internal functions}
\alias{APL.KF}
\alias{APL.KF.gr}
\alias{APL.KF.individualloop}
\alias{CutThirdDim}
\alias{ExtKalmanFilter}
\alias{ExtKalmanSmoother}
\alias{IndividualLL.KF}
\alias{IndividualLL.KF.gr}
\alias{LinKalmanFilter}
\alias{LinKalmanSmoother}
\alias{ModelCheck}
\alias{logit}
\alias{invlogit}


\title{Internal functions in the PSM-package.}
\description{
  Internal functions in the PSM-package.
}
\usage{
APL.KF(THETA, Model, Pop.Data, LB = NULL, UB = NULL, GUIFlag = 0, longOutput = FALSE,
       fast=TRUE,Linear=NULL)
APL.KF.gr(THETA, Model, Pop.Data, LB = NULL, UB = NULL, GradSTEP = 1e-04, GUIFlag = 0,
          fast=TRUE,Linear=NULL)
APL.KF.individualloop(theta, OMEGA, Model, Data, GUIFlag = 0, fast=TRUE,Linear)
CutThirdDim(a)
ExtKalmanFilter(phi, Model, Data, outputInternals = FALSE)
ExtKalmanSmoother(phi, Model, Data)
IndividualLL.KF(eta, theta, OMEGA, Model, Data, fast=TRUE,Linear=NULL)
IndividualLL.KF.gr(eta, theta, OMEGA, Model, Data, GradSTEP = 1e-04, GUIFlag = 0,
                   fast=TRUE,Linear=NULL)
LinKalmanFilter(phi, Model, Data, echo = FALSE, outputInternals = FALSE, fast=TRUE)
LinKalmanSmoother(phi, Model, Data)
ModelCheck(Model, Data, Par,DataHasY=TRUE)
logit(x, xmin, xmax)
invlogit(y, xmin, xmax)
}
%- maybe also 'usage' for other objects documented here.

\details{
  \describe{
    \item{\code{APK.KF}}{
      evaluates the population likelihood function.
    }
    \item{\code{APK.KF.gr}}{
      evaluates the gradient of APL.KF.
    }
    \item{\code{APL.KF.individualloop}}{
      contains the innner loop over individuals for APL.KF.
    }
    \item{\code{CutThirdDim}}{
      removes third and higher dimensions of
  dim-attribute for an array and thus creating a matrix.
    }
    \item{\code{ExtKalmanFilter}}{
      Performs a Extended Kalman filtering.
    }
    \item{\code{ExtKalmanSmoother}}{
      performs a non-linear Kalman smoothing.
    }
    \item{\code{IndividualLL.KF}}{
      evaluates the indivdual neg. log-likelihood function.
    }
    \item{\code{IndividualLL.KF.gr}}{
      evaluates the gradient of the indivdual neg. log-likelihood function.
    }
    \item{\code{LinKalmanFilter}}{
      performs a linear Kalman filtering.
    }
    \item{\code{LinKalmanSmoother}}{
      performs a linear Kalman smoothing.
    }
    \item{\code{ModelCheck}}{
      checks for dimensionalities and model objects. Furthermore it
      tests the Model objects and the dimensions in the Data set.
    }
    \item{\code{logit}}{
      gives logit transformation of a vector.
    }
    \item{\code{invlogit}}{
      gives invlogit transformation of a vector.
    }
  }
}

\author{ Stig B. Mortensen and Soeren Klim }
\seealso{ \code{\link{PSM}}}
%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}

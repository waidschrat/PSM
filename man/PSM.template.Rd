\name{PSM.template}
\alias{PSM.template}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Creates a template for a model in PSM}
\description{
  Creates a template with R-syntax to help setup a model in PSM. 
}
\usage{
  PSM.template(Linear=FALSE,dimX=2,dimY=3,dimU=4,dimEta=5,file="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Linear}{
Boolean. Linear or non-linear model.
  }
  \item{dimX}{
Number of state equations.
  }
  \item{dimY}{
Number of response variables.
  }
  \item{dimU}{
Number of input variables (can be zero).
  }
  \item{dimEta}{
Number of random effects (can be zero).
  }
  \item{file}{
  A character string naming the file to print
  to.  If '""' (the default), \code{PSM.template} prints to the standard
 output connection.
  }

}

\value{
  None (invisible \code{NULL}).
}

\note{
  For further details please also read the package vignette pdf-document
  by writing \code{vignette("PSM")} in R.
}

\references{ Please visit \url{http://www.imm.dtu.dk/psm} or refer to
  the help page for \code{\link{PSM}}. } 

\author{ Stig B. Mortensen and Søren Klim }

\seealso{  \code{\link{PSM}}, \code{\link{PSM.estimate}},
  \code{\link{PSM.smooth}}, \code{\link{PSM.template}}
}

\examples{
# Linear model with input, random effects and dose
PSM.template(Linear=TRUE,dimX=1,dimY=2,dimU=3,dimEta=4)

# Non-linear model without input, random effects and dose
PSM.template(Linear=FALSE,dimX=1,dimY=2,dimU=0,dimEta=0)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}        %Statistical Inference
\keyword{models}       %Statistical Models
\keyword{multivariate} %Multivariate Techniques
\keyword{ts}           %Time Series

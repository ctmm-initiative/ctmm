\name{overlap}
\alias{overlap}
\alias{overlap.ctmm}
\alias{overlap.telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Calculate the overlap between two stationary distributions}

\description{This function calculates a useful measure of similarity between distributions known as the \emph{Bhattacharyya coefficient} in statistics and simply the \emph{fidelity} or \emph{overlap} in quantum and statistical mechanics. It is roughly speaking the ratio of the intersection area to the average individual area. When applied to \code{ctmm} objects, this function returns the overlap of the two Gaussian distributions. When applied to \code{telemetry} objects, this function returns the overlap of their kernel density estimates.}

\usage{
overlap(object1,object2,...)

\S3method{overlap}{ctmm}(object1, object2, level=0.95, ...)

\S3method{overlap}{telemetry}(object1, object2, CTMM1, CTMM2, level=0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object1}{A \code{ctmm} fit or \code{telemetry} object for one animal.}
  \item{object2}{A \code{ctmm} fit or \code{telemetry} object for another animal.}
  \item{level}{The confidence level desired for the output.}
  \item{CTMM1}{The \code{ctmm} fit corresponding to \code{telemetry} \code{object1}.}
  \item{CTMM2}{The \code{ctmm} fit corresponding to \code{telemetry} \code{object2}.}
  \item{...}{Additional arguments relevant for \code{\link{akde}}, such as \code{debias} and \code{res}.}
}

%\details{}

\value{Confidence intervals on the overlap estimate. A value of \code{1} implies that the two distributions are identical, while a value of \code{0} implies that the two distributions share no area in common.}

%% ~put references to the literature/web site here ~
%\references{}

\author{C.H.Fleming}

\note{Uncertainties in \code{CTMM1} and \code{CTMM2} are propagated into the overlap estimate under the approximation that the Bhattacharyya distance is a chi-square random variable.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# Fit a continuous-velocity model with tau ~ c(10 days,1 hour)
# also see help(variogram.fit)
GUESS <- ctmm(tau=c(10*24*60^2,60^2))
FIT1 <- ctmm.fit(buffalo[[1]],GUESS)
FIT2 <- ctmm.fit(buffalo[[2]],GUESS)

# Gaussian overlap between these two buffalo
overlap(FIT1,FIT2)
}}

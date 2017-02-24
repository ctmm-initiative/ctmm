\name{overlap}
\alias{overlap}
%\alias{overlap.ctmm}
%\alias{overlap.telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Calculate the overlap between two stationary distributions}

\description{This function calculates a useful measure of similarity between distributions known as the \emph{Bhattacharyya coefficient} in statistics and simply the \emph{fidelity} or \emph{overlap} in quantum and statistical mechanics. It is roughly speaking the ratio of the intersection area to the average individual area. When applied to \code{ctmm} objects, this function returns the overlap of the two Gaussian distributions. When applied to \code{telemetry} or (aligned) \code{UD} objects with corresponding movement models, this function returns the overlap of their (autocorrelated) kernel density estimates.}

\usage{ overlap(object,CTMM=NULL,level=0.95,...) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{list} of \code{ctmm} fit, \code{telemetry}, or aligned \code{UD} objects to compare.}
  \item{CTMM}{A \code{list} of \code{ctmm} fit objects corresponding to \code{object} list.}
  \item{level}{The confidence level desired for the output.}
  \item{...}{Additional arguments relevant for \code{\link{akde}}, such as \code{res} and \code{weights}.}
}

%\details{}

\value{A table of confidence intervals on the overlap estimate. A value of \code{1} implies that the two distributions are identical, while a value of \code{0} implies that the two distributions share no area in common. \code{ctmm} objects are necessary to provide confidence intervals on the point esitmate.}

%% ~put references to the literature/web site here ~
%\references{}

\author{C. H. Fleming and K. Winner}

\note{Uncertainties in \code{CTMM1} and \code{CTMM2} are propagated into the overlap estimate under the approximation that the Bhattacharyya distance is a chi-square random variable.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{ctmm.fit}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# Fit a continuous-velocity model with tau ~ c(10 days,1 hour)
# also see help(variogram.fit)
GUESS <- ctmm(tau=c(10*24*60^2,60^2))
FITS <- list()
FITS[[1]] <- ctmm.fit(buffalo[[1]],GUESS)
FITS[[2]] <- ctmm.fit(buffalo[[2]],GUESS)
names(FITS) <- names(buffalo[1:2])

# Gaussian overlap between these two buffalo
overlap(FITS)

# AKDE overlap between these two buffalo
overlap(buffalo[1:2],FITS)
}}

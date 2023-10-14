\name{distance}
\alias{distance}
\encoding{UTF-8}
%\alias{overlap.ctmm}
%\alias{overlap.telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Calculate the square distance between two distributions or location estimates}

\description{This function calculates various square distances measures between distributions,
including the, Bhattacharyya distance, Mahalanobis distance, and Euclidean distance.}

\usage{ distance(object,method="Mahalanobis",sqrt=FALSE,level=0.95,debias=TRUE,...) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{list} of \code{ctmm} fit objects or single-location \code{telemetry} objects to compare.}
  \item{method}{Square distance measure to return: \code{"Bhattacharyya"}, \code{"Mahalanobis"}, or \code{"Euclidean"}.}
  \item{sqrt}{Return the linear distance.}
  \item{level}{The confidence level desired for the output.}
  \item{debias}{Approximate debiasing of the square distance.}
  \item{...}{Not currently used.}
}

%\details{}

\value{A list with tables \code{DOF}, containing the effective samples sizes of the estimates, and \code{CI}, containing the confidence intervals of the distance estimates. A value of \code{0} implies that the two distributions have the same mean location, while larger values imply that the two distributions are farther apart. The (square) Euclidean distance has units of square meters, if \code{sqrt=FALSE}. The square Mahalanobis and Bhattacharyya distances are unitless. For the Euclidean distance, only the centroids are compared (in meters if \code{sqrt=TRUE} or square meters if \code{sqrt=FALSE}).
}

%% ~put references to the literature/web site here ~
%\references{}

\author{C. H. Fleming}

\note{The Bhattacharyya distance (BD) is naturally of a squared form and is not further squared.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}}, \code{\link{overlap}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# fit models for first two buffalo
GUESS <- lapply(buffalo[1:2], function(b) ctmm.guess(b,interactive=FALSE) )
# using ctmm.fit here for speed, but you should almost always use ctmm.select
FITS <- lapply(1:2, function(i) ctmm.fit(buffalo[[i]],GUESS[[i]]) )
names(FITS) <- names(buffalo[1:2])

# Mahalanobis distance between these two buffalo
distance(FITS)
}}

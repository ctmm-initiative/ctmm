\name{distance}
\alias{distance}
\encoding{UTF-8}
%\alias{overlap.ctmm}
%\alias{overlap.telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Calculate the square distance between two stationary distributions}

\description{This function calculates various square distances measures between distributions,
including the, Bhattacharyya distance, Mahalanobis distance, and Euclidean distance.}

\usage{ distance(object,method="Mahalanobis",level=0.95,debias=TRUE,...) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{list} of \code{ctmm} fit objects to compare.}
  \item{method}{Square distance measure to return: \code{"Bhattacharyya"}, \code{"Mahalanobis"}, or \code{"Euclidean"}.}
  \item{level}{The confidence level desired for the output.}
  \item{debias}{Approximate debiasing of the square distance.}
  \item{...}{Not currently used.}
}

%\details{}

\value{A table of confidence intervals on the square distance estimate. A value of \code{0} implies that the two distributions have the same mean location, while larger values imply that the two distributions are farther apart. The square Euclidean distance has units of square meters. The square Mahalanobis and Bhattacharyya distances are unitless.
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

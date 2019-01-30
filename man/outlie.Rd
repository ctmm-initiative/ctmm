\name{outlie}
\alias{outlie}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods to facilitate outlier detection.}
\description{Produces a \code{data.frame} of speed and distance estimates to analyze, as well as a plot highlighting potential speed and distance outliers in \code{telemetry} data.
}
\usage{
outlie(data,UERE=10,standardize=FALSE,plot=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{\code{telemetry} object.}
  \item{UERE}{Device-dependent telemetry error in meters. Only necessary for uncalibrated data.}
  \item{standardize}{Standardize speed and distance output by their median absolute deviation (MAD).}
  \item{plot}{Output a plot highlighting high speeds (blue) and distant locations (red).}
  \item{...}{Arguments passed to \code{plot}.}
}

\details{Distances are calculated from the median longitude & latitude, while speeds are calculated over the timesteps. Both estimates account for telemetry error and condition on as few data points as possible. The speed estimates furthermore account for timestamp truncation and assign each timestep's speed to the most likely offending time, based on its other adjacent speed estimate.

If \code{plot=TRUE}, intervals of high speed are highlighted with blue segments, while distant locations are highlighted with red points.}

\value{Returns a \code{data.frame} of distances and speeds. Can also produce a plot as a side effect.}

%\references{}

\author{C. H. Fleming.}

\note{The speed estimates here are tailored for outlier detection and have poor statistical efficiency. The \code{\link{predict}} and \code{\link{speed}} methods are appropriate for estimating speed (after outliers have been removed and a movement model has been selected). }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{as.telemetry}}.  }

\examples{
# Load package and data
library(ctmm)
data(buffalo)

# look for outliers in Cilla
OUT <- outlie(buffalo[[1]])

# look at the distribution of speed esitmates
hist(OUT$speed)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

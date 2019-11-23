\name{outlie}
\alias{outlie}
\alias{outlier}
\alias{plot.outlie}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods to facilitate outlier detection.}
\description{Produces a \code{data.frame} of speed and distance estimates to analyze, as well as a plot highlighting potential speed and distance outliers in \code{telemetry} data.
}
\usage{
outlie(data,UERE=10,standardize=FALSE,plot=TRUE,...)

\S3method{plot}{outlie}(x,level=0.95,units=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{\code{telemetry} object.}
  \item{UERE}{Device-dependent telemetry error in meters (RMS UERE). Only necessary for uncalibrated data.}
  \item{standardize}{Standardize speed and distance output by their median absolute deviation (MAD).}
  \item{plot}{Output a plot highlighting high speeds (blue) and distant locations (red).}
  \item{...}{Arguments passed to \code{plot}.}
  \item{x}{\code{outlie} object to plot.}
  \item{level}{Confidence level for error bars.}
  \item{units}{Convert axes to natural units.}
}

\details{If \code{plot=TRUE} in \code{outlie()}, intervals of high speed are highlighted with blue segments, while distant locations are highlighted with red points.

When plotting the \code{outlie} object itself, `core deviation' denotes distances from the median longitude & latitude, while `minimum speed' denotes the minimum speed required to explain the location estimate's displacement as straight-line motion. Both estimates account for telemetry error and condition on as few data points as possible. The speed estimates furthermore account for timestamp truncation and assign each timestep's speed to the most likely offending time, based on its other adjacent speed estimate.

The output \code{outlie} object contains the above noted speed and distance estimates in a \code{data.frame}, with rows corresponding to those of the input \code{telemetry} object.
}

\value{Returns a \code{data.frame} of distances and speeds. Can also produce a plot as a side effect.}

%\references{}

\author{C. H. Fleming.}

\note{The speed estimates here are tailored for outlier detection and have poor statistical efficiency. The \code{\link{predict}} and \code{\link{speed}} methods are appropriate for estimating speed (after outliers have been removed and a movement model has been selected). }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{as.telemetry}}.  }

\examples{
# Load package and data
library(ctmm)
data(turtle)

# look for outliers in a turtle
OUT <- outlie(turtle[[3]])

# look at the distribution of estimates
plot(OUT)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

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
outlie(data,plot=TRUE,by='d',...)

\S3method{plot}{outlie}(x,level=0.95,units=TRUE,axes=c('d','v'),xlim=NULL,ylim=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{\code{telemetry} object.}
  \item{plot}{Output a plot highlighting high speeds (blue) and distant locations (red).}
  \item{by}{Color and size side-effect plot points by \code{'d'}, \code{'v'}, \code{'dz'}, \code{'vz'}, for distance from center, minimum speed, vertical distance from center, and minimum vertical speed.}
  \item{...}{Arguments passed to \code{plot}.}
  \item{x}{\code{outlie} object to plot.}
  \item{level}{Confidence level for error bars.}
  \item{units}{Convert axes to natural units.}
  \item{axes}{\eqn{x}-\eqn{y} axes to plot. Can be any of \code{'t'}, \code{'d'}, \code{'v'}, \code{'dz'}, \code{'vz'}, for time, distance from center, minimum speed, vertical distance from center, and minimum vertical speed.}
  \item{xlim}{\eqn{x}-axis plotting range in SI units.}
  \item{ylim}{\eqn{y}-axis plotting range in SI units.}
}

\details{If \code{plot=TRUE} in \code{outlie()}, intervals of high speed are highlighted with blue segments, while distant locations are highlighted with red points.

When plotting the \code{outlie} object itself, `median deviation' denotes distances from the geometric median, while `minimum speed' denotes the minimum speed required to explain the location estimate's displacement as straight-line motion. Both estimates account for telemetry error and condition on as few data points as possible. The speed estimates furthermore account for timestamp truncation and assign each timestep's speed to the most likely offending time, based on its other adjacent speed estimate.

The output \code{outlie} object contains the above noted speed and distance estimates in a \code{data.frame}, with rows corresponding to those of the input \code{telemetry} object.
}

\value{Returns an \code{outlie} object, which is a data.frame of distance and speed information. Can also produce a plot as a side effect.}

\references{
C. H. Fleming et al, ``A comprehensive framework for handling location error in animal tracking data'', bioRxiv 2020.06.12.130195 (2020) \doi{10.1101/2020.06.12.130195}.
}

\author{C. H. Fleming.}

\note{The speed estimates here are tailored for outlier detection and have poor statistical efficiency. The \code{\link{predict}} and \code{\link{speed}} methods are appropriate for estimating speed (after outliers have been removed and a movement model has been selected).

In \code{ctmm} v0.6.1 the \code{UERE} argument was deprecated. For uncalibrated data, the initial esitmates used by \code{outlie} are now generated on import and stated by \code{summary(uere(data))}. These values not be reasonable for arbitrary datasets.}

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

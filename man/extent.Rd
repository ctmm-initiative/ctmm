\name{extent}
%\alias{extent}
\alias{extent,telemetry-method}
\alias{extent,ctmm-method}
\alias{extent,UD-method}
\alias{extent,variogram-method}
\alias{extent,list-method}
\alias{extent,data.frame-method}
\alias{extent,matrix-method}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extent}
\description{Functions to calculate the \eqn{(x,y)} plotting extent (or bounding box) of various ctmm objects or list of such objects, for use when plotting multiple ctmm objects.}

\usage{%extent(x,...)
\S4method{extent}{telemetry}(x,level=1,...)

\S4method{extent}{ctmm}(x,level=0.95,level.UD=0.95,...)

\S4method{extent}{UD}(x,level=0.95,level.UD=0.95,...)

\S4method{extent}{variogram}(x,level=0.95,threshold=2,...)

\S4method{extent}{list}(x,...)

\S4method{extent}{data.frame}(x,level=1,...)

\S4method{extent}{matrix}(x,level=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{telemetry}, \code{ctmm}, or \code{UD} object. }
  \item{level}{ For \code{telemetry} objects, this is the fraction of locations bounded, according to two-sided quantiles. For \code{ctmm} and \code{UD} objects, this is confidence level for the magnitude of the utilization area circumscribed by \code{level.UD}. }
  \item{level.UD}{ Confidence level of the UD area. I.e., the 50\% core home range would be given by \code{level.UD=0.50}. }
  \item{threshold}{Limit \code{ylim} to \code{threshold} times the maximum semi-variance, even if the \code{level} confidence intervals exceed this amount.}
  \item{...}{ Optional arguments for future extensions.}
}

\details{Returns a \code{data.frame} with columns \code{x} and \code{y} with rows \code{min} and \code{max}.
See \code{vignette('akde')} for an example of \code{extent} used to plot multiple UDs on the same scale.
}

%\references{}

\author{C. H. Fleming}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot.telemetry}}, \code{\link{plot.variogram}}. }

%\examples{
%# Load package and data
%library(ctmm)
%data(buffalo)

%# Extract movement data for a single animal
%billy <- buffalo[[1]]

%# Plot the data
%plot.telemetry(billy)
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

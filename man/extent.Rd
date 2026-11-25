\name{extent}
\alias{extent}
\alias{extent.telemetry}
\alias{extent.ctmm}
\alias{extent.UD}
\alias{extent.variogram}
\alias{extent.list}
\alias{extent.data.frame}
\alias{extent.matrix}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extent}
\description{Functions to calculate the \eqn{(x,y)} plotting extent (or bounding box) of various ctmm objects or list of such objects, for use when plotting multiple ctmm objects.}

\usage{extent(x,...)

\S3method{extent}{telemetry}(x,...)

\S3method{extent}{ctmm}(x,level=0.95,level.UD=0.95,complete=FALSE,...)

\S3method{extent}{UD}(x,level=0.95,level.UD=0.95,complete=FALSE,...)

\S3method{extent}{variogram}(x,level=0.95,threshold=2,...)

\S3method{extent}{list}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{telemetry}, \code{ctmm}, or \code{UD} object. }
  \item{level}{ For \code{ctmm} and \code{UD} objects, this is confidence level for the magnitude of the utilization area circumscribed by \code{level.UD}. }
  \item{level.UD}{ Coverage level of the UD area. I.e., the 50\% core home range would be given by \code{level.UD=0.50}. }
  \item{complete}{ Also calculate longitude-latitude extent of \code{ctmm} and \code{UD} objects. }
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

\name{projection}
%\alias{projection}
\alias{projection,telemetry-method}
\alias{projection<-,telemetry-method}
\alias{projection,ctmm-method}
\alias{projection,UD-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Projection}
\description{Get the coordinate reference system (CRS) of a ctmm object}

\usage{%projection(x,...)
\S4method{projection}{telemetry}(x,asText=TRUE)

\S4method{projection}{ctmm}(x,asText=TRUE)

\S4method{projection}{UD}(x,asText=TRUE)

\S4method{projection}{telemetry}(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{telemetry}, \code{ctmm}, or \code{UD} object. }
  \item{asText}{ If \code{TRUE}, the projection is returned as text. Otherwise a \code{\link{CRS}} object is returned.}
  \item{value}{Projection to apply.}
}

\details{\code{projection(x)} returns the projection information from ctmm object \code{x}, while \code{projection(x) <- value} applies the projection \code{value} to ctmm object \code{x}.}

%\references{}

\author{C. H. Fleming}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{as.telemetry}}. }

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

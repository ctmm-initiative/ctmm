\name{speed}
\alias{speed}
\alias{speed.ctmm}
\alias{speed.telemetry}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate the average speed of a tracked animal}
\description{
Given a \code{ctmm} movement model and \code{telemetry} data, \code{speed} simulates multiple realizations of the animal's trajectory and estimates the average speed, which is proportional to distance travled. In this way, both tortuosity and telemetry error can be accounted for.
}
\usage{
speed(object,...)

\S3method{speed}{ctmm}(object,data=NULL,level=0.95,prior=TRUE,fast=TRUE,error=0.01,mc.cores=NULL,...)

\S3method{speed}{telemetry}(object,CTMM,level=0.95,prior=TRUE,fast=TRUE,error=0.01,mc.cores=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{ctmm} movement-model or \code{telemetry} object, which requires an additional \code{CTMM} argument.}
  \item{data}{Optional \code{telemetry} object on which the simulations will be conditioned.}
  \item{CTMM}{Movement model object.}
  \item{level}{Confidence level to report on the estimated average speed.}
  \item{prior}{Account for model uncertainty.}
  \item{fast}{Whether or not to invoke the central-limit theorem (see \code{\link{emulate}}).}
  \item{error}{Target (relative) standard error.}
  \item{mc.cores}{How many CPU cores to use (Unix required). Will use all if \code{NULL}.}
  \item{...}{Arguments passed to \code{\link{emulate}}.}
}
\details{
}

\value{ }

\references{
}


\author{ C. H. Fleming. }

\note{ The average speeds estimated here are mean speeds. The speeds reported by \code{\link{summary.ctmm}} are root-mean-square (RMS) speeds. These quantities are sometimes proportional, but not equivalent. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{emulate}}, \code{\link{simulate}} }

%\examples{}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{speed}
\alias{speed}
\alias{speed.ctmm}
\alias{speed.telemetry}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate the average speed of a tracked animal}
\description{
Given a \code{ctmm} movement model and \code{telemetry} data, \code{speed} simulates multiple realizations of the animal's trajectory and estimates the average speed, which is proportional to distance traveled. Both tortuosity (non straight-line motion between the data) and telemetry error can be accounted for. Given only a \code{ctmm} movement model and no data, \code{speed} calculates the average speed of the Gaussian movement process.
}
\usage{
speed(object,...)

\S3method{speed}{ctmm}(object,data=NULL,level=0.95,units=TRUE,prior=TRUE,fast=TRUE,cor.min=0.5,dt.max=NULL,
      error=0.01,cores=1,...)

\S3method{speed}{telemetry}(object,CTMM,level=0.95,units=TRUE,prior=TRUE,fast=TRUE,cor.min=0.5,dt.max=NULL,
      error=0.01,cores=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{ctmm} movement-model or \code{telemetry} object, which requires an additional \code{CTMM} argument.}
  \item{data}{Optional \code{telemetry} object on which the simulations will be conditioned.}
  \item{CTMM}{Movement model object.}
  \item{level}{Confidence level to report on the estimated average speed.}
  \item{units}{Convert result to natural units.}
  \item{prior}{Account for model uncertainty.}
  \item{fast}{Whether or not to invoke the central-limit theorem (see \code{\link{emulate}}).}
  \item{cor.min}{Velocity correlation threshold for skipping gaps.}
  \item{dt.max}{Absolute gap sizes to skip (in seconds), alternative to \code{cor.min}.}
  \item{error}{Target (relative) standard error.}
  \item{cores}{Number of simulations to run in parallel. \code{cores=NULL} will use all cores, while \code{cores<0} will reserve \code{abs(cores)}.}
  \item{...}{Arguments passed to \code{\link{emulate}}.}
}
\details{
The \code{cor.min} or \code{dt.max} arguments are used to constrain the estimate to be derived from simulations near the data, and therefore ensure that the estimate is more reflective of the data than the model.
}

\value{ Returns the estimated mean speed of the sampled trajectory with CIs by default. If \code{level=NULL}, then the ensemble of mean speeds is returned instead. }

%\references{ }

\author{ C. H. Fleming. }

\note{ The average speeds estimated here are mean speeds. The speeds reported by \code{\link{summary.ctmm}} are root-mean-square (RMS) speeds. These quantities are sometimes proportional, but not equivalent. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{emulate}}, \code{\link{simulate}} }

\examples{
# Load package and data
library(ctmm)
data(buffalo)
Cilla <- buffalo$Cilla

GUESS <- ctmm.guess(Cilla,interactive=FALSE)
FIT <- ctmm.fit(Cilla,GUESS)

# stationary Gaussian estimate
speed(FIT)

# conditional estimate
\donttest{speed(FIT,Cilla)}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

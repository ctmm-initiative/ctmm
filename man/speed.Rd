\name{speed}
\alias{speed}
\alias{speed.ctmm}
\alias{speed.telemetry}
\alias{speeds}
\alias{speeds.ctmm}
\alias{speeds.telemetry}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate the average speed of a tracked animal}
\description{
Given a \code{ctmm} movement model and \code{telemetry} data, \code{speed} simulates multiple realizations of the individual's trajectory to estimate the time-averaged speed, which is proportional to distance traveled, while \code{speeds} estimates instantaneous speeds at a specified array of times \code{t}. Both tortuosity (non straight-line motion between the data) and telemetry error can be accounted for. Given only a \code{ctmm} movement model and no data, \code{speed} calculates the average speed of the Gaussian movement process.
}
\usage{
speed(object,...)

\S3method{speed}{ctmm}(object,data=NULL,level=0.95,robust=FALSE,units=TRUE,prior=TRUE,fast=TRUE,
      cor.min=0.5,dt.max=NULL,error=0.01,cores=1,...)

\S3method{speed}{telemetry}(object,CTMM,level=0.95,robust=FALSE,units=TRUE,prior=TRUE,fast=TRUE,cor.min=0.5,
      dt.max=NULL,error=0.01,cores=1,...)

speeds(object,...)

\S3method{speeds}{ctmm}(object,data=NULL,t=NULL,cycle=Inf,level=0.95,robust=FALSE,prior=FALSE,fast=TRUE,
       error=0.01,cores=1,...)

\S3method{speeds}{telemetry}(object,CTMM,t=NULL,cycle=Inf,level=0.95,robust=FALSE,prior=FALSE,fast=TRUE,
       error=0.01,cores=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{ctmm} movement-model or \code{telemetry} object, which requires an additional \code{CTMM} argument.}
  \item{data}{Optional \code{telemetry} object on which the simulations will be conditioned.}
  \item{CTMM}{Movement model object.}
  \item{t}{Array of times to estimate instantaneous speeds at.}
  \item{cycle}{Average over time \code{t} indices modulo \code{cycle}. E.g., for \code{t} sequenced by hours, \code{cycle=24} gives daily the cycle of speeds. (Not yet supported.)}
  \item{level}{Confidence level to report on the estimated average speed.}
  \item{robust}{Use robust statistics for the ensemble average and its confidence intervals (see Details).}
  \item{units}{Convert result to natural units.}
  \item{prior}{Account for model uncertainty.}
  \item{fast}{Whether or not to invoke the central-limit theorem (see \code{\link{emulate}}) when \code{prior=TRUE}.}
  \item{cor.min}{Velocity correlation threshold for skipping gaps.}
  \item{dt.max}{Absolute gap sizes to skip (in seconds), alternative to \code{cor.min}.}
  \item{error}{Target (relative) standard error.}
  \item{cores}{Number of simulations to run in parallel. \code{cores=0} will use all cores, while \code{cores<0} will reserve \code{abs(cores)}.}
  \item{...}{Arguments passed to \code{\link{emulate}}.}
}
\details{
The \code{cor.min} or \code{dt.max} arguments are used to constrain the estimate to be derived from simulations near the data, and therefore ensure that the estimate is more reflective of the data than the model.

If data quality is poor and velocity can barely be resolved, then the sampling distribution may occassionally include impersistent motion and the ensemble mean will be infinite. In these cases \code{robust=TRUE} can be used to report the ensemble median rather than the ensemble mean. The time average of \code{speed}, in either case, is still the mean and the resulting quantity is still proportional to distance traveled. Furthermore, note that medians should be compared to medians and means to means, so the \code{robust} option should be the same for all compared individuals.
}

\value{ Returns the estimated mean speed of the sampled trajectory with CIs by default. If \code{level=NULL}, then the ensemble of mean speeds is returned instead. }

%\references{ }

\author{ C. H. Fleming. }

\note{ The average speeds estimated here are mean speeds. The speeds reported by \code{\link{summary.ctmm}} are root-mean-square (RMS) speeds. These quantities are sometimes proportional, but not equivalent. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{emulate}}, \code{\link{simulate}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)
Gabs <- buffalo$Gabs

GUESS <- ctmm.guess(Gabs,interactive=FALSE)
FIT <- ctmm.fit(Gabs,GUESS)

# stationary Gaussian estimate
speed(FIT)

# conditional estimate
speed(FIT,Gabs)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

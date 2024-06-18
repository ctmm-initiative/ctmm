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
Given a \code{ctmm} movement model and \code{telemetry} data, \code{speed} simulates multiple realizations of the individual's trajectory to estimate the time-averaged speed, which is proportional to distance traveled, while \code{speeds} estimates instantaneous speeds at a specified array of times \code{t}. Both tortuosity (non straight-line motion between the data) and telemetry error can be accounted for. Given only a \code{ctmm} movement model and no data, \code{speed} calculates the mean speed of the Gaussian movement process. All methods are described in Noonan & Fleming et al (2019).
}
\usage{
speed(object,...)

\S3method{speed}{ctmm}(object,data=NULL,t=NULL,level=0.95,robust=FALSE,units=TRUE,prior=TRUE,fast=TRUE,
      cor.min=0.5,dt.max=NULL,error=0.01,cores=1,trace=TRUE,...)

\S3method{speed}{telemetry}(object,CTMM,t=NULL,level=0.95,robust=FALSE,units=TRUE,prior=TRUE,fast=TRUE,
      cor.min=0.5,dt.max=NULL,error=0.01,cores=1,trace=TRUE,...)

speeds(object,...)

\S3method{speeds}{ctmm}(object,data=NULL,t=NULL,cycle=Inf,level=0.95,robust=FALSE,prior=FALSE,fast=TRUE,
       error=0.01,cores=1,trace=TRUE,...)

\S3method{speeds}{telemetry}(object,CTMM,t=NULL,cycle=Inf,level=0.95,robust=FALSE,prior=FALSE,fast=TRUE,
       error=0.01,cores=1,trace=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{ctmm} movement-model or \code{telemetry} object, which requires an additional \code{CTMM} argument.}
  \item{data}{Optional \code{telemetry} object on which the simulations will be conditioned.}
  \item{CTMM}{Movement model object.}
  \item{t}{Array of times to estimate instantaneous speeds at, or range of times to estimate mean speed over.}
  \item{cycle}{Average over time \code{t} indices modulo \code{cycle}. E.g., for \code{t} sequenced by hours, \code{cycle=24} gives daily the cycle of speeds. (Not yet supported.)}
  \item{level}{Confidence level to report on the estimated average speed.}
  \item{robust}{Use robust statistics for the ensemble average and its confidence intervals (see Details).}
  \item{units}{Convert result to natural units.}
  \item{prior}{Account for model parameter uncertainty.}
  \item{fast}{Whether or not to invoke the central-limit theorem when propagating parameter uncertainty (see \code{\link{emulate}}).}
  \item{cor.min}{Velocity correlation threshold for skipping gaps.}
  \item{dt.max}{Absolute gap sizes to skip (in seconds), alternative to \code{cor.min}.}
  \item{error}{Target (relative) standard error.}
  \item{cores}{Number of simulations to run in parallel. \code{cores=0} will use all cores, while \code{cores<0} will reserve \code{abs(cores)}.}
  \item{trace}{Display a progress bar.}
  \item{...}{Arguments passed to \code{\link{emulate}}.}
}
\details{
The \code{cor.min} or \code{dt.max} arguments are used to constrain the estimate to be derived from simulations near the data, and therefore ensure that the estimate is more reflective of the data than the model.

If data quality is poor and velocity can barely be resolved, then the sampling distribution may occassionally include impersistent motion and its mean will be infinite. In these cases \code{robust=TRUE} can be used to report the sampling distribution's median rather than its mean. The time average of \code{speed}, in either case, is still the mean average of times and the resulting quantity is still proportional to distance traveled. Furthermore, note that medians should be compared to medians and means to means, so the \code{robust} option should be the same for all compared individuals.
}

\value{ Returns the estimated mean speed of the sampled trajectory with CIs by default. The \code{DOF} slot corresponds to a scaled-\eqn{\chi} sampling distribution. If \code{level=NULL}, then the ensemble of mean speeds is returned instead. }

\references{
M. J. Noonan, C. H. Fleming, T. S. Akre, J. Drescher-Lehman, E. Gurarie, A.-L. Harrison, R. Kays, Justin Calabrese,
``Scale-insensitive estimation of speed and distance traveled from animal tracking data'',
\href{https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-019-0177-1}{Movement Ecology, 7:35 (2019)}.
}

\author{ C. H. Fleming. }

\note{ The mean speed estimated by \code{speed} is applicable only during the sampling periods. If an individual is diurnal/nocturnal and only tracked during the day/night, then the output of \code{speed} will only be the mean speed during the day/night. For instance, if an individual is tracked the 12 hours per day during which it is active, and \code{speed} reports a mean speed of 10 kilometers per day during those periods, then the average distance traveled per day is only 5 kilometers (from 10 kilometers \code{/} day \code{*} 12 hours). An average of 10 kilometers would only result if the individual were similarly active for 24 hours a day.

The average speeds estimated here are mean speeds. The speeds reported by \code{\link{summary.ctmm}} are root-mean-square (RMS) speeds. These quantities are sometimes proportional, but not equivalent. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{emulate}}, \code{\link{simulate}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)
DATA <- buffalo$Gabs

GUESS <- ctmm.guess(DATA,interactive=FALSE)
# in general, you should use ctmm.select instead
FIT <- ctmm.fit(DATA,GUESS)

# stationary Gaussian estimate
speed(FIT)

# conditional estimate
# you will likely want trace=TRUE
speed(FIT,DATA,trace=FALSE)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

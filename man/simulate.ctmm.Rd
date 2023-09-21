\name{simulate.ctmm}
\alias{predict}
\alias{predict.ctmm}
\alias{predict.telemetry}
\alias{simulate}
\alias{simulate.ctmm}
\alias{simulate.telemetry}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict or simulate from a continuous-time movement model}
\description{
Given a \code{ctmm} movement model (and optional \code{telemetry} data to condition upon) these functions predict or simulate animal locations over a prescribed set of times.
}
\usage{
predict(object,...)

\S3method{predict}{ctmm}(object,data=NULL,VMM=NULL,t=NULL,dt=NULL,res=1,complete=FALSE,...)

\S3method{predict}{telemetry}(object,CTMM=NULL,VMM=NULL,t=NULL,dt=NULL,res=1,complete=FALSE,...)

simulate(object,nsim=1,seed=NULL,...)

\S3method{simulate}{ctmm}(object,nsim=1,seed=NULL,data=NULL,VMM=NULL,t=NULL,dt=NULL,res=1,complete=FALSE,
         precompute=FALSE,...)

\S3method{simulate}{telemetry}(object,nsim=1,seed=NULL,CTMM=NULL,VMM=NULL,t=NULL,dt=NULL,res=1,complete=FALSE,
         precompute=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{ctmm} movement-model or \code{telemetry} object, which requires an additional \code{CTMM} argument.}
  \item{data}{Optional \code{telemetry} object on which the prediction or simulation will be conditioned.}
  \item{CTMM}{A \code{ctmm} movement model in the same format as the output of \code{\link{ctmm.fit}} or \code{\link{variogram.fit}}.}
  \item{VMM}{An optional vertical \code{ctmm} movement model for 3D predictions and simulations.}
  \item{t}{Optional array of numeric time values over which the process will be predicted or simulated.}
  \item{dt}{Timestep to space the prediction or simulation over if \code{data} is specified.}
  \item{res}{Average number of locations to predict or simulate per \code{data} time.}
  \item{complete}{Additionally calculate timestamps and geographic coordinates.}
  \item{nsim}{Generates a list of \code{nsim} simulations.}
  \item{seed}{Optional random seed to fix.}
  \item{precompute}{Precalculate matrices of the Kalman filter (see details).}
  \item{...}{Unused options.}
}
\details{
The prediction or simulation necessarily requires a \code{ctmm} model object.
If a \code{telemetry} \code{data} object is supplied, the output will be conditional on the \code{data} (i.e., simulations that run through the data).
If no \code{data} is provided then the output will be purely Gaussian, and times \code{t} must be provided.
Details of the movement model parameters can be found in \code{\link{ctmm.fit}}.

The \code{t} argument fixes the output times to a specific array of times.
The \code{dt} and \code{res} arguments are relative to the sampling schedule present in the optional \code{telemetry} object.
The same span of time will be used, while \code{dt} will fix the sampling rate absolutely and \code{res} will fix the sampling rate relative to that of the data.

The \code{precompute} option can speed up calculations of multiple simulations of the same model, data, and \emph{irregular} sampling schedule.
First run \code{simulate} with \code{precompute=TRUE} to calculate and store all of the necessary matrices of the Kalman filter.
A simulated \code{telemetry} object will be produced, as usual, and the precomputed objects are stored in the environment.
Subsequent simulations with \code{precompute=-1} will then apply these precomputed matrices for a computational cost savings.
If the sampling schedule is irregular, then this can result in faster simulations.
}

\value{ A simulated animal-tracking \code{telemetry} object with components \code{t}, \code{x}, and \code{y}, or a predicted \code{telemetry} object that also includes \code{x}-\code{y} covariances for the location point estimates \code{x} and \code{y}. }

\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, W. F. Fagan,
``From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales'',
The American Naturalist, 183:5, E154-E167 (2014) \doi{10.1086/675504}.

C. H. Fleming, D. Sheldon, E. Gurarie, W. F. Fagan, S. LaPoint, J. M. Calabrese,
``Kálmán filters for continuous-time movement models'',
Ecological Informatics, 40, 8-21 (2017) \doi{10.1016/j.ecoinf.2017.04.008}.
}


\author{ C. H. Fleming. }

\note{ Predictions are autocorrelated and should not be treated as data. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}} }

\examples{\donttest{
#Load package
library(ctmm)

#prepare simulation parameters
t <- 1:1000
MODEL <- ctmm(tau=c(100,10),sigma=10,mu=c(0,0))

#simulate data
SIM <- simulate(MODEL,t=t)

#plot data with Gaussian model
plot(SIM,CTMM=MODEL)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{simulate.ctmm}
\alias{predict}
\alias{predict.ctmm}
\alias{predict.telemetry}
\alias{simulate}
\alias{simulate.ctmm}
\alias{simulate.telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict or simulate from a continuous-time movement model}
\description{
Given a \code{ctmm} movement model (and optional \code{telemetry} data to condition upon) these functions predict or simulate animal locations over a prescribed set of times.
}
\usage{
predict(object,...)

\S3method{predict}{ctmm}(object,data=NULL,t=NULL,dt=NULL,res=1,...)

\S3method{predict}{telemetry}(object,CTMM=NULL,t=NULL,dt=NULL,res=1,...)

simulate(object,nsim=1,seed=NULL,...)

\S3method{simulate}{ctmm}(object,nsim=1,seed=NULL,data=NULL,t=NULL,dt=NULL,res=1,...)

\S3method{simulate}{telemetry}(object,nsim=1,seed=NULL,CTMM=NULL,t=NULL,dt=NULL,res=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{ctmm} movement-model or \code{telemetry} object, which requires an additional \code{CTMM} argument.}
  \item{data}{Optional \code{telemetry} object on which the prediction or simulation will be conditioned.}
  \item{t}{Array of numeric time values over which the process will be predicted or simulated.}
  \item{dt}{Timestep to space the prediction or simulation over if \code{data} is specified.}
  \item{res}{Average number of locations to predict or simulate per \code{data} time.}
  \item{CTMM}{A \code{ctmm} movement-model in the same format as the output of \code{\link{ctmm.fit}} or \code{\link{variogram.fit}}.}
  \item{nsim}{Not yet supported.}
  \item{seed}{Optional random seed to fix.}
  \item{...}{Unused options.}
}
\details{
The prediction or simulation necessarily requires a \code{ctmm} model object and optionally can be conditioned off of \code{telemetry} data, if specified.
If no data is provided, the simulation will be purely Gaussian.
Details of the movement model parameters can be found in \code{\link{ctmm.fit}}.

The \code{t} argument fixes the output times to a specific array of times.
The \code{dt} and \code{res} arguments are relative to the sampling schedule present in the optional \code{telemetry} object.
The same span of time will be used, while \code{dt} will fix the sampling rate absolutely and \code{res} will fix the sampling rate relative to that of the data.
}

\value{ A simulated animal-tracking \code{telemetry} object with components \code{t}, \code{x}, and \code{y}, or a predicted \code{telemetry} object that also includes \code{x}-\code{y} covariances for the location point estimates \code{x} and \code{y}. }

\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, and W. F. Fagan. (2014). From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales. \href{http://www.jstor.org/discover/10.1086/675504}{The American Naturalist, 183(5), E154-E167}.

C. H. Fleming and Y. Subasi and J. M. Calabrese. (2015). A maximum-entropy description of animal movement. \href{http://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.032107}{Physical Review E, 91, 032107}.

C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese. (2016). Estimating where and how animals travel: An optimal framework for path reconstruction from autocorrelated tracking data. \href{http://onlinelibrary.wiley.com/doi/10.1890/15-1607/full}{Ecology, DOI:10.1890/15-1607}.
}


\author{ C. H. Fleming. }

\note{ Predictions are autocorrelated and should not be treated as data. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}} }

\examples{
#Load package and data
library(ctmm)

#prepare simulation parameters
t <- 1:1000
MODEL <- ctmm(tau=c(100,10),sigma=10,mu=c(0,0))

#simulate data
SIM <- simulate(MODEL,t=t)

#plot data with Gaussian model
plot(SIM,CTMM=MODEL)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

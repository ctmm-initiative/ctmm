\name{simulate.ctmm}
\alias{simulate}
\alias{simulate.ctmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a continuous-time movement model
}
\description{
This function generates an \eqn{O(n)} hidden-Markov simulation of continuous-time movement models described in Fleming et al (2014) and Fleming et al (2015).
}
\usage{
simulate(object,nsim=1,seed=NULL,...)

\S3method{simulate}{ctmm}(object,nsim=1,seed=NULL,data=NULL,t=NULL,dt=NULL,res=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{ctmm} movement-model in the same format as the output of \code{\link{ctmm.fit}} or \code{\link{variogram.fit}}.}
  \item{nsim}{Not yet supported.}
  \item{seed}{Optional random seed to fix.}
  \item{data}{Optional \code{telemetry} object on which the simulation will be conditioned.}
  \item{t}{Array of numeric time values over which the process will be simulated.}
  \item{dt}{Timestep to run simulation over if \code{data} is specified.}
  \item{res}{Average number of locations to simulate per \code{data} time.}
  \item{...}{Unused options.}
}
\details{
The simulation can be conditioned off of \code{data}, if specified. Otherwise, the simulation will be Gaussian. 
Details of the movement model parameters can be found in \code{\link{ctmm.fit}}.
}

\value{ A simulated animal-tracking \code{telemetry} object with components \code{t}, \code{x}, and \code{y}. }

\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, and W. F. Fagan. (2014). From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales. \href{http://www.jstor.org/discover/10.1086/675504}{The American Naturalist, 183(5), E154-E167.}

C. H. Fleming and Y. Subasi and J. M. Calabrese. (2015). A maximum-entropy description of animal movement. \href{http://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.032107}{Physical Review E, 91, 032107.}
}


\author{ C. H. Fleming. }
% \note{ }

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

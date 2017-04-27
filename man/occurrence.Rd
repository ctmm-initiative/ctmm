\name{occurrence}
\alias{occurrence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate a Kriged occurrence disribution estimate}
\description{
This function calculates an occurrence distribution from \code{telemetry} data and a continuous-time movement model.
}
\usage{
occurrence(data,CTMM,H=0,res.time=10,res.space=10,grid=NULL,
cor.min=0.5,dt.max=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ 2D timeseries telemetry data represented as a \code{telemetry} object. }
  \item{CTMM}{A \code{ctmm} movement model from the output of \code{ctmm.fit}.}
  \item{H}{Optional additional banwidth matrix for future use.}
  \item{res.time}{Number of temporal grid points per median timestep.}
  \item{res.space}{Number of grid points along each axis, relative to the average diffusion (per median timestep) from a stationary point.}
  \item{grid}{Optional grid specification with columns labeled \code{x} and \code{y}. Not yet supported.}
  \item{cor.min}{Minimum gap correlation for Kriging interpolation (alternative to \code{dt.max}).}
  \item{dt.max}{Maximum gap size for Kriging interpolation (alternative to \code{cor.min}).}
}

%\details{}

\value{Returns a \code{UD} object containing the sampled grid line locations \code{x} and \code{y}, the probability density and cumulative distribution functions evaluated on the sampled grid locations \code{PDF} & \code{CDF}, the optional bandwidth matrix \code{H}, and the area of each grid cell \code{dA}.
}

\references{
C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese. (2016). Estimating where and how animals travel: An optimal framework for path reconstruction from autocorrelated tracking data. \href{http://onlinelibrary.wiley.com/doi/10.1890/15-1607/full}{Ecology, DOI:10.1890/15-1607}.
}

\author{ C. H. Fleming. }

\note{Large gaps have a tendency to blow up the estimate, and can be avoided with the \code{cor.min} or \code{dt.max} arguments.

In the case of coarse grids, the value of \code{PDF} in a grid cell actually corresponds to the average probability density over the entire rectangular cell. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{raster,UD-method}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)
Cilla <- buffalo$Cilla

GUESS <- ctmm.guess(Cilla,interactive=FALSE)
FIT <- ctmm.fit(Cilla,GUESS)

# Compute occurence distribution
UD <- occurrence(Cilla,FIT)

# Plot occurrence UD
plot(UD,col.level=NA)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

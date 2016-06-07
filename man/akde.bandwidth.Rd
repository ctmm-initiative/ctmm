\name{akde.bandwidth}
\alias{akde.bandwidth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the optimal bandwidth matrix of movement data
}
\description{
This function calculates the optimal bandwidth matrix (kernel covariance) for a two-dimensional animal tracking dataset, given an autocorrelated movement model (Fleming et al, 2015). This optimal bandwidth can fully take into account all autocorrelation in the data, assuming it is captured by the movement model.
}
\usage{
akde.bandwidth(data,CTMM,fast=NULL,dt=NULL,verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ 2D timeseries telemetry data represented as a \code{telemetry} object. }
  \item{CTMM}{A \code{ctmm} movement model from the output of \code{ctmm.fit}.}
  \item{fast}{ Use the direct algorithm if \code{FALSE} or the FFT algorithm if \code{TRUE}. The default behavior switches from slow to fast at \code{length(data$t)==100}.}
  \item{dt}{ Optional lag bin width for the FFT algorithm.}
  \item{verbose}{ Optionally return the effective sample size \code{DOF.H} along with the bandwidth matrix \code{H}.}
}
\details{
One of two algorithms is used to count the time lags. The direct method is always exact but comes at an \eqn{O(n^2)} computational cost. The fast \eqn{O(n \log n)} algorithm is based upon the FFT method of Marcotte (1996). This FFT method is exact for evenly scheduled data, and can account for missing observations.

If no \code{dt} is specified, the median sampling interval is used if \code{fast=TRUE}. This is typically a good assumption for most data, even when there are gaps. By decreasing \code{dt}, the FFT method can be made arbitrarily precise with no ill effects.
}

\value{Returns a bandwidth \code{matrix} object, which is to be the optimal covariance matrix of the individual kernels of the kernel density estimate.}

\references{
C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese (2015).
Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
\href{http://www.esajournals.org/doi/abs/10.1890/14-2010.1}{Ecology, 96(5), 1182-1188}.

D. Marcotte. (1996). Fast variogram computation with FFT. \emph{Computers and Geosciences} 22(10), 1175-1186.
}

\author{ C. H. Fleming. }

\note{ To obtain a bandwidth scalar representing the variance of each kernel, a \code{ctmm} object with \code{isotropic=TRUE} is required.  In this case, \code{akde.bandwidth} will return bandwidth matrix with identical variances along its diagonal. Note that this will provide an inaccurate estimate for very eccentric distributions}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{ctmm.fit}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)
cilla <- buffalo[[1]]

# Fit a continuous-velocity model with tau ~ c(10 days, 1 hour)
# see help(variogram.fit)
GUESS <- ctmm(tau=c(10*24*60^2,60^2))
FIT <- ctmm.fit(cilla,GUESS)

# Optimize bandwidth matrix
H <- akde.bandwidth(cilla,FIT)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

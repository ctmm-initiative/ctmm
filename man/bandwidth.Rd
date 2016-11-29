\name{bandwidth}
\alias{bandwidth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the optimal bandwidth matrix of movement data
}
\description{
This function calculates the optimal bandwidth matrix (kernel covariance) for a two-dimensional animal tracking dataset, given an autocorrelated movement model (Fleming et al, 2015). This optimal bandwidth can fully take into account all autocorrelation in the data, assuming it is captured by the movement model.
}
\usage{
bandwidth(data,CTMM,VMM=NULL,weights=FALSE,fast=TRUE,
   dt=NULL,precision=1/2,PC="Markov",verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ 2D timeseries telemetry data represented as a \code{telemetry} object. }
  \item{CTMM}{A \code{ctmm} movement model as from the output of \code{ctmm.fit}.}
  \item{VMM}{An optional vertical \code{ctmm} object for 3D bandwidth calculation.}
  \item{weights}{By default, the weights are taken to be uniform. \code{weights=TRUE} will optimize the weights, while a \code{numeric} array will fix the weights.}
  \item{fast}{Use FFT algorithms for weight optimization.}
  \item{dt}{ Optional lag bin width for the FFT algorithm.}
  \item{precision}{Fraction of maximum possible digits of precision to target in weight optimization. \code{precision=1/2} results in about 7 decimal digits of precision if the preconditioner is stable.}
  \item{PC}{Preconditioner to use: can be "Markov", "circulant", "IID", or "direct".}
  \item{verbose}{ Optionally return the optimal \code{weights}, effective sample size \code{DOF.H}, and other information along with the bandwidth matrix \code{H}.}
}
\details{
The \code{weights} argument can be used to correct fixed sampling bias with a prescribed \code{numeric} array or to correct temporal sampling bias caused by autocorrelation with \code{weights=TRUE}.
\code{weights=TRUE} will optimize \code{n=length(data$t)} weights via constrained & preconditioned conjugate gradient algorithms.
These algorithms have a few options that should be considered if the data are very irregular.

\code{FFT=TRUE} grids the data with grid width \code{dt} and applies FFT algorithms, for a computational cost as low as \eqn{O(n \log n)} with only \eqn{O(n)} function evaluations.
If no \code{dt} is specified, the minimum sampling interval is used, but \strong{if the data are very irregular, then \code{dt} may need to be several times smaller} to avoid slow down.
A \code{dt} Insufficient to resolve the details of the sampling schedule will cause an excessive, \eqn{O(n)} number of feasibility assessments.
%With the sub-optimal \code{PC="IID"} preconditioner, \code{dt} can be relaxed to the median value, but the computational cost is then increased to \eqn{O(n^2 \log n)}.

\code{FFT=FALSE} uses exact times and has a computational cost as low as \eqn{O(n^2)}, including \eqn{O(n^2)} function evaluations. With \code{PC="direct"} this method will produce a result that is exact to within machine precision, but with a computational cost of \eqn{O(n^3)}.

}

\value{Returns a bandwidth \code{matrix} object, which is to be the optimal covariance matrix of the individual kernels of the kernel density estimate.}

\references{
T. F. Chan (1988).
An Optimal Circulant Preconditioner for Toeplitz Systems.
SIAM Journal on Scientific and Statistical Computing, 9(4), 766-771.

C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese (2015).
Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
\href{http://www.esajournals.org/doi/abs/10.1890/14-2010.1}{Ecology, 96(5), 1182-1188}.

D. Marcotte. (1996). Fast variogram computation with FFT. \emph{Computers and Geosciences} 22(10), 1175-1186.
}

\author{ C. H. Fleming. }

\note{ To obtain a bandwidth scalar representing the variance of each kernel, a \code{ctmm} object with \code{isotropic=TRUE} is required.  In this case, \code{bandwidth} will return bandwidth matrix with identical variances along its diagonal. Note that this will provide an inaccurate estimate for very eccentric distributions.}

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
H <- bandwidth(cilla,FIT,verbose=TRUE)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

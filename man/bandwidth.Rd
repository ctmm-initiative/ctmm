\name{bandwidth}
\alias{bandwidth}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the optimal bandwidth matrix of movement data}
\description{This function calculates the optimal bandwidth matrix (kernel covariance) for a two-dimensional animal tracking dataset, given an autocorrelated movement model (Fleming et al, 2015). This optimal bandwidth can fully take into account all autocorrelation in the data, assuming it is captured by the movement model.}
\usage{
bandwidth(data,CTMM,VMM=NULL,weights=FALSE,fast=NULL,dt=NULL,PC="Markov",error=0.01,
          precision=1/2,verbose=FALSE,trace=FALSE,dt.plot=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ 2D timeseries telemetry data represented as a \code{telemetry} object. }
  \item{CTMM}{A \code{ctmm} movement model as from the output of \code{ctmm.fit}.}
  \item{VMM}{An optional vertical \code{ctmm} object for 3D bandwidth calculation.}
%  \item{UD}{A list of \code{UD} objects for calculating population weights and bandwidth [IN DEVELOPMENT].}
%  \item{h}{Fixed bandwidth value.}
  \item{weights}{By default, the weights are taken to be uniform, whereas \code{weights=TRUE} will optimize the weights.}
  \item{fast}{Use FFT algorithms for weight optimization. \code{fast=NULL} will attempt to intelligently decide between the fast and exact algorithms based on computational complexity.}
  \item{dt}{Optional lag bin width for the FFT algorithm.}
  \item{PC}{Preconditioner to use: can be "Markov", "circulant", "IID", or "direct".}
  \item{error}{Maximum grid error for FFT algorithm, if \code{dt} is not specified.}
  \item{precision}{Fraction of maximum possible digits of precision to target in weight optimization. \code{precision=1/2} results in about 7 decimal digits of precision if the preconditioner is stable.}
  \item{verbose}{ Optionally return the optimal \code{weights}, effective sample size \code{DOF.H}, and other information along with the bandwidth matrix \code{H}.}
  \item{trace}{Produce tracing information on the progress of weight optimization.}
  \item{dt.plot}{Execute a diagnostic \code{\link{dt.plot}} with a red line at \code{dt}, if \code{weights=TRUE}.}
  \item{...}{Arguments passed to \code{\link{mean.ctmm}}.}
}
\details{
The \code{weights=TRUE} argument can be used to correct temporal sampling bias caused by autocorrelation.
\code{weights=TRUE} will optimize \code{n=length(data$t)} weights via constrained & preconditioned conjugate gradient algorithms.
These algorithms have a few options that should be considered if the data are very irregular.

\code{fast=TRUE} is an approximation that discretizes the data with timestep \code{dt} and applies FFT algorithms, for a computational cost as low as \eqn{O(n \log n)} with only \eqn{O(n)} function evaluations.
If no \code{dt} is specified, then a choice of \code{dt} will be automated with a message.
\strong{If the data contain some very tiny time intervals}, say 1 second among hourly sampled data, then the default \code{dt} setting can create an excessively high-resolution discretization of time, which will cause slowdown. In this case \code{CTMM} should contain a location-error model and \code{dt} should be increased to a larger fraction of the most-frequent sampling intervals.
\strong{If the data are irregular (permitting gaps), then \code{dt} may need to be several times smaller} than the median to avoid slow down.
In this case, try setting \code{trace=TRUE} and decreasing \code{dt} below the median until the interations speed up and the number of feasibility assessments becomes less than \eqn{O(n)}.
%With the sub-optimal \code{PC="IID"} preconditioner, \code{dt} can be relaxed to the median value, but the computational cost is then increased to \eqn{O(n^2 \log n)}.

\code{fast=FALSE} uses exact time spacing and has a computational cost as low as \eqn{O(n^2)}, including \eqn{O(n^2)} function evaluations. With \code{PC="direct"} this method will produce a result that is exact to within machine precision, but with a computational cost of \eqn{O(n^3)}. \strong{\code{fast=FALSE,PC='direct'} is often the fastest method with small datasets}, where \eqn{n \le O}(1,000), but scales terribly with larger datasets.
}

\value{Returns a bandwidth \code{matrix} object, which is to be the optimal covariance matrix of the individual kernels of the kernel density estimate.}

\references{
T. F. Chan,
``An Optimal Circulant Preconditioner for Toeplitz Systems'',
SIAM Journal on Scientific and Statistical Computing, 9:4, 766-771 (1988) \doi{10.1137/0909051}.

D. Marcotte, ``Fast variogram computation with FFT'', Computers and Geosciences 22:10, 1175-1186 (1996) \doi{10.1016/S0098-3004(96)00026-X}.

C. H. Fleming, W. F. Fagan, T. Mueller, K. A. Olson, P. Leimgruber, J. M. Calabrese,
``Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator'',
Ecology, 96:5, 1182-1188 (2015) \doi{10.1890/14-2010.1}.

C. H. Fleming, D. Sheldon, W. F. Fagan, P. Leimgruber, T. Mueller, D. Nandintsetseg, M. J. Noonan, K. A. Olson, E. Setyawan, A. Sianipar, J. M. Calabrese,
``Correcting for missing and irregular data in home-range estimation'',
Ecological Applications, 28:4, 1003-1010 (2018) \doi{10.1002/eap.1704}.
}

\author{ C. H. Fleming. }

\note{ To obtain a bandwidth scalar representing the variance of each kernel, a \code{ctmm} object with \code{isotropic=TRUE} is required.  In this case, \code{bandwidth} will return bandwidth matrix with identical variances along its diagonal. Note that forcing \code{isotropic=TRUE} will provide an inaccurate estimate for very eccentric distributions.

In v1.0.1 the default \code{fast}, \code{dt}, \code{PC} arguments depend on the sample size, with \code{fast=FALSE}, \code{PC="Direct"} for small sample sizes, \code{fast=FALSE}, \code{PC="Markov"} for moderate sample sizes, and \code{fast=TRUE}, \code{PC="Markov"} for large sample sizes, where \code{dt} is taken to be the integer fraction of the median sampling interval closest to the minimum sampling interval.

In v0.6.2 the default \code{dt} was increased form the minimum time difference to a small quantile no less than \code{error} times the median.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{ctmm.fit}} }

%\examples{\donttest{
%# Load package and data
%library(ctmm)
%data(buffalo)
%cilla <- buffalo[[1]]
%
%# Fit a continuous-velocity model with tau ~ c(10 days, 1 hour)
%# see help(variogram.fit)
%GUESS <- ctmm(tau=c(10*24*60^2,60^2))
%FIT <- ctmm.fit(cilla,GUESS)
%
%# Optimize bandwidth matrix
%H <- bandwidth(cilla,FIT,verbose=TRUE)
%}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

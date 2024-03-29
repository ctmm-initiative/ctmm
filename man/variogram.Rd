\name{variogram}
\alias{variogram}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate an empirical variogram from movement data}
\description{
This function calculates the empirical variogram of multi-dimensional tracking data for visualizing stationary (time-averaged) autocorrelation structure.
One of two algorithms is used. The slow \eqn{O(n^2)} algorithm is based upon Fleming & Calabrese et al (2014), but with interval-weights instead of lag-weights and an iterative algorithm to adjust for calibrated errors.
Additional modifications have also been included to accommodate drift in the sampling rate.
The fast \eqn{O(n \log n)} algorithm is based upon the FFT method of Marcotte (1996), with some tweaks to better handle irregularly sampled data.
Both methods reduce to the unbiased \dQuote{method of moments} estimator in the case of evenly \emph{scheduled} data, even with missing observations, but they produce slightly different outputs for irregularly sampled data.
}
\usage{
variogram(data,dt=NULL,fast=TRUE,res=1,CI="Markov",error=FALSE,axes=c("x","y"),
          precision=1/8,trace=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ \code{telemetry} data object of the 2D timeseries data. }
  \item{dt}{ Lag bin width. An ordered array will yield a progressive coarsening of the lags. Defaults to the median sampling interval.}
  \item{fast}{ Use the interval-weighted algorithm if \code{FALSE} or the FFT algorithm if \code{TRUE}. The slow algorithm outputs a progress bar.}
  \item{res}{Increase the discretization resolution for irregularly sampled data with \code{res>1}. Decreases bias at the cost of smoothness.}
  \item{CI}{Argument for confidence-interval estimation. Can be \code{"IID"} to consider all unique lags as independent, \code{"Markov"} to consider only non-overlapping lags as independent, or \code{"Gauss"} for an exact calculation (see Details below).}
  \item{error}{Adjust for the effect of calibrated errors.}
  \item{axes}{Array of axes to calculate an average (isotropic) variogram for.}
  %\item{max.lag}{Maximum lag of the variogram to calculate when \code{fast=FALSE}, so to minimize computational cost.}
  \item{precision}{Fraction of machine precision to target when adjusting for telemetry error (\code{fast=FALSE} with calibrated errors). \code{precision=1/8} returns about 2 decimal digits of precision.}
  \item{trace}{Display a progress bar if \code{fast=FALSE}.}
}
\details{
If no \code{dt} is specified, the median sampling interval is used. This is typically a good assumption for most data, even when there are gaps. A \code{dt} coarser than the sampling interval may bias the variogram (particuarly if \code{fast=TRUE}) and so this should be reserved for poor data quality.

For irregularly sampled data, it may be useful to provide an array of time-lag bin widths to progressively coarsen the variogram. I.e., if you made the very bad choice of changing your sampling interval on the fly from \code{dt1} to \code{dt2}, where \code{dt1} \eqn{<} \code{dt2}, the an appropriate choice would be \code{dt=c(dt1,dt2)}. On the other hand, if your sampling is itself a noisy process, then you might want to introduce larger and larger \code{dt} components as the visual appearance of the variogram breaks down with increasing lags.
Alternatively, you might try the \code{fast=FALSE} option or aggregating multiple individuals with \code{\link{mean.variogram}}.

With irregularly sampled data, different size lags must be aggregated together, and with current fast methods there is a tradeoff between bias and smoothness. The default settings produce a relatively smooth estimate, while increasing \code{res} (or setting \code{fast=FALSE}) will produce a less biased estimate, which is very useful for \code{\link{correlogram}}.

In conventional variogram regression treatments, all lags are considered as independent (\code{CI="IID"}) for the purposes of confidence-interval estimation, even if they overlap in time. However, in high resolution datasets this will produce vastly underestimated confidence intervals. Therefore, the default \code{CI="Markov"} behavior is to consider only the maximum number of non-overlapping lags in calculating confidence intervals, though this is a crude approximation and is overly conservative at large lags. \code{CI="Gauss"} implements exact confidence intervals under the assumption of a stationary Gaussian process, but this algorithm is \eqn{O(n^2 \log n)} even when \code{fast=TRUE}.

If \code{fast=FALSE} and the tracking data are calibrated (see \code{\link{uere}}), then with \code{error=TRUE} the variogram of the movement process (sans the telemetry-error process) is estimated using an iterative maximum-likelihood esitmator that downweights more erroneous location estimates (Fleming et al, 2020). The variogram is targeted to have \code{precision} fraction of machine precision. If the data are very irregular and location errors are very homoskedastic, then this algorithm can be slow to converge at time lags where there are few data pairs.
If \code{fast=TRUE} and \code{error=TRUE}, then the estimated contribution to the variogram from location error is subtracted on a per lag basis, which is less ideal for heteroskedastic errors.
}

\value{
Returns a variogram object (class variogram) which is a dataframe containing the time-lag, \code{lag}, the semi-variance estimate at that lag, \code{SVF}, and the approximate number of degrees of freedom associated with that semi-variance, \code{DOF}, with which its confidence intervals can be estimated.
}
\references{
D. Marcotte, ``Fast variogram computation with FFT'', Computers and Geosciences 22:10, 1175-1186 (1996) \doi{10.1016/S0098-3004(96)00026-X}.

C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, W. F. Fagan,
``From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales'',
The American Naturalist, 183:5, E154-E167 (2014) \doi{10.1086/675504}.

C. H. Fleming et al, ``A comprehensive framework for handling location error in animal tracking data'',
bioRxiv (2020) \doi{10.1101/2020.06.12.130195}.
}

\author{ C. H. Fleming and J. M. Calabrese. }
\note{ Prior to \code{ctmm} v0.3.6, \code{fast=FALSE} used the lag-weighted esitmator of Fleming et al (2014). Lag weights have been abandoned in favor of interval weights, which are less sensitive to sampling irregularity. The same weighting formulas are used, but with \code{dt} instead of the current lag. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{vignette("variogram")}, \code{\link{correlogram}}, \code{\link{mean.variogram}}, \code{\link{plot.variogram}}, \code{\link{variogram.fit}}. }

\examples{
#Load package and data
library(ctmm)
data(buffalo)

#Extract movement data for a single animal
DATA <- buffalo$Cilla

#Calculate variogram
SVF <- variogram(DATA)

#Plot the variogram with 50\% and 95\% CIs
plot(SVF,level=c(0.5,0.95))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

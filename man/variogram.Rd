\name{variogram}
\alias{variogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate an empirical variogram from movement data
}
\description{
This function calculates the isotropic empirical variogram of two-dimensional tracking data for visualizing stationary (time-averaged) autocorrelation structure.
One of two algorithms is used. The slow \eqn{O(n^2)} algorithm is based upon the lag-weighted variogram of Fleming et al (2014).
Additional modifications have also been included to accommodate any drift in the sampling rate.
The fast \eqn{O(n \log n)} algorithm is based upon the FFT method of Marcotte (1996), with some tweaks to better handle irregularly sampled data.
Both methods reduce to the standard \dQuote{method of moments} estimator in the case of evenly \emph{scheduled} data, even with missing observations, but they produce slightly different outputs for irregularly sampled data.
}
\usage{
variogram(data,dt=NULL,fast=TRUE,CI="Markov",axes=c("x","y"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ \code{telemetry} data object of the 2D timeseries data. }
  \item{dt}{ Lag bin width. An ordered array will yield a progressive coarsening of the lags. Defaults to the median sampling interval.}
  \item{fast}{ Use the lag-weighted algorithm if \code{FALSE} or the FFT algorithm if \code{TRUE}. The slow algorithm outputs a progress bar.}
  \item{CI}{Defaults to only consider non-overlapping lags independent. All unique lags are considered independent with "IID".}
  \item{axes}{Array of axes to calculate an average (isotropic) variogram for.}
}
\details{
If no \code{dt} is specified, the median sampling interval is used. This is typically a good assumption for most data, even when there are gaps. A \code{dt} coarser than the sampling interval may bias the variogram (particuarly if \code{fast=TRUE}) and so this should be reserved for poor data quality.

For irregularly sampled data, it may be useful to provide an array of time-lag bin widths to progressively coarsen the variogram. I.e., if you made the very bad choice of changing your sampling interval on the fly from \code{dt1} to \code{dt2}, where \code{dt1} \eqn{<} \code{dt2}, the an appropriate choice would be \code{dt=c(dt1,dt2)}. On the other hand, if your sampling is itself a noisy process, then you might want to introduce larger and larger \code{dt} components as the visual appearance of the variogram breaks down with increasing lags.
Alternatively, you might try the \code{fast=FALSE} option or aggregating multiple individuals with \code{\link{mean.variogram}}.

In standard variogram regression treatments, all lags are considered as independent (\code{CI="IID"}) for the purposes of confidence-interval estimation, even if they overlap in time. However, in high resolution datasets this will produce vastly underestimated confidence intervals. Therefore, the default \code{CI="Markov"} behavior is to consider only the maximum number of non-overlapping lags in calculating confidence intervals.
}
\value{
Returns a variogram object (class variogram) which is a dataframe containing the time-lag, \code{lag}, the semi-variance estimate at that lag, \code{SVF}, and the approximate number of degrees of freedom associated with that semi-variance, \code{DOF}, with which its confidence intervals can be estimated.
}
\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, and W. F. Fagan. (2014). From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales. \href{http://www.jstor.org/discover/10.1086/675504}{The American Naturalist, 183(5), E154-E167}.

D. Marcotte. (1996). Fast variogram computation with FFT. \emph{Computers and Geosciences} 22(10), 1175-1186.
}

\author{ C. H. Fleming and J. M. Calabrese. }
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{vignette("variogram")}, \code{\link{mean.variogram}}, \code{\link{plot.variogram}}, \code{\link{variogram.fit}}. }

\examples{
#Load package and data
library(ctmm)
data(buffalo)

#Extract movement data for a single animal
cilla <- buffalo[[1]]

#Calculate variogram
SVF <- variogram(cilla)

#Plot the variogram
plot(SVF)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{summary.ctmm}
\alias{summary.ctmm}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize a continuous-time movement model}
\description{ This function returns a list of biologically interesting parameters in human readable format, as derived from a continuous-time movement model.}
\usage{
\S3method{summary}{ctmm}(object,level=0.95,level.UD=0.95,units=TRUE,IC=NULL,MSPE=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{ctmm} movement-model object from the output of \code{ctmm.fit}. }
  \item{level}{ Confidence level for parameter estimates. }
  \item{level.UD}{ Coverage level for the Gaussian home-range area. }
  \item{units}{Convert result to natural units.}
  \item{IC}{Information criteria for sorting lists of \code{ctmm} objects. Can be \code{"AICc"}, \code{"AIC"}, \code{"BIC"}, \code{"LOOCV"}, \code{"HSCV"}, or none (\code{NA}). AICc is approximate.}
  \item{MSPE}{Sort models with the same autocovariance structure by the mean square predictive error of \code{"position"}, \code{"velocity"}, or not (\code{NA}).}
  \item{...}{Unused options.}
}

%\details{}

\value{
If summary is called with a single \code{ctmm} object output from \code{\link{ctmm.fit}}, then a list is returned with the effective sample sizes of various parameter estimates (\code{DOF}) and a parameter estimate table \code{CI}, with low, point, and high estimates for the following possible parameters:
\describe{
\item{\code{tau}}{The autocorrelation timescales. \code{tau position} is also the home-range crossing timescale.}
\item{\code{area}}{The Gaussian home-range area, where the point estimate has a significance level of \code{level.UD}. I.e., the core home range is where the animal is located 50\% of the time with \code{level.UD=0.50}.
This point estimate itself is subject to uncertainty, and is given confidence intervals derived from \code{level}.
This Gaussian estimate differs from the kernel density estimate of \code{\link{summary.UD}}. The Gaussian estimate has more statistical efficiency, but is less related to space use for non-Gaussian processes.}
\item{\code{speed}}{The Gaussian root-mean-square (RMS) velocity, which is a convenient measure of average speed but not the conventional measure of average speed (see \code{\link{speed}}).}
}
Furthermore, if summary is called on a population-level model, then population-level standard deviations (\code{SD}) and coefficients of variation (\code{CoV}) are also returned.

If summary is called on a list of \code{ctmm} objects output from \code{\link{ctmm.select}}, then a table is returned with the model names and IC differences for comparison across autocovariance structures. The mean square prediction error (MSPE) is also returned for comparison across trend structures (with autocovariance structure fixed). For the model names, "IID" denotes the uncorrelated bi-variate Gaussian model, "OU" denotes the continuous-position Ornstein-Uhlenbeck model, "OUF" denotes the continuous-velocity Ornstein-Uhlenbeck-F model, "OUf" denotes the OUF model where the two autocorrelation timescales cannot be statistically distinguished.
}

%\references{}

\author{ C. H. Fleming. }

\note{ Confidence intervals on the autocorrelation timescales assume they are sufficiently greater than zero and less than infinity.

\code{IC="LOOCV"} can only be attempted if also specified during \code{\link{ctmm.select}}, as this argument requires additional calculations.

Prior to \code{ctmm} v0.6.2, timescale confidence intervals were constructed from normal and inverse-normal sampling distributions, whereas v0.6.2 onward uses gamma and inverse-gamma sampling distributions.

In \code{ctmm} v0.5.1 onward the MSPE is averaged over all possible times instead of over all sampled times.

In \code{ctmm} v0.3.4 the speed estimate was fixed to be the RMS velocity and not \eqn{1/\sqrt{2}} times the RMS velocity.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}}, \code{\link{ctmm.select}}. }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# Extract movement data for a single animal
DATA <- buffalo$Cilla

# fit model
GUESS <- ctmm.guess(DATA,interactive=FALSE)
FIT <- ctmm.fit(DATA,GUESS)

# Tell us something interpretable
summary(FIT)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

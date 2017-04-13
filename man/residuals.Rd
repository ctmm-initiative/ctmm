\name{residuals.ctmm}
\alias{residuals.ctmm}
\alias{residuals.telemetry}
\alias{correlogram}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate model fit residuals and assess their autocorrelation}
\description{
These functions calculate the residuals of a CTMM model, which should be standardized and IID if the model correctly specified.
A correlogram method is also provided to assess remaining autocorrelation.
This function is analogous to \code{acf}, but can handle missing data.
}
\usage{
\S3method{residuals}{ctmm}(object,data,...)

\S3method{residuals}{telemetry}(object,CTMM,...)

correlogram(data,dt=NULL,axes=c("x","y"))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{ctmm} model object or \code{telemetry} data object for calculating residuals.}
  \item{data}{\code{telemetry} data object or \code{data.frame} with time column \code{t} and data columns \code{axes}.}
  \item{CTMM}{\code{ctmm} model object.}
  \item{...}{Unused arguments.}
  \item{dt}{ Lag bin width. An ordered array will yield a progressive coarsening of the lags. Defaults to the median sampling interval.}
  \item{axes}{Array of axes to calculate an average (isotropic) correlogram for.}
}
\details{
Given a \code{telemetry} dataset and \code{ctmm} model, \code{residuals} calculates the standardized residuals of the Kalman filter, which can be tested for independence. The residuals object can then be plotted with \code{plot} or fed into the \code{correlogram} method to test independence. Output of the correlogram can then be plotted as well, though \code{zoom} is much more useful. This function is more-or-less equivalent to \code{acf}, however it is coded to handle multiple dimensions and missing data. Another difference from \code{acf} is that the confidence intervals here are calculated per lag, using a Fisher transformation of the correlation, which is more accurate far from lag zero.
}

\value{
\code{residuals} return a residual object (class \code{telemetry}, but flagged as residual) and \code{correlogram} returns a correlogram object (class \code{variogram}, but flagged as an ACF).
}
%\references{}

\author{ C. H. Fleming}

\note{If the sampling schedule is irregular, permitting gaps, then the correlogram may not look good even if the model is correctly specified. In this case the correlogram of the residuals should be compared to the correlogram of simulated residuals, using "data" simulated from the fit model and with the same sampling schedule.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ \code{\link{plot.periodogram}}. }

%\examples{}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

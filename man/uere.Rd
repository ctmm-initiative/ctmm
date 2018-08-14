\name{uere}
\alias{uere}
\alias{uere<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate UERE from calibration data}

\description{Functions for estimating and assigning the User Equivalent Range Error (UERE) of a GPS device from calibration data.}

\usage{
uere(data,override=FALSE,precision=1/2,trace=FALSE)

uere(data) <- value
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Telemetry object or list of telemetry objects, preferably with DOP columns.}
  \item{override}{Return a UERE estimate even for data that are calibrated.}
  \item{precision}{Fraction of maximum possible digits of precision to target in categorical error fitting. \code{precision=1/2} results in about 7 decimal digits of precision.}
  %\item{integrate}{Consider HDOP and VDOP as being specified by the same UERE value.}
  \item{trace}{Print extra details, such as the UERE's 95\% confidence interval.}
  \item{value}{UERE value(s) to assign to telemetry data (see details).}
}

\details{Often times GPS animal tracking devices return HDOP values but do not specifiy the device's UERE necessary to transform the HDOP values into absolute errors.
\code{uere()} allows users to estimate the UERE from calibration data, where the device was left fixed over a period of time.
Otherwise, when \code{error=TRUE} in \code{\link{ctmm}}, \code{\link{ctmm.fit}} will estimate the UERE simultaneously with the movement model, which is less reliable than using calibration data.

The calibration UERE can then be applied to tracking data during import, with \code{\link{as.telemetry}}, or after import with the \code{uere()<-} assignment method.

Ideally, device HDOP and VDOP values are scaled to the same UERE.
However, in many cases not all DOP values are supplied and some must be approximated.
In such cases, different UERE values need to be estimated for each set of DOP values, via the \code{axes} argument of \code{uere}.
A \code{c} combination of properly named UERE values can then be used in assignment via \code{uere<-}.
Proper names include \code{"horizontal"}, \code{"vertical"}, and \code{"speed"}.
}

\value{The UERE estimate.}

%\references{}

\author{C. H. Fleming}

\note{The GPS device should be fixed during calibraiton.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{as.telemetry}}, \code{\link{residuals.telemetry}}.  }

\examples{
# Load package and data
library(ctmm)
data(turtle)

# the first two datasets are calibration data
names(turtle)

# estimate UERE from calibration data
UERE <- uere(turtle[1:2])

# apply UERE to entire dataset
uere(turtle) <- UERE

# calculate residuals of calibration data
RES <- lapply(turtle[1:2],residuals)

# scatter plot of residuals with 50\% and 95\% quantiles
plot(RES,col.DF=NA,level.UD=c(0.50,0.95))

# check calibration data for autocorrelation using fast=FALSE because samples are small
ACFS <- lapply(RES,function(R){correlogram(R,fast=FALSE,dt=10 \%#\% 'min')})

# pooling ACFs
ACF <- mean(ACFS)

plot(ACF)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{internal}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

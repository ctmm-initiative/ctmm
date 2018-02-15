\name{uere}
\alias{uere}
\alias{uere<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate UERE from calibration data}

\description{Functions for estimating and assigning the User Equivalent Range Error (UERE) of a GPS device from calibration data.}

\usage{
uere(data,override=FALSE,integrate=FALSE,diagnostic=FALSE)

uere(data) <- value
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Telemetry object or list of telemetry objects, preferably with DOP columns.}
  \item{override}{Return a UERE estimate even for data that are calibrated.}
  \item{integrate}{Consider HDOP and VDOP as being specified by the same UERE value.}
  \item{diagnostic}{Return extra details, such as the UERE's 95\% confidence interval and a plot.}
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

\seealso{ \code{\link{as.telemetry}}.  }

%\examples{\donttest{
%# Load package and data
%library(ctmm)
%data(buffalo)

%# Extract movement data for a single animal
%billy <- buffalo[[1]]

%# Plot the data
%plot(billy)
%}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{internal}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

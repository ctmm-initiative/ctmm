\name{uere}
\alias{uere}
\alias{uere<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate UERE from calibration data}

\description{Functions for estimating and assigning the User Equivalent Range Error (UERE) of a GPS device from calibration data.}

\usage{
uere(data,diagnostic=FALSE)

uere(data) <- value
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Telemetry object or list of telemetry objects with Movebank column \code{GPS.HDOP}.}
  \item{diagnostic}{Return extra details, such as the UERE's 95\% confidence interval and a plot.}
  \item{value}{UERE value to assign to telemetry data.}
}

\details{Often times GPS animal tracking devices return HDOP values but do not specifiy the device's UERE necessary to transform the HDOP values into absolute errors.
\code{uere()} allows users to estimate the UERE from calibration data, where the device was left fixed over a period of time.
Otherwise, when \code{error=TRUE} in \code{\link{ctmm}}, \code{\link{ctmm.fit}} will estimate the UERE simultaneously with the movement model, which is less reliable than using calibration data.

The calibration UERE can then be applied to tracking data during import, with \code{\link{as.telemetry}}, or after import with the \code{uere()<-} assignment method.}

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

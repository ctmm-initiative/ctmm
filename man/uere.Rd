\name{uere}
\alias{uere}
\alias{uere.fit}
\alias{uere<-}
\alias{summary.UERE}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate UERE from calibration data}

\description{Functions for estimating and assigning the User Equivalent Range Error (UERE) of a GPS device from calibration data.}

\usage{
uere(data)

uere(data) <- value

uere.fit(data,precision=1/2)

\S3method{summary}{UERE}(object,level=0.95,...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{\code{telemetry} object or list of telemetry objects, preferably with DOP columns.}
  \item{value}{UERE value(s) to assign to telemetry data (see details).}
  \item{precision}{Fraction of maximum possible digits of precision to target in categorical error fitting. \code{precision=1/2} results in about 7 decimal digits of precision.}
  \item{object}{\code{UERE} object to summarize.}
  \item{level}{Confidence level for UERE estimate confidence intervals.}
  \item{...}{Further arguments are ignored.}
}

\details{Often times GPS animal tracking devices return HDOP values but do not specifiy the device's UERE necessary to transform the HDOP values into absolute errors.
\code{uere.fit()} allows users to estimate the UERE from calibration data, where the device was left fixed over a period of time.
The calibration UERE can then be applied to tracking data with the \code{uere()<-} assignment method.
Otherwise, when \code{error=TRUE} in \code{\link{ctmm}}, \code{\link{ctmm.fit}} will estimate the UERE simultaneously with the movement model, which is less reliable than using calibration data.

\code{summary()} applied to single \code{UERE} object will return UERE parameter estimates and confidence intervals, while \code{summary()} applied to a list of \code{UERE} objects will return a model-selection table, with AICc and reduced Z squared (goodness of fit) values.
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
UERE <- uere.fit(turtle[1:2])
# inspect UERE estimate
summary(UERE)

# assign UERE to entire dataset
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

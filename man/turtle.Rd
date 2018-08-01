\name{turtle}
\alias{turtle}
\docType{data}

\title{Wood turtle GPS and calibration dataset from Working Land and Seascapes.}

\description{\code{x-y} projected GPS data from 2 calibration runs and 2 wood turtles.
Please contact Tom Akre (akret@si.edu) if you want to publish with these data.}

\usage{data("turtle")}

\format{ A list of 4 \code{telemetry} objects.}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

%\references{}

\seealso{ \code{\link{as.telemetry}}, \code{\link{plot.telemetry}}, \code{\link{uere}}, \code{\link{buffalo}}, \code{\link{coati}}, \code{\link{gazelle}}, \code{\link{wolf}}. }

\examples{
# Load package and data
library(ctmm)
data("turtle")

# Plot a turtle's locations
plot(turtle[[3]])
}

\keyword{datasets}

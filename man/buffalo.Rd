\name{buffalo}
\alias{buffalo}
\encoding{UTF-8}
\docType{data}

\title{African buffalo GPS dataset from Kruger National Park, South Africa.}

\description{GPS data on six African buffalo. When using this dataset, please cite the original article by Getz et al (2007) and the Movebank data package (Cross et al, 2016).}

\usage{data("buffalo")}

\format{ A list of 6 \code{telemetry} objects.}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
W. M. Getz, S. Fortmann-Roe, P. C. Cross, A. J. Lyons, S. J. Ryan, C. C. Wilmers.
LoCoH: Nonparameteric kernel methods for constructing home ranges and utilization distributions.
PLoS ONE 2:2, e207 (2007).

P. C. Cross, J. A. Bowers, C. T. Hay, J. Wolhuter, P. Buss, M. Hofmeyr, J. T. du Toit, W. M. Getz.
Data from: Nonparameteric kernel methods for constructing home ranges and utilization distributions.
Movebank Data Repository. DOI:10.5441/001/1.j900f88t (2016).
}

\note{In \code{ctmm} v0.3.2 the erroneous location fix \code{606} was removed from \code{buffalo[[4]]} "Pepper".}

\seealso{ \code{\link{as.telemetry}}, \code{\link{plot.telemetry}}, \code{\link{coati}}, \code{\link{gazelle}}, \code{\link{jaguar}}, \code{\link{pelican}}, \code{\link{turtle}}, \code{\link{wolf}}. }

\examples{
# Load package and data
library(ctmm)
data("buffalo")

# Extract movement data for a single animal
Cilla <- buffalo$Cilla

# Plot all sampled locations
plot(Cilla)
}

\keyword{datasets}

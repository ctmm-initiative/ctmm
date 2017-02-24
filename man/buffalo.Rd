\name{buffalo}
\alias{buffalo}
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
Getz W. M., Fortmann-Roe S., Cross P. C., Lyons A. J., Ryan S. J., Wilmers C. C. (2007).
LoCoH: Nonparameteric kernel methods for constructing home ranges and utilization distributions.
PLoS ONE 2(2): e207. DOI:10.1371/journal.pone.0000207

Cross P. C., Bowers J. A., Hay C. T., Wolhuter J., Buss P., Hofmeyr M., du Toit J. T., Getz W. M. (2016).
Data from: Nonparameteric kernel methods for constructing home ranges and utilization distributions.
Movebank Data Repository. DOI:10.5441/001/1.j900f88t
}

\note{In \code{ctmm} v0.3.2 the erroneous location fix \code{606} was removed from \code{buffalo[[4]]} "Pepper".}

\seealso{ \code{\link{as.telemetry}}, \code{\link{gazelle}}, \code{\link{plot.telemetry}}. }

\examples{
# Load package and data
library(ctmm)
data("buffalo")

# Extract movement data for a single animal
cilla <- buffalo[[1]]

# Plot all sampled locations
plot(cilla)
}

\keyword{datasets}
